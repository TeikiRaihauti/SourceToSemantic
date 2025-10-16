def _nint(x):
    import math
    if x >= 0:
        return int(math.floor(x + 0.5))
    else:
        return -int(math.floor(-x + 0.5))


def SOILT_EPIC(B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,
               TAVG, TMAX, TMIN, WetDay, WFT, WW,
               TMA, ST, X2_PREV):
    """
    EPIC soil temperature by layer.

    Parameters
    - B: float
    - BCV: float
    - CUMDPT: float (mm)
    - DP: float (mm)
    - DSMID: sequence length NLAYR (mm)
    - NLAYR: int
    - PESW: float (cm)
    - TAV: float
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - WetDay: int (0 or 1)
    - WFT: float
    - WW: float
    - TMA: sequence length 5
    - ST: sequence length NLAYR
    - X2_PREV: float

    Returns (TMA_new, SRFTEMP, ST_new, X2_AVG, X2_PREV_new)
    """
    import math

    # Ensure local copies (purity)
    TMA_new = [float(v) for v in TMA]
    ST_new = [float(v) for v in ST]

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Surface forcing temperature X2
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory (replicate Fortran order)
    TMA_new[0] = X2
    for k in range(4, 0, -1):
        TMA_new[k] = TMA_new[k - 1]
    X2_AVG = sum(TMA_new) / 5.0

    # Cover effect and surface temperature
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)

    X1 = TAV - X3
    LAG = 0.5

    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST_new[L] = LAG * ST_new[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_new = X2_AVG

    return TMA_new, SRFTEMP, ST_new, X2_AVG, X2_PREV_new


def STEMP_EPIC_Initialize(ISWWAT,
                          BD, DLAYR, DS, DUL, LL, NLAYR,
                          SW,
                          TAVG, TMAX, TMIN, TAV,
                          ORGC_MULCHMASS, WATER_SNOW):
    """
    Initialization for EPIC soil temperature routine (SEASINIT block).

    Inputs:
    - ISWWAT: 'Y' or 'N'
    - BD: list length NLAYR (g/cm3)
    - DLAYR: list length NLAYR (cm)
    - DS: list length NLAYR (cm), cumulative depth by layer
    - DUL: list length NLAYR (cm3/cm3)
    - LL: list length NLAYR (cm3/cm3)
    - NLAYR: int
    - SW: list length NLAYR (cm3/cm3) initial soil water
    - TAVG, TMAX, TMIN, TAV: float
    - ORGC_MULCHMASS: float (kg/ha)
    - WATER_SNOW: float (mm)

    Returns:
    - CUMDPT: float (mm)
    - DSMID: list length NLAYR (mm)
    - TDL: float (cm)
    - TMA: list length 5
    - NDays: int
    - WetDay: list length 30 of ints
    - X2_PREV: float
    - SRFTEMP: float
    - ST: list length NLAYR
    """
    import math

    # Copy inputs
    BD = [float(v) for v in BD[:NLAYR]]
    DLAYR = [float(v) for v in DLAYR[:NLAYR]]
    DS = [float(v) for v in DS[:NLAYR]]
    DUL = [float(v) for v in DUL[:NLAYR]]
    LL = [float(v) for v in LL[:NLAYR]]
    SWI = [float(v) for v in SW[:NLAYR]]

    # Depth calculations
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to midpoint (cm*10*0.5)
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    TMA = []
    for _ in range(5):
        TMA.append(_nint(TAVG * 10000.0) / 10000.0)

    ST = [float(TAVG)] * NLAYR

    # 30-day wet day memory
    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    # Soil cover effect
    CV = (ORGC_MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = WATER_SNOW / (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) if (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Spin-up
    SRFTEMP = float(TAVG)
    X2_PREV = 0.0
    for _ in range(8):
        TMA, SRFTEMP, ST, _, X2_PREV = SOILT_EPIC(
            B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,
            TAVG, TMAX, TMIN, 0, WFT, WW,
            TMA, ST, X2_PREV
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(ISWWAT,
               BD, DLAYR, DS, DUL, LL, NLAYR,
               SW,
               TAVG, TMAX, TMIN, TAV,
               WEATHER_RAIN, MGMT_DEPIR, PLANT_BIOMAS, ORGC_MULCHMASS, WATER_SNOW,
               CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST):
    """
    Daily EPIC soil temperature rate calculations (RATE block).

    Inputs:
    - ISWWAT: 'Y' or 'N'
    - BD, DLAYR, DS, DUL, LL: lists length NLAYR
    - NLAYR: int
    - SW: list length NLAYR
    - TAVG, TMAX, TMIN, TAV: float
    - WEATHER_RAIN: float (mm)
    - MGMT_DEPIR: float (mm)
    - PLANT_BIOMAS: float (kg/ha)
    - ORGC_MULCHMASS: float (kg/ha)
    - WATER_SNOW: float (mm)
    - CUMDPT: float (mm)
    - DSMID: list length NLAYR (mm)
    - TDL: float (cm) state variable (will be updated per original code)
    - TMA: list length 5
    - NDays: int
    - WetDay: list length 30
    - X2_PREV: float
    - SRFTEMP: float
    - ST: list length NLAYR

    Returns updated tuple:
    (TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST)
    """
    import math

    # Copy arrays to avoid side effects
    BD = [float(v) for v in BD[:NLAYR]]
    DLAYR = [float(v) for v in DLAYR[:NLAYR]]
    DS = [float(v) for v in DS[:NLAYR]]
    DUL = [float(v) for v in DUL[:NLAYR]]
    LL = [float(v) for v in LL[:NLAYR]]
    SW = [float(v) for v in SW[:NLAYR]]
    DSMID = [float(v) for v in DSMID[:NLAYR]]
    TMA = [float(v) for v in TMA[:5]]
    ST = [float(v) for v in ST[:NLAYR]]
    WetDay = [int(v) for v in WetDay[:30]] + [0] * max(0, 30 - len(WetDay))

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # Note: intentionally not reset per original code
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # 30-day memory of wet days
    if NDays == 30:
        # shift positions 1..29 to 0..28
        WetDay = WetDay[1:30] + [0]
        wet_index = 29
    else:
        NDays = NDays + 1
        wet_index = NDays - 1

    if WEATHER_RAIN + MGMT_DEPIR > 1.0e-6:
        WetDay[wet_index] = 1
    else:
        WetDay[wet_index] = 0

    NWetDays = sum(WetDay)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function
    CV = (PLANT_BIOMAS + ORGC_MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = WATER_SNOW / (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) if (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Update soil temperatures
    TMA, SRFTEMP, ST, _, X2_PREV = SOILT_EPIC(
        B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,
        TAVG, TMAX, TMIN, WetDay[wet_index], WFT, WW,
        TMA, ST, X2_PREV
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def test_STEMP_EPIC_basic():
    """
    Basic test derived from ASKEE_DSSAT_EPIC.for
    Does not execute at import; call explicitly in a test harness.
    """
    # Inputs similar to ASKEE program
    ISWWAT = 'Y'
    NLAYR = 4
    BD = [1.6] * NLAYR
    DLAYR = [10.0] * NLAYR
    DS = [10.0, 20.0, 30.0, 40.0]
    DUL = [0.3] * NLAYR
    LL = [0.2] * NLAYR

    SW = [0.2] * NLAYR
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0

    WEATHER_RAIN = 0.0
    MGMT_DEPIR = 0.0
    PLANT_BIOMAS = 0.0
    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0

    # Initialize
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC_Initialize(
        ISWWAT,
        BD, DLAYR, DS, DUL, LL, NLAYR,
        SW,
        TAVG, TMAX, TMIN, TAV,
        ORGC_MULCHMASS, WATER_SNOW
    )

    # Simple assertions on initialization
    assert len(DSMID) == NLAYR
    assert len(TMA) == 5
    assert len(ST) == NLAYR
    assert len(WetDay) == 30

    # Run a few daily steps
    for _ in range(5):
        TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC(
            ISWWAT,
            BD, DLAYR, DS, DUL, LL, NLAYR,
            SW,
            TAVG, TMAX, TMIN, TAV,
            WEATHER_RAIN, MGMT_DEPIR, PLANT_BIOMAS, ORGC_MULCHMASS, WATER_SNOW,
            CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
        )

    # Post-step assertions
    assert NDays == 5
    assert isinstance(SRFTEMP, float)
    assert all(isinstance(x, float) for x in ST)
    assert all(isinstance(x, float) for x in TMA)