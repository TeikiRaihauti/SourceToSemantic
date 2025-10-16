import math


def SOILT_EPIC(B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV, TAVG, TMAX, TMIN,
               WetDay, WFT, WW, TMA, ST, X2_PREV):
    """
    EPIC soil temperature routine that updates soil temperature by layer.

    Inputs:
    - B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV, TAVG, TMAX, TMIN, WetDay, WFT, WW
    InOut:
    - TMA (list of 5), ST (list length NLAYR), X2_PREV

    Returns:
    - TMA (updated)
    - SRFTEMP
    - ST (updated)
    - X2_AVG
    - X2_PREV (updated)
    """
    # Ensure list copies to avoid in-place mutation of caller inputs
    TMA = list(TMA)
    ST = list(ST)

    # Water content ratio (dimensionless)
    # WC = PESW(cm) / (WW(dimensionless) * CUMDPT(mm)) * (10 mm/cm)
    # Guard against division by zero and negative WW
    denom = WW * CUMDPT if WW * CUMDPT != 0.0 else 1e-12
    WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Surface temperature driver X2 (Potter & Williams, 1994)
    if WetDay > 0:
        # Eqn. 2: X2 = WFT*(TAVG - TMIN) + TMIN
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        # Eqn. 1: X2 = WFT*(TMAX - TAVG) + TAVG + 2.
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update TMA (note: preserves original Fortran order of operations)
    old_TMA = TMA[:]
    TMA[0] = X2
    for K in range(4, 0, -1):  # 5->2 in Fortran 1-based corresponds to 4->1 zero-based
        TMA[K] = old_TMA[K - 1]
    X2_AVG = sum(TMA) / 5.0

    # Eqn. 4: X3 = (1 - BCV)*X2_AVG + BCV*X2_PREV
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    # SRFTEMP = min(X2_AVG, X3)
    SRFTEMP = min(X2_AVG, X3)

    # Eqn. 6 (partial): X1 = TAV - X3
    X1 = TAV - X3

    # Layer temperatures
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0  # Eqn. 8
        # Eqn. 7: F attenuation
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        # Eqn. 6 full: T(L) = LAG*T(L) + (1-LAG)*(F*X1 + X3)
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_AVG, X2_PREV


def STEMP_EPIC(DYNAMIC,
               ISWWAT,
               BD, DLAYR, DS, DUL, LL, NLAYR,
               SW,
               TAVG, TMAX, TMIN, TAV,
               RAIN, DEPIR, BIOMAS, MULCHMASS, SNOW,
               CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST):
    """
    Determines soil temperature by layer using the EPIC method.

    This function combines the seasonal initialization (SEASINIT) and the daily
    rate computation (RATE) from the original Fortran subroutine. Use DYNAMIC to
    select the branch:
      - DYNAMIC == 2 or 'SEASINIT' -> seasonal initialization
      - DYNAMIC == 3 or 'RATE'     -> daily rate computation

    Inputs:
    - DYNAMIC: int or str ('SEASINIT' or 'RATE')
    - ISWWAT: 'Y' or 'N'
    - BD, DLAYR, DS, DUL, LL: lists of length NLAYR
    - NLAYR: number of soil layers
    - SW: list of length NLAYR
    - TAVG, TMAX, TMIN, TAV: floats (degC)
    - RAIN: precipitation (mm)
    - DEPIR: irrigation depth (mm)
    - BIOMAS: plant biomass (kg/ha)
    - MULCHMASS: mulch mass (kg/ha)
    - SNOW: snow cover (mm)

    InOut state variables (provide current values; updated values returned):
    - CUMDPT: cumulative profile depth (mm)
    - DSMID: list length NLAYR; depth to midpoint of each layer (mm)
    - TDL: total water at DUL (cm) - preserved behavior from Fortran
    - TMA: list length 5; prior 5 days of surface temp driver
    - NDays: number of days tracked in 30-day wet-day memory
    - WetDay: list length 30; binary wet-day history
    - X2_PREV: previous day's X2 average (degC)
    - SRFTEMP: surface temperature (degC)
    - ST: list length NLAYR; soil temperature by layer (degC)

    Returns:
    - CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    """
    # Copy arrays to avoid mutating caller-owned lists
    BD = list(BD)
    DLAYR = list(DLAYR)
    DS = list(DS)
    DUL = list(DUL)
    LL = list(LL)
    SW = list(SW)
    DSMID = list(DSMID)
    TMA = list(TMA)
    WetDay = list(WetDay)
    ST = list(ST)

    def fortran_nint(x):
        return math.floor(x + 0.5) if x >= 0.0 else -math.floor(-x + 0.5)

    # Determine dynamic branch
    dynamic_is_seasinit = (DYNAMIC == 2) or (isinstance(DYNAMIC, str) and DYNAMIC.upper() == 'SEASINIT')
    dynamic_is_rate = (DYNAMIC == 3) or (isinstance(DYNAMIC, str) and DYNAMIC.upper() == 'RATE')

    if dynamic_is_seasinit:
        # Initialization
        SWI = SW[:]  # initial soil water equal to provided SW
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL = 0.0
        CUMDPT = 0.0
        for L in range(NLAYR):
            DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm depth to midpoint of layer
            CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm profile depth
            TBD = TBD + BD[L] * DLAYR[L]
            TLL = TLL + LL[L] * DLAYR[L]
            TSW = TSW + SWI[L] * DLAYR[L]
            TDL = TDL + DUL[L] * DLAYR[L]

        if ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)  # cm
        else:
            PESW = max(0.0, TDL - TLL)  # cm

        # Average bulk density (g/cm3)
        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX              # mm
        WW = 0.356 - 0.144 * ABD               # volumetric fraction
        B = math.log(500.0 / DP)

        # Initialize TMA with TAVG (rounded to 4 decimals as in Fortran NINT)
        for I in range(5):
            TMA[I] = fortran_nint(TAVG * 10000.0) / 10000.0
        # Initialize ST for all layers
        for L in range(NLAYR):
            ST[L] = TAVG

        # Wet-day memory
        WFT = 0.1
        WetDay = [0] * 30
        NDays = 0

        # Soil cover function (mulch + snow)
        CV = (MULCHMASS) / 1000.0  # t/ha
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
        BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
        BCV = max(BCV1, BCV2)

        # Spin-up calls
        for _ in range(8):
            TMA, SRFTEMP, ST, _, X2_PREV = SOILT_EPIC(
                B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,
                TAVG, TMAX, TMIN, 0, WFT, WW, TMA, ST, X2_PREV
            )

    elif dynamic_is_rate:
        # Daily rate calculations
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        for L in range(NLAYR):
            TBD = TBD + BD[L] * DLAYR[L]
            TDL = TDL + DUL[L] * DLAYR[L]  # Preserved from Fortran (no reset)
            TLL = TLL + LL[L] * DLAYR[L]
            TSW = TSW + SW[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX   # mm
        WW = 0.356 - 0.144 * ABD    # volumetric fraction
        B = math.log(500.0 / DP)

        if ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)  # cm
        else:
            PESW = max(0.0, TDL - TLL)  # cm

        # Update 30-day wet-day memory with rainfall + irrigation
        if NDays == 30:
            # Shift left
            for i in range(29):
                WetDay[i] = WetDay[i + 1]
            # Keep NDays at 30
        else:
            NDays = NDays + 1

        if (RAIN + DEPIR) > 1.0e-6:
            WetDay[NDays - 1] = 1
        else:
            WetDay[NDays - 1] = 0

        NWetDays = sum(WetDay[:NDays]) if NDays > 0 else 0
        WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

        # Soil cover function (plant biomass + mulch + snow)
        CV = (BIOMAS + MULCHMASS) / 1000.0  # t/ha
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
        BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
        BCV = max(BCV1, BCV2)

        # Update temperatures
        TMA, SRFTEMP, ST, _, X2_PREV = SOILT_EPIC(
            B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,
            TAVG, TMAX, TMIN, WetDay[NDays - 1], WFT, WW, TMA, ST, X2_PREV
        )

    # Return updated state
    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def test_STEMP_EPIC_from_ASKEE():
    """
    Basic sanity test derived from ASKEE_DSSAT_EPIC.for driver program.
    Does not assert exact numerical values but checks shape and execution.
    """
    # Setup based on ASKEE
    ISWWAT = 'Y'
    NLAYR = 4
    BD = [1.6] * NLAYR
    DLAYR = [10.0] * NLAYR  # cm
    DS = [10.0, 20.0, 30.0, 40.0]  # cm cumulative
    DUL = [0.3] * NLAYR
    LL = [0.2] * NLAYR
    SW = [0.2] * NLAYR

    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0

    RAIN = 0.0
    DEPIR = 0.0
    BIOMAS = 0.0
    MULCHMASS = 0.0
    SNOW = 0.0

    # Initial states
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TDL = 0.0
    TMA = [0.0] * 5
    NDays = 0
    WetDay = [0] * 30
    X2_PREV = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NLAYR

    # Seasonal initialization
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC(
        'SEASINIT',
        ISWWAT,
        BD, DLAYR, DS, DUL, LL, NLAYR,
        SW,
        TAVG, TMAX, TMIN, TAV,
        RAIN, DEPIR, BIOMAS, MULCHMASS, SNOW,
        CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    )

    # One daily rate step
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC(
        'RATE',
        ISWWAT,
        BD, DLAYR, DS, DUL, LL, NLAYR,
        SW,
        TAVG, TMAX, TMIN, TAV,
        RAIN, DEPIR, BIOMAS, MULCHMASS, SNOW,
        CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    )

    # Sanity checks
    assert isinstance(CUMDPT, float)
    assert isinstance(TDL, float)
    assert isinstance(SRFTEMP, float)
    assert isinstance(NDays, int)
    assert len(DSMID) == NLAYR
    assert len(ST) == NLAYR
    assert len(TMA) == 5
    assert len(WetDay) == 30