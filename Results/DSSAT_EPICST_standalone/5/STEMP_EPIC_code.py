def STEMP_EPIC_initialize(
    iswitch_ISWWAT: str,
    soilprop_BD: list[float],
    soilprop_DLAYR: list[float],
    soilprop_DS: list[float],
    soilprop_DUL: list[float],
    soilprop_LL: list[float],
    soilprop_NLAYR: int,
    SW: list[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    ORGC_MULCHMASS: float,
    WATER_SNOW: float,
    X2_PREV: float = 0.0,
) -> tuple[
    float,           # CUMDPT
    list[float],     # DSMID
    float,           # TDL
    list[float],     # TMA (length 5)
    int,             # NDays
    list[int],       # WetDay (length 30)
    float,           # X2_PREV
    float,           # SRFTEMP
    list[float],     # ST (length NLAYR)
]:
    """
    Initialize EPIC soil temperature state.

    Inputs:
    - iswitch_ISWWAT: str ('Y' if water balance simulated, else other)
    - soilprop_BD: list[float], soil bulk density by layer (g/cm3)
    - soilprop_DLAYR: list[float], layer thickness (cm)
    - soilprop_DS: list[float], cumulative depth to bottom of each layer (cm)
    - soilprop_DUL: list[float], drained upper limit (cm3/cm3)
    - soilprop_LL: list[float], lower limit (cm3/cm3)
    - soilprop_NLAYR: int, number of soil layers
    - SW: list[float], initial soil water by layer (cm3/cm3)
    - TAVG: float, daily average air temperature (deg C)
    - TMAX: float, daily max air temperature (deg C)
    - TMIN: float, daily min air temperature (deg C)
    - TAV: float, long-term average temperature (deg C)
    - ORGC_MULCHMASS: float, mulch mass on surface (kg/ha)
    - WATER_SNOW: float, snow water equivalent on surface (mm)
    - X2_PREV: float, previous smoothed surface temperature driver (deg C)

    Returns:
    - CUMDPT: float, profile depth (mm)
    - DSMID: list[float], depth to midpoint of soil layers (mm)
    - TDL: float, sum of (DUL*DLAYR) across layers (cm)
    - TMA: list[float], 5-day memory of surface driver
    - NDays: int, number of days recorded in wet-day memory (0 at init)
    - WetDay: list[int], 30-day memory of wet days (0/1 flags)
    - X2_PREV: float, updated previous driver after spin-up (deg C)
    - SRFTEMP: float, current surface soil temperature (deg C)
    - ST: list[float], soil temperature by layer (deg C)
    """
    import math

    def nint(x: float) -> int:
        return math.floor(x + 0.5) if x >= 0.0 else -math.floor(-x + 0.5)

    NLAYR = soilprop_NLAYR
    BD = soilprop_BD
    DLAYR = soilprop_DLAYR
    DS = soilprop_DS
    DUL = soilprop_DUL
    LL = soilprop_LL

    # Initialize accumulators and depths
    SWI = SW[:]  # initial water per layer
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to midpoint
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm cumulative depth
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if iswitch_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm (use DUL if water not simulated)

    ABD = TBD / DS[NLAYR - 1]  # g/cm3
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    # Initialize moving average array TMA
    TMA_val = nint(TAVG * 10000.0) / 10000.0
    TMA = [TMA_val for _ in range(5)]
    X2_AVG = TMA[0] * 5.0

    # Initialize soil temperature profile
    ST = [TAVG for _ in range(NLAYR)]

    # Wet-day memory
    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    # Soil cover function (mulch / snow)
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG  # will be updated by SOILT_EPIC spin-up
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
            B=B,
            BCV=BCV,
            CUMDPT=CUMDPT,
            DP=DP,
            DSMID=DSMID,
            NLAYR=NLAYR,
            PESW=PESW,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            TMIN=TMIN,
            WetDay=0,
            WFT=WFT,
            WW=WW,
            TMA=TMA,
            ST=ST,
            X2_PREV=X2_PREV,
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    iswitch_ISWWAT: str,
    soilprop_BD: list[float],
    soilprop_DLAYR: list[float],
    soilprop_DS: list[float],
    soilprop_DUL: list[float],
    soilprop_LL: list[float],
    soilprop_NLAYR: int,
    SW: list[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    weather_RAIN: float,
    mgmt_DEPIR: float,
    plant_BIOMAS: float,
    orgc_MULCHMASS: float,
    water_SNOW: float,
    CUMDPT: float,
    DSMID: list[float],
    TDL: float,
    TMA: list[float],
    NDays: int,
    WetDay: list[int],
    X2_PREV: float,
    SRFTEMP: float,
    ST: list[float],
) -> tuple[
    float,           # CUMDPT (unchanged)
    list[float],     # DSMID (unchanged)
    float,           # TDL (updated per original code behavior)
    list[float],     # TMA
    int,             # NDays
    list[int],       # WetDay
    float,           # X2_PREV
    float,           # SRFTEMP
    list[float],     # ST
]:
    """
    Daily EPIC soil temperature process (RATE step).

    Inputs:
    - iswitch_ISWWAT: str ('Y' if water balance simulated, else other)
    - soilprop_BD: list[float], soil bulk density by layer (g/cm3)
    - soilprop_DLAYR: list[float], layer thickness (cm)
    - soilprop_DS: list[float], cumulative depth to bottom of each layer (cm)
    - soilprop_DUL: list[float], drained upper limit (cm3/cm3)
    - soilprop_LL: list[float], lower limit (cm3/cm3)
    - soilprop_NLAYR: int, number of soil layers
    - SW: list[float], soil water by layer for current day (cm3/cm3)
    - TAVG: float, daily average air temperature (deg C)
    - TMAX: float, daily max air temperature (deg C)
    - TMIN: float, daily min air temperature (deg C)
    - TAV: float, long-term average temperature (deg C)
    - weather_RAIN: float, rainfall (mm/d)
    - mgmt_DEPIR: float, irrigation depth (mm/d)
    - plant_BIOMAS: float, plant biomass (kg/ha)
    - orgc_MULCHMASS: float, mulch mass (kg/ha)
    - water_SNOW: float, snow water equivalent (mm)

    InOut States:
    - CUMDPT: float, profile depth (mm) [unchanged here]
    - DSMID: list[float], depth to midpoint of layers (mm) [unchanged]
    - TDL: float, accumulator per original code (cm) [updated without reset]
    - TMA: list[float], 5-day memory of surface driver
    - NDays: int, number of days stored in wet-day memory (<=30)
    - WetDay: list[int], length 30, wet-day flags (0/1)
    - X2_PREV: float, previous smoothed surface temperature driver (deg C)
    - SRFTEMP: float, surface soil temperature (deg C)
    - ST: list[float], soil temperature by layer (deg C)

    Returns updated states in the same order as described above.
    """
    import math

    NLAYR = soilprop_NLAYR
    BD = soilprop_BD
    DLAYR = soilprop_DLAYR
    DS = soilprop_DS
    DUL = soilprop_DUL
    LL = soilprop_LL

    # Accumulate per original logic (note: TDL is not reset in the original RATE code)
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    if isswitch_ISWWAT := iswitch_ISWWAT:
        pass  # no-op; preserve variable name usage
    if iswitch_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    # 30-day memory of wet days (rain + irrigation)
    RAIN = weather_RAIN
    DEPIR = mgmt_DEPIR
    if NDays == 30:
        # shift left by one
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
        # NDays remains 30
    else:
        NDays = NDays + 1
    # set today's wet-day flag
    WetDay_index = NDays - 1 if NDays > 0 else 0
    WetDay[WetDay_index] = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    NWetDays = sum(WetDay)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function
    BIOMAS = plant_BIOMAS
    MULCHMASS = orgc_MULCHMASS
    SNOW = water_SNOW
    CV = (BIOMAS + MULCHMASS) / 1000.0
    denom1 = (CV + math.exp(5.3396 - 2.3951 * CV))
    BCV1 = CV / denom1 if denom1 != 0.0 else 0.0
    denom2 = (SNOW + math.exp(2.303 - 0.2197 * SNOW))
    BCV2 = SNOW / denom2 if denom2 != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Call SOILT_EPIC to update temperatures
    WetDayFlag = WetDay[WetDay_index] if NDays > 0 else 0
    TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
        B=B,
        BCV=BCV,
        CUMDPT=CUMDPT,
        DP=DP,
        DSMID=DSMID,
        NLAYR=NLAYR,
        PESW=PESW,
        TAV=TAV,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        WetDay=WetDayFlag,
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def SOILT_EPIC(
    B: float,
    BCV: float,
    CUMDPT: float,
    DP: float,
    DSMID: list[float],
    NLAYR: int,
    PESW: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    TMIN: float,
    WetDay: int,
    WFT: float,
    WW: float,
    TMA: list[float],
    ST: list[float],
    X2_PREV: float,
) -> tuple[
    list[float],  # TMA
    float,        # SRFTEMP
    list[float],  # ST
    float,        # X2_AVG
    float,        # X2_PREV
]:
    """
    EPIC soil temperature by layer.

    Inputs:
    - B: float, damping coefficient (ln(500/DP))
    - BCV: float, soil cover factor (0-1)
    - CUMDPT: float, total profile depth (mm)
    - DP: float, damping depth (mm)
    - DSMID: list[float], depth to midpoint per layer (mm)
    - NLAYR: int, number of layers
    - PESW: float, plant-extractable soil water (cm)
    - TAV: float, long-term average temperature (deg C)
    - TAVG: float, daily average air temperature (deg C)
    - TMAX: float, daily maximum air temperature (deg C)
    - TMIN: float, daily minimum air temperature (deg C)
    - WetDay: int, flag (1 if wet day, else 0)
    - WFT: float, fraction of wet days over last 30 days
    - WW: float, volumetric water content bound (dimensionless)
    - TMA: list[float], 5-element moving average driver
    - ST: list[float], current soil temperature profile (deg C)
    - X2_PREV: float, previous smoothed driver

    Returns:
    - TMA: updated 5-element driver list
    - SRFTEMP: float, surface temperature (deg C)
    - ST: updated soil temperature by layer
    - X2_AVG: float, average driver over last 5 entries
    - X2_PREV: updated previous driver (set to X2_AVG)
    """
    import math

    # Water content effect (ratio)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0  # PESW[cm], CUMDPT[mm]
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Surface temperature driver (Potter & Williams, 1994)
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update TMA (preserve original Fortran order)
    TMA[0] = X2
    for K in range(4, 0, -1):  # indices 4..1 correspond to Fortran 5..2
        TMA[K] = TMA[K - 1]
    X2_AVG = sum(TMA) / 5.0

    # Surface smoothing with cover
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)

    # Layer temperatures
    X1 = TAV - X3
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_AVG, X2_PREV