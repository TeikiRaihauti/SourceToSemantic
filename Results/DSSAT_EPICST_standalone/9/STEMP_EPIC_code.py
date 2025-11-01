def STEMP_EPIC_initialize(
    ISWWAT,
    BD,
    DLAYR,
    DS,
    DUL,
    LL,
    NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_TAMP,
    ORGC_MULCHMASS,
    WATER_SNOW,
    X2_PREV
):
    """
    Initialization for EPIC soil temperature routine (SEASINIT branch).

    Inputs:
    - ISWWAT: 'Y' if water simulation is on, otherwise 'N'
    - BD: list of bulk densities per layer (g/cm3), length >= NLAYR
    - DLAYR: list of layer thicknesses (cm), length >= NLAYR
    - DS: list of cumulative depths (cm) per layer, length >= NLAYR
    - DUL: list of drained upper limit (vol/vol), length >= NLAYR
    - LL: list of lower limit (vol/vol), length >= NLAYR
    - NLAYR: number of soil layers
    - SW: list of initial soil water (vol/vol), length >= NLAYR
    - TAVG: average daily temperature (degC)
    - TMAX: maximum daily temperature (degC)
    - TMIN: minimum daily temperature (degC)
    - TAV: average annual soil temperature (degC)
    - WEATHER_TAMP: annual amplitude of air temperature (degC) [not used in core calc]
    - ORGC_MULCHMASS: mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)
    - X2_PREV: previous day's surface temperature memory (degC)

    Returns:
    - CUMDPT: cumulative profile depth (mm)
    - DSMID: list of depths to midpoint of each layer (mm), length NLAYR
    - TDL: total water content at DUL over profile (cm)
    - TMA: 5-day memory array for surface temp driver, length 5
    - NDays: number of wet-day memory days initialized (int)
    - WetDay: 30-length list of wet-day flags (0/1)
    - X2_PREV: updated previous surface temperature (degC)
    - SRFTEMP: surface temperature (degC)
    - ST: list of soil temperatures by layer (degC), length NLAYR
    """
    import math

    def nint(x):
        return int(math.floor(x + 0.5)) if x >= 0.0 else int(math.ceil(x - 0.5))

    # Initialize depth-related variables
    CUMDPT = 0.0
    DSMID = [0.0 for _ in range(NLAYR)]
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to midpoint of layer
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm cumulative depth
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    # Potential extractable soil water (cm)
    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Soil thermal parameters
    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD   # volumetric fraction
    B = math.log(500.0 / DP)

    # Initialize TMA memory
    TMA_val = nint(TAVG * 10000.0) / 10000.0
    TMA = [TMA_val for _ in range(5)]
    # X2_AVG initial placeholder (not used before updated in SOILT_EPIC)
    # Initialize soil layer temperatures
    ST = [TAVG for _ in range(NLAYR)]

    # Wet-day memory and fraction
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Soil cover function (mulch + snow)
    CV = (ORGC_MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    SNOW = WATER_SNOW
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG  # initialize

    # Spin-up calls
    for _ in range(8):
        SRFTEMP, ST, X2_AVG, X2_PREV, TMA = SOILT_EPIC(
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
            X2_PREV=X2_PREV
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    ISWWAT,
    BD,
    DLAYR,
    DS,
    DUL,
    LL,
    NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_RAIN,
    MGMT_DEPIR,
    PLANT_BIOMAS,
    ORGC_MULCHMASS,
    WATER_SNOW,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    NDays,
    WetDay,
    X2_PREV,
    SRFTEMP,
    ST
):
    """
    Daily EPIC soil temperature routine (RATE branch).

    Inputs:
    - ISWWAT: 'Y' if water simulation is on, otherwise 'N'
    - BD: list of bulk densities per layer (g/cm3), length >= NLAYR
    - DLAYR: list of layer thicknesses (cm), length >= NLAYR
    - DS: list of cumulative depths (cm) per layer, length >= NLAYR
    - DUL: list of drained upper limit (vol/vol), length >= NLAYR
    - LL: list of lower limit (vol/vol), length >= NLAYR
    - NLAYR: number of soil layers
    - SW: list of current soil water (vol/vol), length >= NLAYR
    - TAVG: average daily temperature (degC)
    - TMAX: maximum daily temperature (degC)
    - TMIN: minimum daily temperature (degC)
    - TAV: average annual soil temperature (degC)
    - WEATHER_RAIN: rainfall (mm)
    - MGMT_DEPIR: irrigation (mm)
    - PLANT_BIOMAS: plant biomass (kg/ha)
    - ORGC_MULCHMASS: mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)
    - CUMDPT: cumulative profile depth (mm)
    - DSMID: list of depths to midpoint of each layer (mm), length NLAYR
    - TDL: total water content at DUL over profile (cm) [state]
    - TMA: 5-day memory array for surface temp driver, length 5 [state]
    - NDays: number of days recorded in wet-day memory [state]
    - WetDay: length-30 list of wet-day flags [state]
    - X2_PREV: previous day's smoothed surface temp (degC) [state]
    - SRFTEMP: previous surface temperature (degC) [state]
    - ST: previous soil temperatures by layer (degC) [state]

    Returns:
    - TDL: updated total water content at DUL over profile (cm)
    - TMA: updated 5-day memory array
    - NDays: updated number of days in wet-day memory (<= 30)
    - WetDay: updated wet-day memory list
    - X2_PREV: updated previous smoothed surface temp (degC)
    - SRFTEMP: updated surface temperature (degC)
    - ST: updated soil temperatures by layer (degC)
    """
    import math

    # Sums over profile
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # Note: matches Fortran RATE block behavior
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    # Soil thermal parameters
    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    # Potential extractable soil water (cm)
    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Wet-day memory update (30-day window)
    RAIN = WEATHER_RAIN
    DEPIR = MGMT_DEPIR
    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
    else:
        NDays = NDays + 1
    WetDay_val = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    WetDay[NDays - 1] = WetDay_val
    NWetDays = sum(WetDay[:NDays])
    WFT = float(NWetDays) / float(NDays)

    # Soil cover function (biomass + mulch + snow)
    BIOMAS = PLANT_BIOMAS
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (BIOMAS + MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP, ST, X2_AVG, X2_PREV, TMA = SOILT_EPIC(
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
        WetDay=WetDay[NDays - 1],
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def SOILT_EPIC(
    B,
    BCV,
    CUMDPT,
    DP,
    DSMID,
    NLAYR,
    PESW,
    TAV,
    TAVG,
    TMAX,
    TMIN,
    WetDay,
    WFT,
    WW,
    TMA,
    ST,
    X2_PREV
):
    """
    EPIC soil temperature by layer.

    Inputs:
    - B: exponential decay factor
    - BCV: soil cover fraction (0-1)
    - CUMDPT: cumulative soil depth (mm)
    - DP: damping depth (mm)
    - DSMID: list of mid-layer depths (mm), length NLAYR
    - NLAYR: number of soil layers
    - PESW: potential extractable soil water over profile (cm)
    - TAV: average annual soil temperature (degC)
    - TAVG: average daily temperature (degC)
    - TMAX: daily maximum temperature (degC)
    - TMIN: daily minimum temperature (degC)
    - WetDay: 1 if wet day, else 0
    - WFT: fraction of wet days (0-1)
    - WW: volumetric water content parameter
    - TMA: memory array (5 elements)
    - ST: list of soil temperatures by layer
    - X2_PREV: previous smoothed surface temperature (degC)

    Returns:
    - SRFTEMP: updated surface temperature (degC)
    - ST: updated soil temperatures by layer (degC)
    - X2_AVG: average driver temperature (degC)
    - X2_PREV: updated smoothed surface temperature memory (degC)
    - TMA: updated memory array
    """
    import math

    # Water content ratio
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Compute X2 based on wet/dry day classification
    if WetDay > 0:
        # Eqn. 2: X2 = WFT*(TAVG - TMIN) + TMIN
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        # Eqn. 1 (modified): X2 = WFT*(TMAX - TAVG) + TAVG + 2.
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update TMA memory (replicating Fortran order)
    TMA_new = list(TMA)
    TMA_new[0] = X2
    for K in range(4, 0, -1):  # indices 4->1 correspond to Fortran 5->2
        TMA_new[K] = TMA_new[K - 1]
    X2_AVG = sum(TMA_new) / 5.0

    # Surface cover mixing
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)

    # Temperature offset from annual mean
    X1 = TAV - X3

    # Layer temperature updates with lag
    LAG = 0.5
    ST_new = list(ST)
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0  # avoid division by zero
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        ST_new[L] = LAG * ST_new[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_new = X2_AVG

    return SRFTEMP, ST_new, X2_AVG, X2_PREV_new, TMA_new


def test_STEMP_EPIC_example():
    """
    Derived test from ASKEE_DSSAT_EPIC.for program.
    Runs initialization and a few daily steps to ensure outputs are reasonable.
    """
    # Inputs analogous to ASKEE
    ISWWAT = 'Y'
    NLAYR = 4
    BD = [1.6 for _ in range(NLAYR)]
    DLAYR = [10.0 for _ in range(NLAYR)]
    DS = [10.0, 20.0, 30.0, 40.0]
    DUL = [0.3 for _ in range(NLAYR)]
    LL = [0.2 for _ in range(NLAYR)]
    SW = [0.2 for _ in range(NLAYR)]

    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0
    WEATHER_TAMP = 10.0
    WEATHER_RAIN = 0.0

    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0
    MGMT_DEPIR = 0.0
    PLANT_BIOMAS = 0.0
    X2_PREV = 0.0

    # Initialize
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC_initialize(
        ISWWAT=ISWWAT,
        BD=BD,
        DLAYR=DLAYR,
        DS=DS,
        DUL=DUL,
        LL=LL,
        NLAYR=NLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        WEATHER_TAMP=WEATHER_TAMP,
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW,
        X2_PREV=X2_PREV
    )

    # Simple assertions on initialization
    assert len(DSMID) == NLAYR
    assert len(ST) == NLAYR
    assert len(TMA) == 5
    assert len(WetDay) == 30

    # Run a few daily steps
    for _ in range(5):
        TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC(
            ISWWAT=ISWWAT,
            BD=BD,
            DLAYR=DLAYR,
            DS=DS,
            DUL=DUL,
            LL=LL,
            NLAYR=NLAYR,
            SW=SW,
            TAVG=TAVG,
            TMAX=TMAX,
            TMIN=TMIN,
            TAV=TAV,
            WEATHER_RAIN=WEATHER_RAIN,
            MGMT_DEPIR=MGMT_DEPIR,
            PLANT_BIOMAS=PLANT_BIOMAS,
            ORGC_MULCHMASS=ORGC_MULCHMASS,
            WATER_SNOW=WATER_SNOW,
            CUMDPT=CUMDPT,
            DSMID=DSMID,
            TDL=TDL,
            TMA=TMA,
            NDays=NDays,
            WetDay=WetDay,
            X2_PREV=X2_PREV,
            SRFTEMP=SRFTEMP,
            ST=ST
        )

    # Basic sanity checks after updates
    assert len(ST) == NLAYR
    assert isinstance(SRFTEMP, float)
    assert len(TMA) == 5
    assert 1 <= NDays <= 30
    assert all(isinstance(w, int) for w in WetDay)