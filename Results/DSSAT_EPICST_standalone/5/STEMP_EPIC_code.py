def initialize_STEMP_EPIC(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    ORGC_MULCHMASS,
    WATER_SNOW
):
    """
    Initialization for EPIC soil temperature routine (SEASINIT branch).

    Inputs:
    - ISWWAT: 'Y' if water simulation enabled, else 'N'
    - SOILPROP_BD: list/sequence of bulk density by layer (g/cm3)
    - SOILPROP_DLAYR: list/sequence of layer thickness (cm)
    - SOILPROP_DS: list/sequence of cumulative depth to bottom of each layer (cm)
    - SOILPROP_DUL: list/sequence of water content at drained upper limit (cm3/cm3)
    - SOILPROP_LL: list/sequence of lower limit water content (cm3/cm3)
    - NLAYR: integer number of soil layers to consider
    - SW: list/sequence of initial volumetric water content by layer (cm3/cm3)
    - TAVG: average daily air temperature (degC)
    - TMAX: maximum daily air temperature (degC)
    - TMIN: minimum daily air temperature (degC)
    - TAV: average annual soil temperature (degC)
    - ORGC_MULCHMASS: mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)

    Returns tuple:
    (CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST)

    Where:
    - CUMDPT: cumulative soil profile depth (mm)
    - DSMID: list of depth to midpoint of each soil layer (cm)
    - TDL: total water content of soil at drained upper limit (cm)
    - TMA: list of last 5 days' surface temperature memory
    - NDays: number of days stored in WetDay memory (int, starts at 0)
    - WetDay: list of length 30 with indicators of wet days (0/1)
    - X2_PREV: previous day's average surface temperature estimator (degC)
    - SRFTEMP: current surface temperature (degC)
    - ST: list of soil temperatures by layer (degC)
    """
    from math import exp, log

    # Ensure inputs are lists with appropriate lengths
    BD = list(SOILPROP_BD)[:NLAYR]
    DLAYR = list(SOILPROP_DLAYR)[:NLAYR]
    DS = list(SOILPROP_DS)[:NLAYR]
    DUL = list(SOILPROP_DUL)[:NLAYR]
    LL = list(SOILPROP_LL)[:NLAYR]
    SWI = list(SW)[:NLAYR]

    # Initialize cumulative depth and depths to midpoints
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm depth to mid-point
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm profile depth

    # Sums over layers
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)      # cm
    else:
        PESW = max(0.0, TDL - TLL)      # cm

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX           # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / DP)

    # Initialize TMA memory and initial soil temperature profile
    tavg_round = round(TAVG, 4)
    TMA = [tavg_round for _ in range(5)]
    X2_AVG = TMA[0] * 5.0

    ST = [TAVG for _ in range(NLAYR)]

    # Wet day memory
    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    # Soil cover function (BCV)
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if (CV + exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + exp(2.303 - 0.2197 * SNOW)) if (SNOW + exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    # Initialize X2_PREV to a reasonable value (use TAVG)
    X2_PREV = TAVG
    SRFTEMP = TAVG

    # Spin-up calls to stabilize temperature memory
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
            B=B, BCV=BCV, CUMDPT=CUMDPT, DP=DP, DSMID=DSMID, NLAYR=NLAYR, PESW=PESW, TAV=TAV,
            TAVG=TAVG, TMAX=TMAX, TMIN=TMIN, WetDay=0, WFT=WFT, WW=WW,
            TMA=TMA, ST=ST, X2_PREV=X2_PREV
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
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
    ST
):
    """
    Daily EPIC soil temperature update (RATE branch).

    Inputs:
    - ISWWAT: 'Y' if water simulation enabled, else 'N'
    - SOILPROP_BD: list bulk density by layer (g/cm3)
    - SOILPROP_DLAYR: list layer thickness (cm)
    - SOILPROP_DS: list cumulative depth to bottom of each layer (cm)
    - SOILPROP_DUL: list water content at DUL (cm3/cm3)
    - SOILPROP_LL: list lower limit water content (cm3/cm3)
    - NLAYR: integer number of soil layers
    - SW: list volumetric water content by layer (cm3/cm3)
    - TAVG, TMAX, TMIN: daily temperatures (degC)
    - TAV: average annual soil temperature (degC)
    - WEATHER_RAIN: daily rainfall (mm)
    - MGMT_DEPIR: irrigation depth (mm)
    - PLANT_BIOMAS: aboveground biomass (kg/ha)
    - ORGC_MULCHMASS: mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)
    - CUMDPT: cumulative soil profile depth (mm) [from initialize]
    - DSMID: list depth to midpoint of each soil layer (cm) [from initialize]
    - TDL: running total drained upper limit water content (cm) [state, accumulates per original code]
    - TMA: list of last 5 days' surface temperature memory
    - NDays: number of days stored in WetDay memory (int)
    - WetDay: list of length 30 with indicators of wet days (0/1)
    - X2_PREV: previous day's avg surface temperature estimator (degC)
    - ST: list of soil temperatures by layer (degC)

    Returns tuple:
    (TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST)
    """
    from math import exp, log

    BD = list(SOILPROP_BD)[:NLAYR]
    DLAYR = list(SOILPROP_DLAYR)[:NLAYR]
    DS = list(SOILPROP_DS)[:NLAYR]
    DUL = list(SOILPROP_DUL)[:NLAYR]
    LL = list(SOILPROP_LL)[:NLAYR]
    SW = list(SW)[:NLAYR]
    DSMID = list(DSMID)[:NLAYR]
    ST = list(ST)[:NLAYR]
    TMA = list(TMA)[:5]

    # Sums over layers
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        # Note: The original code does not reset TDL to zero daily; it accumulates across days.
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)      # cm
    else:
        PESW = max(0.0, TDL - TLL)      # cm

    # Wet day memory management (30-day rolling window)
    RAIN = WEATHER_RAIN
    DEPIR = MGMT_DEPIR

    if WetDay is None or len(WetDay) != 30:
        WetDay = [0] * 30
    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
        # NDays remains 30
        idx = 29
    else:
        NDays = NDays + 1
        idx = NDays - 1

    WetDay[idx] = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    NWetDays = sum(WetDay)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function (BCV)
    BIOMAS = PLANT_BIOMAS
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (BIOMAS + MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if (CV + exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + exp(2.303 - 0.2197 * SNOW)) if (SNOW + exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    # Update soil temperatures
    TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
        B=B, BCV=BCV, CUMDPT=CUMDPT, DP=DP, DSMID=DSMID, NLAYR=NLAYR, PESW=PESW, TAV=TAV,
        TAVG=TAVG, TMAX=TMAX, TMIN=TMIN, WetDay=WetDay[idx], WFT=WFT, WW=WW,
        TMA=TMA, ST=ST, X2_PREV=X2_PREV
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
    - B, DP, WW: derived soil thermal parameters
    - BCV: soil cover function (0-1)
    - CUMDPT: cumulative soil profile depth (mm)
    - DSMID: list depth to midpoint of each layer (cm)
    - NLAYR: number of layers
    - PESW: potential extractable soil water over profile (cm)
    - TAV: average annual soil temperature (degC)
    - TAVG, TMAX, TMIN: daily temps (degC)
    - WetDay: 1 if wet day, else 0
    - WFT: fraction of wet days in recent window
    - WW: volumetric water fraction function
    - TMA: list memory of previous 5 days surface temps (in/out)
    - ST: list soil temps by layer (in/out)
    - X2_PREV: previous day's average surface temp estimator (in/out)

    Returns tuple:
    (TMA, SRFTEMP, ST, X2_AVG, X2_PREV)
    """
    from math import exp

    # Ensure lists are copies
    DSMID = list(DSMID)[:NLAYR]
    ST = list(ST)[:NLAYR]
    TMA = list(TMA)[:5]

    # Water content ratio (dimensionless)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Determine surface temp estimator X2
    if WetDay > 0:
        # Eqn 2: X2 = WFT*(TAVG - TMIN) + TMIN
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        # Eqn 1: X2 = WFT*(TMAX - TAVG) + TAVG + 2.
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update temperature memory (note original order: set then shift)
    TMA[0] = X2
    for K in range(4, 0, -1):  # 5->2 in 1-based indexing; here 4->1 zero-based
        TMA[K] = TMA[K - 1]

    X2_AVG = sum(TMA) / 5.0

    # Blend with previous day considering cover (Eqn 4)
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    # Surface temperature (min of X2_AVG and X3)
    SRFTEMP = min(X2_AVG, X3)

    # Layer factor and update soil temperature (Eqns 6-8)
    LAG = 0.5
    X1 = TAV - X3
    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    # Update previous day's average for next step
    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_AVG, X2_PREV


def test_STEMP_EPIC_initialize_example():
    """
    Test derived from ASKEE_DSSAT_EPIC initialization scenario.
    Does not assert numerical values; returns outputs to allow inspection.
    """
    ISWWAT = 'Y'
    NLAYR = 4
    SOILPROP_BD = [1.6] * NLAYR
    SOILPROP_DLAYR = [10.0] * NLAYR
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * NLAYR
    SOILPROP_LL = [0.2] * NLAYR
    SW = [0.2] * NLAYR
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0
    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0

    return initialize_STEMP_EPIC(
        ISWWAT=ISWWAT,
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        NLAYR=NLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW
    )


def test_STEMP_EPIC_daily_example():
    """
    Test derived from ASKEE_DSSAT_EPIC daily scenario for one day.
    Returns outputs to allow inspection.
    """
    # Initialize
    init = test_STEMP_EPIC_initialize_example()
    (CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP_init, ST) = init

    # Daily inputs
    ISWWAT = 'Y'
    NLAYR = 4
    SOILPROP_BD = [1.6] * NLAYR
    SOILPROP_DLAYR = [10.0] * NLAYR
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * NLAYR
    SOILPROP_LL = [0.2] * NLAYR
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

    return STEMP_EPIC(
        ISWWAT=ISWWAT,
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
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
        ST=ST
    )