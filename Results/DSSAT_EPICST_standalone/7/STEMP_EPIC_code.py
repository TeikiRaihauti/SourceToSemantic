from math import exp, log

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
    SRFTEMP,
    ST,
    X2_AVG,
    X2_PREV,
):
    """
    Supporting function implementing EPIC soil temperature by layer.
    Mirrors the Fortran SOILT_EPIC subroutine logic.

    Parameters:
    - B, BCV, CUMDPT, DP, NLAYR, PESW, TAV, TAVG, TMAX, TMIN, WetDay, WFT, WW: scalars (floats or ints)
    - DSMID: list of floats (length >= NLAYR)
    - TMA: list of 5 floats (in/out)
    - SRFTEMP: float (ignored input, computed as output)
    - ST: list of floats (length >= NLAYR) (in/out)
    - X2_AVG: float (ignored input, computed as output)
    - X2_PREV: float (in/out)

    Returns:
    - TMA (updated list of 5 floats)
    - SRFTEMP (float)
    - ST (updated list of floats length NLAYR)
    - X2_AVG (float)
    - X2_PREV (updated float)
    """
    # Ensure copies so inputs are not mutated unexpectedly
    TMA_loc = list(TMA)
    ST_loc = list(ST)

    # Water content ratio (WC)
    denom = WW * CUMDPT
    if denom <= 0.0:
        # Guard against division by zero; Fortran code assumes positive denom given valid inputs
        WC = 0.01
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Daily driving temperature, Potter & Williams (1994)
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day moving array (note: follows Fortran order and side-effects)
    TMA_loc[0] = X2
    for K in range(4, 0, -1):  # 5..2 in Fortran 1-based
        TMA_loc[K] = TMA_loc[K - 1]
    X2_AVG_loc = sum(TMA_loc) / 5.0

    # Surface litter/cover temperature mix
    X3 = (1.0 - BCV) * X2_AVG_loc + BCV * X2_PREV

    SRFTEMP_loc = min(X2_AVG_loc, X3)

    # Offset from annual average soil temperature
    X1 = TAV - X3

    # Layer temperatures
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD > 0.0 else 0.0
        F = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        ST_loc[L] = LAG * ST_loc[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_loc = X2_AVG_loc

    return TMA_loc, SRFTEMP_loc, ST_loc, X2_AVG_loc, X2_PREV_loc


def STEMP_EPIC_init(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_TAMP,
    ORGC_MULCHMASS,
    WATER_SNOW,
):
    """
    Initialization function corresponding to the SEASINIT branch of STEMP_EPIC.
    Computes initial soil thermal state and supporting variables.

    Inputs:
    - ISWWAT: 'Y' if water simulation is on, else other
    - SOILPROP_BD: list of bulk density by layer (g/cm3)
    - SOILPROP_DLAYR: list of layer thickness (cm)
    - SOILPROP_DS: list of cumulative depths (cm); SOILPROP_DS[NLAYR-1] used
    - SOILPROP_DUL: list of drained upper limit by layer (cm3/cm3)
    - SOILPROP_LL: list of lower limit by layer (cm3/cm3)
    - SOILPROP_NLAYR: number of layers (int)
    - SW: list of initial water content by layer (cm3/cm3)
    - TAVG: average daily air temperature (C)
    - TMAX, TMIN: daily max/min air temperature (C)
    - TAV: average annual air temperature (C)
    - WEATHER_TAMP: annual amplitude (degC) [not used by EPIC method]
    - ORGC_MULCHMASS: surface mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)

    Returns:
    - CUMDPT: cumulative soil profile depth (mm)
    - DSMID: list of depths to midpoints of layers (mm)
    - TDL: total water at DUL over profile (cm)
    - TMA: list of last 5 days of X2 (degC proxy)
    - NDays: count of days in wet-day memory (int)
    - WetDay: list of length 30, 0/1 flags
    - X2_PREV: previous surface temperature driver (degC)
    - SRFTEMP: surface temperature (degC)
    - ST: list of soil temperatures by layer (degC)
    """
    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    SWI = list(SW)

    # Initialize cumulative depth and layer midpoints (mm)
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to layer midpoint
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm cumulative depth
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Average bulk density over the profile
    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / DP)

    # Initialize TMA with TAVG to 4 decimals, and ST to TAVG
    TMA = [round(TAVG, 4) for _ in range(5)]
    X2_AVG = TMA[0] * 5.0
    ST = [TAVG for _ in range(NLAYR)]

    # 30-day wet-day memory
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Surface cover factor (mulch or snow)
    CV = (ORGC_MULCHMASS) / 1000.0
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if CV + exp(5.3396 - 2.3951 * CV) != 0.0 else 0.0
    BCV2 = WATER_SNOW / (WATER_SNOW + exp(2.303 - 0.2197 * WATER_SNOW)) if WATER_SNOW + exp(2.303 - 0.2197 * WATER_SNOW) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG  # placeholder initial value
    X2_PREV = 0.0   # will be set within SOILT_EPIC; initial value not used if BCV=0

    # Spin-up calls to stabilize soil temperatures and TMA
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
            SRFTEMP=SRFTEMP,
            ST=ST,
            X2_AVG=X2_AVG,
            X2_PREV=X2_PREV,
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
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
    ST,
):
    """
    Main daily biophysical process function corresponding to RATE branch of STEMP_EPIC.

    Inputs:
    - ISWWAT: 'Y' if water is simulated
    - SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL: lists per layer
    - SOILPROP_NLAYR: int
    - SW: list of current volumetric water content per layer (cm3/cm3)
    - TAVG, TMAX, TMIN, TAV: daily/annual temperatures (C)
    - WEATHER_RAIN: daily rainfall (mm)
    - MGMT_DEPIR: irrigation depth (mm)
    - PLANT_BIOMAS: biomass (kg/ha)
    - ORGC_MULCHMASS: surface mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)
    - CUMDPT: cumulative profile depth (mm)
    - DSMID: list of depths to layer midpoints (mm)
    - TDL: running total DUL water over profile (cm); note: replicates Fortran behavior (accumulates daily)
    - TMA: list of 5 floats (state memory)
    - NDays: int
    - WetDay: list of length 30 with 0/1 history
    - X2_PREV: float state
    - ST: list of soil temps per layer (degC)

    Returns:
    - TDL: updated
    - TMA: updated list length 5
    - NDays: updated int
    - WetDay: updated list length 30
    - X2_PREV: updated float
    - SRFTEMP: surface temperature (degC)
    - ST: updated soil temperature profile (list)
    """
    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    SW_now = list(SW)
    ST_loc = list(ST)
    TMA_loc = list(TMA)
    WetDay_loc = list(WetDay)

    # Accumulate bulk density, water contents across profile
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # Note: cumulative across days, preserving Fortran behavior
        TLL += LL[L] * DLAYR[L]
        TSW += SW_now[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Update wet-day memory (rainfall + irrigation)
    if NDays == 30:
        # Shift entries left; keep length 30
        for i in range(29):
            WetDay_loc[i] = WetDay_loc[i + 1]
        WetDay_loc[29] = 0
    else:
        NDays += 1

    wet_today = 1 if (WEATHER_RAIN + MGMT_DEPIR) > 1.0e-6 else 0
    WetDay_loc[NDays - 1] = wet_today  # Python 0-based index

    NWetDays = sum(WetDay_loc)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Surface cover function
    CV = (PLANT_BIOMAS + ORGC_MULCHMASS) / 1000.0
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if CV + exp(5.3396 - 2.3951 * CV) != 0.0 else 0.0
    BCV2 = WATER_SNOW / (WATER_SNOW + exp(2.303 - 0.2197 * WATER_SNOW)) if WATER_SNOW + exp(2.303 - 0.2197 * WATER_SNOW) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Call supporting temperature routine
    SRFTEMP = 0.0  # placeholder input; will be computed by SOILT_EPIC
    X2_AVG = 0.0
    TMA_loc, SRFTEMP, ST_loc, X2_AVG, X2_PREV = SOILT_EPIC(
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
        WetDay=wet_today,
        WFT=WFT,
        WW=WW,
        TMA=TMA_loc,
        SRFTEMP=SRFTEMP,
        ST=ST_loc,
        X2_AVG=X2_AVG,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA_loc, NDays, WetDay_loc, X2_PREV, SRFTEMP, ST_loc


# ----------------------------- Tests ----------------------------------

def test_STEMP_EPIC_basic():
    """
    Basic test derived from ASKEE_DSSAT_EPIC.for program.
    Constructs inputs and runs initialization and one daily step.
    """
    # Soil properties (4 layers)
    SOILPROP_NLAYR = 4
    SOILPROP_BD = [1.6] * SOILPROP_NLAYR
    SOILPROP_DLAYR = [10.0] * SOILPROP_NLAYR  # cm
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]    # cm cumulative
    SOILPROP_DUL = [0.3] * SOILPROP_NLAYR
    SOILPROP_LL = [0.2] * SOILPROP_NLAYR

    # Initial conditions and weather
    ISWWAT = 'Y'
    SW = [0.2] * SOILPROP_NLAYR
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0
    WEATHER_TAMP = 10.0
    WEATHER_RAIN = 0.0

    # Management and cover
    MGMT_DEPIR = 0.0
    PLANT_BIOMAS = 0.0
    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0

    # Initialize
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC_init(
        ISWWAT=ISWWAT,
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        SOILPROP_NLAYR=SOILPROP_NLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        WEATHER_TAMP=WEATHER_TAMP,
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW,
    )

    # Sanity checks after init
    assert len(DSMID) == SOILPROP_NLAYR
    assert len(TMA) == 5
    assert len(WetDay) == 30
    assert len(ST) == SOILPROP_NLAYR

    # Daily rate step
    TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC(
        ISWWAT=ISWWAT,
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        SOILPROP_NLAYR=SOILPROP_NLAYR,
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
        ST=ST,
    )

    # Sanity checks after one day
    assert NDays == 1
    assert WetDay[0] in (0, 1)
    assert isinstance(SRFTEMP, float)
    assert all(isinstance(x, float) for x in ST)