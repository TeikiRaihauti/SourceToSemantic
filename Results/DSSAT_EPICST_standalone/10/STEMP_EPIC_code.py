import math

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
    X2_PREV,
):
    # Inputs:
    #   B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV, TAVG, TMAX, TMIN, WetDay, WFT, WW
    # InOut:
    #   TMA (list length 5), ST (list length NLAYR), X2_PREV (float)
    # Output:
    #   SRFTEMP (float), X2_PREV (float), TMA (updated), ST (updated)

    # Convert to local mutable copies where needed
    # WC calculation (dimensionless)
    # PESW in cm, CUMDPT in mm, WW dimensionless
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Daily surface temperature forcing according to EPIC method
    if WetDay > 0:
        # Potter & Williams (1994) Eqn. 2
        # X2 = WFT*(TX - TMN) + TMN
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        # Potter & Williams (1994) Eqn. 1 (without ST0 factor, +2 offset used)
        # X2 = WFT*(TMX - TX) + TX + 2
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory TMA (most recent in position 0)
    # Shift previous values towards end
    for k in range(4, 0, -1):
        TMA[k] = TMA[k - 1]
    TMA[0] = X2

    X2_AVG = sum(TMA) / 5.0

    # Soil cover effect
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    # Surface litter temperature (degC)
    SRFTEMP = min(X2_AVG, X3)

    # Thermal wave lag factor
    LAG = 0.5

    # Temperature offset
    X1 = TAV - X3

    # Update soil layer temperatures
    for L in range(NLAYR):
        ZD = DSMID[L] / DD  # depth index
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    # Store average forcing as "previous" for next timestep
    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_PREV


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
    MULCHMASS,
    SNOW,
    X2_PREV,
):
    # Initialize soil thermal state for EPIC soil temperature computation.
    # Parameters:
    #   ISWWAT: 'Y' or 'N'
    #   BD, DLAYR, DS, DUL, LL, SW: lists of length >= NLAYR
    #   NLAYR: int
    #   TAVG, TMAX, TMIN, TAV: floats
    #   MULCHMASS, SNOW: floats
    #   X2_PREV: float (initial previous-day surface temperature forcing)
    #
    # Returns:
    #   CUMDPT (float, mm),
    #   DSMID (list, length NLAYR, mm),
    #   TDL (float, cm),
    #   TMA (list, length 5),
    #   NDays (int),
    #   WetDay (list, length 30 ints 0/1),
    #   X2_PREV (float, updated),
    #   SRFTEMP (float),
    #   ST (list, length NLAYR, degC)

    # Initialize cumulative profile depth and DSMID (mm)
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0

    # Initial soil water profile to use for PESW at initialization
    SWI = list(SW)

    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to midpoint of layer
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm profile depth increment
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm (use DUL when water not simulated)

    ABD = TBD / DS[NLAYR - 1]  # g/cm3
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD   # volumetric fraction
    B = math.log(500.0 / DP)

    # Initialize 5-day memory with air TAVG
    TMA = [round(TAVG, 4) for _ in range(5)]
    # ST layers initial at TAVG
    ST = [TAVG for _ in range(NLAYR)]

    # Wet-day memory (30-day), fraction WFT
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Soil cover function (t/ha from kg/ha)
    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV > 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if SNOW > 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG  # will be updated by spin-up calls

    # Spin-up 8 steps to initialize temperature memory and layers
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_PREV = SOILT_EPIC(
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
            X2_PREV=X2_PREV,
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC_rate(
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
    RAIN,
    DEPIR,
    BIOMAS,
    MULCHMASS,
    SNOW,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    NDays,
    WetDay,
    X2_PREV,
    SRFTEMP,
    ST,
):
    # Perform daily EPIC soil temperature update (RATE phase).
    # Parameters:
    #   ISWWAT: 'Y' or 'N'
    #   BD, DLAYR, DS, DUL, LL, SW: lists length >= NLAYR
    #   NLAYR: int
    #   TAVG, TMAX, TMIN, TAV, RAIN, DEPIR, BIOMAS, MULCHMASS, SNOW: floats
    #   State In:
    #     CUMDPT (mm), DSMID (mm list), TDL (cm), TMA (5 list),
    #     NDays (int), WetDay (30 list), X2_PREV (float), SRFTEMP (float),
    #     ST (NLAYR list)
    #
    # Returns updated state:
    #   TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: The original code accumulates TDL without resetting in RATE.
    # This behavior is preserved here for fidelity.
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX   # mm
    WW = 0.356 - 0.144 * ABD    # volumetric
    B = math.log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    # Update 30-day wet day memory
    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
    else:
        NDays += 1

    if (RAIN + DEPIR) > 1.0e-6:
        WetDay[NDays - 1] = 1
    else:
        WetDay[NDays - 1] = 0

    NWetDays = sum(WetDay[:NDays])
    WFT = float(NWetDays) / float(NDays)

    # Soil cover function (t/ha from kg/ha)
    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV > 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if SNOW > 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Call soil temperature profile update
    TMA, SRFTEMP, ST, X2_PREV = SOILT_EPIC(
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
        SRFTEMP=SRFTEMP,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def test_STEMP_EPIC_runs():
    # This test reproduces the ASKEE demonstration program scenario.
    # It initializes a simple 4-layer soil and runs initialization
    # followed by a single RATE update.

    # Switches
    ISWWAT = 'Y'

    # Soil properties (4 layers)
    NLAYR = 4
    BD = [1.6, 1.6, 1.6, 1.6]
    DLAYR = [10.0, 10.0, 10.0, 10.0]  # cm
    DS = [10.0, 20.0, 30.0, 40.0]     # cm (cumulative)
    DUL = [0.3, 0.3, 0.3, 0.3]
    LL = [0.2, 0.2, 0.2, 0.2]

    # Initial water content
    SW = [0.2, 0.2, 0.2, 0.2]

    # Weather and temps
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0

    # Cover and water inputs
    MULCHMASS = 0.0
    SNOW = 0.0
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
        MULCHMASS=MULCHMASS,
        SNOW=SNOW,
        X2_PREV=X2_PREV,
    )

    # Daily rate inputs
    RAIN = 0.0
    DEPIR = 0.0
    BIOMAS = 0.0

    # One daily update
    TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC_rate(
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
        RAIN=RAIN,
        DEPIR=DEPIR,
        BIOMAS=BIOMAS,
        MULCHMASS=MULCHMASS,
        SNOW=SNOW,
        CUMDPT=CUMDPT,
        DSMID=DSMID,
        TDL=TDL,
        TMA=TMA,
        NDays=NDays,
        WetDay=WetDay,
        X2_PREV=X2_PREV,
        SRFTEMP=SRFTEMP,
        ST=ST,
    )

    # Return a subset of results for potential assertions by a test harness
    return {
        "CUMDPT": CUMDPT,
        "DSMID": DSMID,
        "TDL": TDL,
        "TMA": TMA,
        "NDays": NDays,
        "WetDay": WetDay,
        "X2_PREV": X2_PREV,
        "SRFTEMP": SRFTEMP,
        "ST": ST,
    }