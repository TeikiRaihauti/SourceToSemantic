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
) -> tuple[list[float], float, list[float], float, float]:
    """
    SOILT_EPIC
    Determines soil temperature by layer (EPIC method).

    Inputs:
    - B: float, dimensionless shape factor (ln(500/DP)).
    - BCV: float, soil cover factor from mulch/canopy/snow (0..1).
    - CUMDPT: float, profile depth (mm).
    - DP: float, damping depth (mm).
    - DSMID: list[float], depth to midpoints of layers (mm), length >= NLAYR.
    - NLAYR: int, number of soil layers.
    - PESW: float, plant extractable soil water in profile (cm).
    - TAV: float, annual average air temperature (C).
    - TAVG: float, daily average air temperature (C).
    - TMAX: float, daily maximum air temperature (C).
    - TMIN: float, daily minimum air temperature (C).
    - WetDay: int, 1 if wet day (rain or irrigation), else 0.
    - WFT: float, fraction of wet days over last 30 days (0..1).
    - WW: float, water content parameter = 0.356 - 0.144*ABD (dimensionless).
    - TMA: list[float], moving 5-day memory variable (length 5), in/out.
    - ST: list[float], soil temperature by layer (C), length >= NLAYR, in/out.
    - X2_PREV: float, previous day's X2_AVG (C), in/out.

    Returns:
    - TMA: updated list[float], length 5.
    - SRFTEMP: float, surface temperature (C).
    - ST: updated list[float], soil temperature by layer (C), length NLAYR.
    - X2_AVG: float, updated 5-day average of X2 (C).
    - X2_PREV: float, updated to current X2_AVG (C).
    """
    import math

    # Water content ratio, PESW in cm; WW*CUMDPT in mm -> convert with *10
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Surface drivers (Potter & Williams, 1994)
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory (ordering matches original Fortran sequence)
    # Note: This sequence duplicates X2 in first two positions, preserving original behavior.
    TMA = list(TMA)
    TMA[0] = X2
    for K in range(4, 0, -1):
        TMA[K] = TMA[K - 1]
    X2_AVG = sum(TMA) / 5.0

    # Surface temperature with cover effect
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)

    # Temperature lag and damping with depth
    LAG = 0.5
    X1 = TAV - X3
    ST = list(ST)
    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_AVG, X2_PREV


def stemp_epic_initialize(
    ISWWAT: str,
    BD: list[float],
    DLAYR: list[float],
    DS: list[float],
    DUL: list[float],
    LL: list[float],
    NLAYR: int,
    SW: list[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    MULCHMASS: float,
    SNOW: float,
    X2_PREV: float = 0.0,
) -> tuple[
    float,  # CUMDPT
    list[float],  # DSMID
    float,  # TDL
    list[float],  # TMA
    int,  # NDays
    list[int],  # WetDay
    float,  # X2_PREV
    float,  # SRFTEMP
    list[float],  # ST
    float,  # X2_AVG
]:
    """
    Seasonal initialization for EPIC soil temperature.

    Inputs:
    - ISWWAT: str, 'Y' if water balance simulated, else 'N'.
    - BD: list[float], bulk density per layer (g/cm3), length >= NLAYR.
    - DLAYR: list[float], layer thickness (cm), length >= NLAYR.
    - DS: list[float], cumulative depth to bottom of each layer (cm), length >= NLAYR.
    - DUL: list[float], drained upper limit per layer (cm3/cm3), length >= NLAYR.
    - LL: list[float], lower limit per layer (cm3/cm3), length >= NLAYR.
    - NLAYR: int, number of soil layers.
    - SW: list[float], initial volumetric water content per layer (cm3/cm3), length >= NLAYR.
    - TAVG: float, today's average air temperature (C).
    - TMAX: float, today's maximum air temperature (C).
    - TMIN: float, today's minimum air temperature (C).
    - TAV: float, annual average air temperature (C).
    - MULCHMASS: float, surface mulch mass (kg/ha).
    - SNOW: float, snow water equivalent (mm).
    - X2_PREV: float, previous X2_AVG (C), default 0.0.

    Returns:
    - CUMDPT: float, profile depth to bottom (mm).
    - DSMID: list[float], depth to midpoint of each layer (mm), length NLAYR.
    - TDL: float, sum(DUL*DLAYR) (cm) initialized.
    - TMA: list[float], 5-day memory vector initialized.
    - NDays: int, initialized to 0.
    - WetDay: list[int], length 30, initialized to zeros.
    - X2_PREV: float, updated after spin-up.
    - SRFTEMP: float, surface temperature after spin-up (C).
    - ST: list[float], soil temperature by layer after spin-up (C), length NLAYR.
    - X2_AVG: float, last computed 5-day average of X2 (C).
    """
    import math

    # Initializations
    SWI = list(SW)
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0

    # Depths and water sums
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to layer midpoint
        CUMDPT += DLAYR[L] * 10.0  # mm to bottom of profile
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    # Plant extractable soil water (cm)
    if ISWWAT.upper() == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Soil parameters
    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    # Temperature memory (rounded to 4 decimals as in Fortran NINT(*1e4)/1e4)
    ta_round = round(TAVG, 4)
    TMA = [ta_round] * 5
    X2_AVG = TMA[0] * 5.0  # Preserves original assignment (not used before recomputation)

    # Initialize soil temperatures
    ST = [TAVG] * NLAYR

    # Wet-day memory
    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    # Soil cover factor from mulch and snow
    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV + math.exp(5.3396 - 2.3951 * CV) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if SNOW + math.exp(2.303 - 0.2197 * SNOW) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Spin-up (8 iterations) with dry-day assumption
    SRFTEMP = TAVG
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

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, X2_AVG


def stemp_epic_process(
    ISWWAT: str,
    BD: list[float],
    DLAYR: list[float],
    DS: list[float],
    DUL: list[float],
    LL: list[float],
    NLAYR: int,
    SW: list[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    RAIN: float,
    DEPIR: float,
    BIOMAS: float,
    MULCHMASS: float,
    SNOW: float,
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
    float,  # TDL (updated)
    list[float],  # TMA
    int,  # NDays
    list[int],  # WetDay
    float,  # X2_PREV
    float,  # SRFTEMP
    list[float],  # ST
    float,  # X2_AVG
]:
    """
    Daily process for EPIC soil temperature.

    Inputs:
    - ISWWAT: str, 'Y' if water balance simulated, else 'N'.
    - BD, DLAYR, DS, DUL, LL: lists[float], soil properties per layer, length >= NLAYR.
    - NLAYR: int, number of soil layers.
    - SW: list[float], daily volumetric water content per layer (cm3/cm3).
    - TAVG, TMAX, TMIN: float, daily temperatures (C).
    - TAV: float, annual average air temperature (C).
    - RAIN: float, daily precipitation (mm).
    - DEPIR: float, daily irrigation (mm).
    - BIOMAS: float, live biomass (kg/ha).
    - MULCHMASS: float, mulch mass (kg/ha).
    - SNOW: float, snow water equivalent (mm).
    - CUMDPT: float, profile depth to bottom (mm) from initialization.
    - DSMID: list[float], depth midpoints per layer (mm) from initialization.
    - TDL: float, accumulated DUL sum (cm), in/out.
    - TMA: list[float], 5-day memory vector, in/out.
    - NDays: int, wet-day memory counter (<=30), in/out.
    - WetDay: list[int], length 30, wet-day flags, in/out.
    - X2_PREV: float, previous X2_AVG (C), in/out.
    - SRFTEMP: float, previous surface temperature (C), in/out.
    - ST: list[float], previous soil temperature by layer (C), in/out.

    Returns:
    - TDL: float, updated (note: accumulates per original code when ISWWAT == 'N').
    - TMA: list[float], updated 5-day memory.
    - NDays: int, updated wet-day count.
    - WetDay: list[int], updated wet-day flags (length 30).
    - X2_PREV: float, updated X2_AVG (C).
    - SRFTEMP: float, updated surface temperature (C).
    - ST: list[float], updated soil temperatures (C), length NLAYR.
    - X2_AVG: float, updated 5-day average X2 (C).
    """
    import math

    # Layer aggregates
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL not reset per original code; it accumulates if ISWWAT == 'N'
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

    if ISWWAT.upper() == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Update 30-day wet-day memory
    WetDay = list(WetDay)
    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
    else:
        NDays += 1
    WetDay_today = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    WetDay[(NDays - 1) if NDays < 30 else 29] = WetDay_today
    NWetDays = sum(WetDay[:NDays])
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover factor including biomass and mulch and snow
    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV + math.exp(5.3396 - 2.3951 * CV) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if SNOW + math.exp(2.303 - 0.2197 * SNOW) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Update soil temperatures
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
        WetDay=WetDay_today,
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, X2_AVG


def test_stemp_epic_basic() -> tuple[
    float,  # CUMDPT
    list[float],  # DSMID
    float,  # TDL_end
    list[float],  # ST_end
    float,  # SRFTEMP_end
]:
    """
    Basic test exercising initialization and a few daily steps using
    representative inputs derived from ASKEE_DSSAT_EPIC.for.

    Returns key outputs to allow external assertions without executing here.
    """
    # Soil setup (4 layers)
    NLAYR = 4
    BD = [1.6] * NLAYR
    DLAYR = [10.0] * NLAYR  # cm
    DS = [10.0, 20.0, 30.0, 40.0]  # cm
    DUL = [0.30] * NLAYR
    LL = [0.20] * NLAYR
    SW = [0.25] * NLAYR

    # Weather and cover
    ISWWAT = "Y"
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0
    RAIN = 0.0
    DEPIR = 0.0
    BIOMAS = 0.0
    MULCHMASS = 0.0
    SNOW = 0.0

    # Initialize
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, X2_AVG = stemp_epic_initialize(
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
        X2_PREV=0.0,
    )

    # Run a few days
    for _ in range(5):
        TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, X2_AVG = stemp_epic_process(
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

    return CUMDPT, DSMID, TDL, ST, SRFTEMP