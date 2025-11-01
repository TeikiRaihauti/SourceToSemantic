from typing import List, Sequence, Tuple
import math


def soilt_epic(
    B: float,
    BCV: float,
    CUMDPT: float,
    DP: float,
    DSMID: Sequence[float],
    NLAYR: int,
    PESW: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    TMIN: float,
    WetDay: int,
    WFT: float,
    WW: float,
    TMA: Sequence[float],
    ST: Sequence[float],
    X2_PREV: float,
) -> Tuple[List[float], float, List[float], float]:
    """
    EPIC soil temperature routine for daily update by layer.

    Inputs:
    - B: float
    - BCV: float (soil cover factor)
    - CUMDPT: float (mm)
    - DP: float (mm)
    - DSMID: Sequence[float] (mm mid-depth per layer), length >= NLAYR
    - NLAYR: int
    - PESW: float (cm plant extractable soil water or proxy)
    - TAV: float (long-term average air temperature, deg C)
    - TAVG: float (daily average air temperature, deg C)
    - TMAX: float (deg C)
    - TMIN: float (deg C)
    - WetDay: int (0/1: wet day flag)
    - WFT: float (fraction of wet days in recent memory)
    - WW: float (volumetric water fraction, dimensionless)
    - TMA: Sequence[float], length 5 (temperature memory array)
    - ST: Sequence[float], soil temperature by layer (deg C), length >= NLAYR
    - X2_PREV: float (previous X2_AVG)

    Returns:
    - TMA_out: List[float], updated TMA length 5
    - SRFTEMP: float, updated surface temperature
    - ST_out: List[float], updated soil temperature by layer (deg C), length NLAYR
    - X2_PREV_out: float, updated previous X2_AVG
    """
    # Local copies for mutation
    TMA_out = list(TMA)
    ST_out = list(ST)

    # Ensure denominators are safe
    denom = WW * CUMDPT
    if denom <= 0.0:
        denom = 1.0

    WC = max(0.01, PESW) / denom * 10.0  # ratio

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # X2 depends on wet/dry status
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # EPIC memory update (note order matches Fortran: assign, then shift)
    TMA_out[0] = X2
    for k in range(4, 0, -1):  # 4..1
        TMA_out[k] = TMA_out[k - 1]

    X2_AVG = sum(TMA_out) / 5.0

    # Surface cover blending
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)

    # Profile update
    X1 = TAV - X3
    LAG = 0.5
    for l in range(NLAYR):
        ZD = DSMID[l] / DD if DD > 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST_out[l] = LAG * ST_out[l] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_out = X2_AVG

    return TMA_out, SRFTEMP, ST_out, X2_PREV_out


def stemp_epic_initialize(
    ISWWAT: str,
    BD: Sequence[float],
    DLAYR: Sequence[float],
    DS: Sequence[float],
    DUL: Sequence[float],
    LL: Sequence[float],
    NLAYR: int,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    MULCHMASS: float,
    SNOW: float,
) -> Tuple[float, List[float], float, List[float], int, List[int], float, float, List[float]]:
    """
    Seasonal initialization for EPIC soil temperature.

    Inputs:
    - ISWWAT: str ('Y' if soil water simulated)
    - BD: Sequence[float], bulk density by layer (g/cm3), length >= NLAYR
    - DLAYR: Sequence[float], layer thickness (cm), length >= NLAYR
    - DS: Sequence[float], cumulative depth (cm), length >= NLAYR (DS[NLAYR-1] total depth)
    - DUL: Sequence[float], drained upper limit by layer (cm3/cm3), length >= NLAYR
    - LL: Sequence[float], lower limit by layer (cm3/cm3), length >= NLAYR
    - NLAYR: int, number of soil layers
    - SW: Sequence[float], initial soil water by layer (cm3/cm3), length >= NLAYR
    - TAVG: float, daily average air temperature (deg C)
    - TMAX: float, daily max air temperature (deg C)
    - TMIN: float, daily min air temperature (deg C)
    - TAV: float, long-term average air temperature (deg C)
    - MULCHMASS: float, surface mulch mass (kg/ha)
    - SNOW: float, snow water equivalent (mm)

    Returns:
    - CUMDPT: float (mm), total profile depth
    - DSMID: List[float] (mm), depth to mid-point of each layer
    - TDL: float, sum over layers of DUL*DLAYR (cm)
    - TMA: List[float] length 5, initialized temperature memory
    - NDays: int, number of wet-day records stored (initially 0)
    - WetDay: List[int] length 30, wet-day memory (0/1)
    - X2_PREV: float, initialized previous X2_AVG
    - SRFTEMP: float, initialized surface temperature (deg C)
    - ST: List[float] length NLAYR, initialized soil temperature by layer (deg C)
    """
    # Initialize accumulators
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    DSMID: List[float] = [0.0] * NLAYR
    for l in range(NLAYR):
        DSMID[l] = CUMDPT + DLAYR[l] * 5.0  # mm to midpoint
        CUMDPT += DLAYR[l] * 10.0  # mm cumulative depth
        TBD += BD[l] * DLAYR[l]
        TLL += LL[l] * DLAYR[l]
        TSW += SW[l] * DLAYR[l]
        TDL += DUL[l] * DLAYR[l]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    total_depth_cm = DS[NLAYR - 1] if NLAYR > 0 else 1.0
    ABD = TBD / total_depth_cm if total_depth_cm > 0.0 else 0.0
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if ABD + 686.0 * math.exp(-5.63 * ABD) != 0 else 0.0
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP) if DP > 0.0 else 0.0

    # Initialize TMA with TAVG rounded similar to Fortran NINT(TAVG*10000)/10000
    tavg_rounded = round(TAVG, 4)
    TMA: List[float] = [tavg_rounded for _ in range(5)]
    X2_AVG = TMA[0] * 5.0  # as in Fortran

    ST: List[float] = [TAVG for _ in range(NLAYR)]

    # Wet-day memory
    WFT = 0.1
    WetDay: List[int] = [0] * 30
    NDays = 0

    # Soil cover function
    CV = (MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    # Spin-up calls (8 iterations)
    X2_PREV = 0.0
    SRFTEMP = TAVG
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_PREV = soilt_epic(
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


def stemp_epic_rate(
    ISWWAT: str,
    BD: Sequence[float],
    DLAYR: Sequence[float],
    DS: Sequence[float],
    DUL: Sequence[float],
    LL: Sequence[float],
    NLAYR: int,
    SW: Sequence[float],
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
    DSMID: Sequence[float],
    TDL: float,
    TMA: Sequence[float],
    NDays: int,
    WetDay: Sequence[int],
    X2_PREV: float,
    SRFTEMP: float,
    ST: Sequence[float],
) -> Tuple[float, List[float], int, List[int], float, float, List[float]]:
    """
    Daily rate calculations for EPIC soil temperature.

    Inputs:
    - ISWWAT: str ('Y' if soil water simulated)
    - BD, DLAYR, DS, DUL, LL: sequences by layer, length >= NLAYR
    - NLAYR: int
    - SW: Sequence[float], soil water by layer (cm3/cm3)
    - TAVG, TMAX, TMIN, TAV: floats (deg C)
    - RAIN: float (mm)
    - DEPIR: float (mm) (irrigation depth equivalent)
    - BIOMAS: float (kg/ha)
    - MULCHMASS: float (kg/ha)
    - SNOW: float (mm)
    - CUMDPT: float (mm), cumulative profile depth
    - DSMID: Sequence[float] (mm), depth to mid-point per layer
    - TDL: float, running sum used when ISWWAT != 'Y'
    - TMA: Sequence[float], length 5 (temperature memory)
    - NDays: int, number of wet-day records currently stored (<= 30)
    - WetDay: Sequence[int], length 30, wet-day memory (0/1)
    - X2_PREV: float
    - SRFTEMP: float
    - ST: Sequence[float], soil temperature by layer (deg C)

    Returns (updated):
    - TDL_out: float
    - TMA_out: List[float], length 5
    - NDays_out: int
    - WetDay_out: List[int], length 30
    - X2_PREV_out: float
    - SRFTEMP_out: float
    - ST_out: List[float], length NLAYR
    """
    # Local accumulators
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL_out = TDL
    for l in range(NLAYR):
        TBD += BD[l] * DLAYR[l]
        TDL_out += DUL[l] * DLAYR[l]
        TLL += LL[l] * DLAYR[l]
        TSW += SW[l] * DLAYR[l]

    total_depth_cm = DS[NLAYR - 1] if NLAYR > 0 else 1.0
    ABD = TBD / total_depth_cm if total_depth_cm > 0.0 else 0.0
    denom_fx = ABD + 686.0 * math.exp(-5.63 * ABD)
    FX = ABD / denom_fx if denom_fx != 0.0 else 0.0
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP) if DP > 0.0 else 0.0

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL_out - TLL)  # cm

    # Wet-day memory update (30-day window)
    WetDay_out = list(WetDay) if len(WetDay) == 30 else ([0] * 30)
    if NDays == 30:
        for i in range(29):
            WetDay_out[i] = WetDay_out[i + 1]
        # NDays remains 30
    else:
        NDays += 1  # increase count up to 30

    WetDay_out[NDays - 1] = 1 if (RAIN + DEPIR) > 1e-6 else 0
    NWetDays = sum(WetDay_out)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function
    CV = (BIOMAS + MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Call SOILT_EPIC
    TMA_out, SRFTEMP_out, ST_out, X2_PREV_out = soilt_epic(
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
        WetDay=WetDay_out[NDays - 1],
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return TDL_out, TMA_out, NDays, WetDay_out, X2_PREV_out, SRFTEMP_out, ST_out


def test_stemp_epic_basic():
    """
    Basic test derived from ASKEE_DSSAT_EPIC.for setup.
    Does not assert specific numeric values; checks shapes and execution.
    """
    # Soil setup
    NLAYR = 4
    BD = [1.6, 1.6, 1.6, 1.6]
    DLAYR = [10.0, 10.0, 10.0, 10.0]  # cm
    DS = [10.0, 20.0, 30.0, 40.0]  # cm cumulative depth
    DUL = [0.3, 0.3, 0.3, 0.3]
    LL = [0.2, 0.2, 0.2, 0.2]
    SW = [0.2, 0.2, 0.2, 0.2]

    # Weather and management
    ISWWAT = 'Y'
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
    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = stemp_epic_initialize(
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
    )

    assert isinstance(CUMDPT, float)
    assert len(DSMID) == NLAYR
    assert len(TMA) == 5
    assert len(WetDay) == 30
    assert len(ST) == NLAYR

    # One daily rate step
    TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = stemp_epic_rate(
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

    assert isinstance(SRFTEMP, float)
    assert len(TMA) == 5
    assert len(WetDay) == 30
    assert len(ST) == NLAYR