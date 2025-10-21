from typing import List, Sequence, Tuple
import math


def SOILT_EPIC(
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
    SRFTEMP: float,
    ST: Sequence[float],
    X2_AVG: float,
    X2_PREV: float,
) -> Tuple[List[float], float, List[float], float, float]:
    """
    SOILT_EPIC
    Determines soil temperature by layer.

    Inputs:
    - B: float
    - BCV: float
    - CUMDPT: float
    - DP: float
    - DSMID: Sequence[float], length NLAYR
    - NLAYR: int
    - PESW: float
    - TAV: float
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - WetDay: int (0/1)
    - WFT: float
    - WW: float
    - TMA: Sequence[float], length 5
    - SRFTEMP: float
    - ST: Sequence[float], length NLAYR
    - X2_AVG: float
    - X2_PREV: float

    Returns:
    - TMA: List[float], updated length 5
    - SRFTEMP: float
    - ST: List[float], updated length NLAYR
    - X2_AVG: float
    - X2_PREV: float
    """
    # Copy inputs to avoid mutating caller data
    TMA_new = list(TMA)
    ST_new = list(ST)

    # Water content effect
    # WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    denom = WW * CUMDPT
    if denom <= 0.0:
        # Guard against divide-by-zero; preserve behavior by using a large WC
        WC = 1.0e6
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    # Surface temperature drivers (Potter & Williams, 1994)
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # EPIC memory of surface temperature (note: order preserved as in original code)
    TMA_new[0] = X2
    for K in range(4, 0, -1):
        TMA_new[K] = TMA_new[K - 1]
    X2_AVG_new = sum(TMA_new) / 5.0

    # Soil cover effects
    X3 = (1.0 - BCV) * X2_AVG_new + BCV * X2_PREV

    # Surface temperature
    SRFTEMP_new = min(X2_AVG_new, X3)

    # Profile temperature update
    X1 = TAV - X3
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        ST_new[L] = LAG * ST_new[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_new = X2_AVG_new

    return TMA_new, SRFTEMP_new, ST_new, X2_AVG_new, X2_PREV_new


def STEMP_EPIC_initialize(
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
    X2_PREV: float = 0.0,
) -> Tuple[float, List[float], float, List[float], int, List[int], float, float, List[float]]:
    """
    STEMP_EPIC seasonal initialization.

    Inputs:
    - ISWWAT: str ('Y' if water simulated)
    - BD: Sequence[float], length NLAYR
    - DLAYR: Sequence[float], length NLAYR (cm)
    - DS: Sequence[float], length NLAYR (cm to bottom of layer)
    - DUL: Sequence[float], length NLAYR
    - LL: Sequence[float], length NLAYR
    - NLAYR: int
    - SW: Sequence[float], length NLAYR (initial soil water)
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - TAV: float
    - MULCHMASS: float (kg/ha)
    - SNOW: float (mm)
    - X2_PREV: float (initial), default 0.0

    Returns:
    - CUMDPT: float (mm)
    - DSMID: List[float], length NLAYR (mm)
    - TDL: float
    - TMA: List[float], length 5
    - NDays: int
    - WetDay: List[int], length 30
    - X2_PREV: float
    - SRFTEMP: float
    - ST: List[float], length NLAYR
    """
    # Initialize cumulative depths and sums
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0

    SWI = list(SW)

    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm to midpoint of layer
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm cumulative depth
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT.upper() == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1] if DS[NLAYR - 1] != 0.0 else 0.0
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP) if DP != 0.0 else 0.0

    TMA = [round(TAVG, 4) for _ in range(5)]
    X2_AVG = TMA[0] * 5.0  # per original code prior to calling SOILT_EPIC

    ST = [TAVG for _ in range(NLAYR)]

    # Wet day memory
    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    # Soil cover function
    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

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
            SRFTEMP=SRFTEMP,
            ST=ST,
            X2_AVG=X2_AVG,
            X2_PREV=X2_PREV,
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
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
    STEMP_EPIC daily rate calculations.

    Inputs:
    - ISWWAT: str ('Y' if water simulated)
    - BD: Sequence[float], length NLAYR
    - DLAYR: Sequence[float], length NLAYR (cm)
    - DS: Sequence[float], length NLAYR (cm)
    - DUL: Sequence[float], length NLAYR
    - LL: Sequence[float], length NLAYR
    - NLAYR: int
    - SW: Sequence[float], length NLAYR
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - TAV: float
    - RAIN: float
    - DEPIR: float
    - BIOMAS: float
    - MULCHMASS: float
    - SNOW: float
    - CUMDPT: float
    - DSMID: Sequence[float], length NLAYR (mm)
    - TDL: float (in/out)
    - TMA: Sequence[float], length 5 (in/out)
    - NDays: int (in/out)
    - WetDay: Sequence[int], length 30 (in/out)
    - X2_PREV: float (in/out)
    - SRFTEMP: float (in/out)
    - ST: Sequence[float], length NLAYR (in/out)

    Returns:
    - TDL: float
    - TMA: List[float], updated length 5
    - NDays: int
    - WetDay: List[int], updated length 30
    - X2_PREV: float
    - SRFTEMP: float
    - ST: List[float], updated length NLAYR
    """
    TMA_new = list(TMA)
    WetDay_new = list(WetDay)
    ST_new = list(ST)

    # Aggregates across layers
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL is intentionally not reset here to preserve original behavior.
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1] if DS[NLAYR - 1] != 0.0 else 0.0
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP) if DP != 0.0 else 0.0

    if ISWWAT.upper() == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # 30-day wet day memory
    if NDays == 30:
        for i in range(29):
            WetDay_new[i] = WetDay_new[i + 1]
    else:
        NDays += 1
    if RAIN + DEPIR > 1.0e-6:
        WetDay_new[NDays - 1] = 1
    else:
        WetDay_new[NDays - 1] = 0
    NWetDays = sum(WetDay_new[:NDays])
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function
    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Daily soil temperature update
    X2_AVG_dummy = 0.0  # not returned to caller; internal to update
    TMA_new, SRFTEMP_new, ST_new, X2_AVG_dummy, X2_PREV_new = SOILT_EPIC(
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
        WetDay=WetDay_new[NDays - 1],
        WFT=WFT,
        WW=WW,
        TMA=TMA_new,
        SRFTEMP=SRFTEMP,
        ST=ST_new,
        X2_AVG=X2_AVG_dummy,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA_new, NDays, WetDay_new, X2_PREV_new, SRFTEMP_new, ST_new


def test_STEMP_EPIC_basic() -> None:
    """
    Basic test derived from ASKEE_DSSAT_EPIC.for program.

    Constructs a minimal setup, runs initialization and a few daily updates.
    Verifies shapes and basic execution.
    """
    # Inputs
    ISWWAT = 'Y'
    NLAYR = 4
    BD = [1.6] * NLAYR
    DLAYR = [10.0] * NLAYR  # cm
    DS = [10.0, 20.0, 30.0, 40.0]  # cm
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
        X2_PREV=0.0,
    )

    assert isinstance(CUMDPT, float)
    assert len(DSMID) == NLAYR
    assert len(TMA) == 5
    assert isinstance(NDays, int)
    assert len(WetDay) == 30
    assert isinstance(X2_PREV, float)
    assert isinstance(SRFTEMP, float)
    assert len(ST) == NLAYR

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

    # Basic post-conditions
    assert len(TMA) == 5
    assert len(WetDay) == 30
    assert len(ST) == NLAYR