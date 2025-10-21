from typing import List, Tuple


def SOILT_EPIC(
    B: float,
    BCV: float,
    CUMDPT: float,
    DP: float,
    DSMID: List[float],
    NLAYR: int,
    PESW: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    TMIN: float,
    WetDay: int,
    WFT: float,
    WW: float,
    TMA: List[float],
    ST: List[float],
    X2_PREV: float,
) -> Tuple[List[float], float, List[float], float, float]:
    """
    SOILT_EPIC
    Determines soil temperature by layer (EPIC method).

    Inputs:
    - B: float
    - BCV: float
    - CUMDPT: float
    - DP: float
    - DSMID: List[float] (mm)
    - NLAYR: int
    - PESW: float (cm)
    - TAV: float
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - WetDay: int (0/1)
    - WFT: float
    - WW: float
    - TMA: List[float] (length 5), in/out
    - ST: List[float] (length NLAYR), in/out
    - X2_PREV: float, in/out

    Returns:
    - TMA: updated List[float] (length 5)
    - SRFTEMP: float
    - ST: updated List[float] (length NLAYR)
    - X2_AVG: float
    - X2_PREV: updated float
    """
    # Ensure copies for purity
    TMA_local = [float(v) for v in TMA]
    ST_local = [float(v) for v in ST]

    LAG = 0.5

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = pow(2.718281828459045, B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day temperature memory (preserving original order/behavior)
    TMA_local[0] = X2
    for K in range(4, 0, -1):
        TMA_local[K] = TMA_local[K - 1]
    X2_AVG = sum(TMA_local) / 5.0

    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)
    X1 = TAV - X3

    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + pow(2.718281828459045, -0.8669 - 2.0775 * ZD))
        ST_local[L] = LAG * ST_local[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_out = X2_AVG

    return TMA_local, SRFTEMP, ST_local, X2_AVG, X2_PREV_out


def STEMP_EPIC_initialize(
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    SW: List[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    MULCHMASS: float,
    SNOW: float,
) -> Tuple[
    float, List[float], float, List[float], int, List[int], float, float, List[float]
]:
    """
    Seasonal initialization for EPIC soil temperature.

    Inputs:
    - ISWWAT: str ('Y' if water simulated)
    - BD: List[float], bulk density by layer
    - DLAYR: List[float], layer thickness (cm)
    - DS: List[float], cumulative depth to bottom of each layer (cm)
    - DUL: List[float], drained upper limit by layer
    - LL: List[float], lower limit by layer
    - NLAYR: int, number of layers
    - SW: List[float], soil water by layer
    - TAVG: float, daily average air temperature
    - TMAX: float, daily max air temperature
    - TMIN: float, daily min air temperature
    - TAV: float, annual average air temperature
    - MULCHMASS: float, surface mulch mass (kg/ha)
    - SNOW: float, snow water equivalent (mm)

    Returns:
    - CUMDPT: float, cumulative profile depth (mm)
    - DSMID: List[float], depth to mid-point of each layer (mm)
    - TDL: float, sum(DUL*DLAYR)
    - TMA: List[float] length 5, temperature memory
    - NDays: int, initialized to 0
    - WetDay: List[int] length 30, initialized to zeros
    - X2_PREV: float
    - SRFTEMP: float
    - ST: List[float], soil temperature by layer
    """
    SWI = [float(v) for v in SW]

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0 for _ in range(NLAYR)]
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0
        CUMDPT = CUMDPT + DLAYR[L] * 10.0
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * pow(2.718281828459045, -5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = 1.0
    if DP != 0.0:
        B = 1.6094379124341003 - (1.0) * (DP / 500.0)  # placeholder
    # Fortran uses natural log: B = ALOG(500.0/DP)
    # Implement directly to avoid precision issues:
    import math

    B = math.log(500.0 / DP)

    TMA = [round(TAVG, 4) for _ in range(5)]
    X2_AVG = TMA[0] * 5.0
    ST = [TAVG for _ in range(NLAYR)]

    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + pow(2.718281828459045, 5.3396 - 2.3951 * CV)) if (CV + pow(2.718281828459045, 5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + pow(2.718281828459045, 2.303 - 0.2197 * SNOW)) if (SNOW + pow(2.718281828459045, 2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG
    X2_PREV = 0.0
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
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    SW: List[float],
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
    DSMID: List[float],
    TDL: float,
    TMA: List[float],
    NDays: int,
    WetDay: List[int],
    X2_PREV: float,
    ST: List[float],
) -> Tuple[float, List[float], int, List[int], float, float, List[float], float, float, float]:
    """
    Daily soil temperature update (EPIC method).

    Inputs:
    - ISWWAT: str ('Y' if water simulated)
    - BD, DLAYR, DS, DUL, LL: List[float] by layer
    - NLAYR: int
    - SW: List[float] soil water by layer
    - TAVG, TMAX, TMIN, TAV: float temperatures
    - RAIN: float (mm)
    - DEPIR: float (mm), irrigation depth
    - BIOMAS: float (kg/ha)
    - MULCHMASS: float (kg/ha)
    - SNOW: float (mm)
    - CUMDPT: float (mm), from initialization
    - DSMID: List[float] (mm), from initialization
    - TDL: float, in/out (accumulates per original code)
    - TMA: List[float] length 5, in/out
    - NDays: int, in/out
    - WetDay: List[int] length 30, in/out
    - X2_PREV: float, in/out
    - ST: List[float] soil temperature by layer, in/out

    Returns (updated):
    - TDL: float
    - TMA: List[float]
    - NDays: int
    - WetDay: List[int]
    - X2_PREV: float
    - SRFTEMP: float
    - ST: List[float]
    - ABD: float
    - DP: float
    - WW: float
    """
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # preserves original behavior
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * pow(2.718281828459045, -5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    import math

    B = math.log(500.0 / DP)

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
    else:
        NDays += 1
    WetDay_today = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    WetDay_index = min(NDays, 30) - 1
    if WetDay_index >= 0:
        WetDay[WetDay_index] = WetDay_today
    NWetDays = sum(WetDay[:NDays])
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + pow(2.718281828459045, 5.3396 - 2.3951 * CV)) if (CV + pow(2.718281828459045, 5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + pow(2.718281828459045, 2.303 - 0.2197 * SNOW)) if (SNOW + pow(2.718281828459045, 2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    TMA_out, SRFTEMP, ST_out, _, X2_PREV_out = SOILT_EPIC(
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

    return TDL, TMA_out, NDays, WetDay, X2_PREV_out, SRFTEMP, ST_out, ABD, DP, WW