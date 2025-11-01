from typing import List, Tuple


def _fortran_nint(x: float) -> int:
    # Fortran NINT: nearest integer, ties away from zero
    if x >= 0.0:
        return int((x + 0.5) // 1)
    else:
        return -int(((-x) + 0.5) // 1)


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
) -> Tuple[float, List[float], List[float], float, float]:
    """
    SOILT_EPIC
    Inputs:
      B: float
      BCV: float
      CUMDPT: float (mm)
      DP: float (mm)
      DSMID: list[float] (mm)
      NLAYR: int
      PESW: float (cm)
      TAV: float
      TAVG: float
      TMAX: float
      TMIN: float
      WetDay: int (0 or 1)
      WFT: float
      WW: float
      TMA: list[float] length 5 (InOut)
      ST: list[float] length NLAYR (InOut)
      X2_PREV: float (InOut)
    Returns:
      SRFTEMP: float
      ST: list[float]
      TMA: list[float]
      X2_AVG: float
      X2_PREV: float
    """
    import math

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Shift memory with EPIC convention
    TMA = TMA[:]  # copy
    TMA[0] = X2
    for K in range(4, 0, -1):
        TMA[K] = TMA[K - 1]
    X2_AVG = sum(TMA) / 5.0

    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)
    X1 = TAV - X3

    LAG = 0.5
    ST = ST[:]
    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return SRFTEMP, ST, TMA, X2_AVG, X2_PREV


def STEMP_EPIC_SEASINIT(
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
    X2_PREV: float,
) -> Tuple[float, List[float], float, List[float], int, List[int], float, float, List[float]]:
    """
    Seasonal initialization for EPIC soil temperature.

    Inputs:
      ISWWAT: str ('Y' or 'N')
      BD: list[float] bulk density by layer
      DLAYR: list[float] layer thickness (cm)
      DS: list[float] cumulative depth (cm) by layer
      DUL: list[float]
      LL: list[float]
      NLAYR: int
      SW: list[float] initial water content by layer
      TAVG: float
      TMAX: float
      TMIN: float
      TAV: float
      MULCHMASS: float (kg/ha)
      SNOW: float (mm)
      X2_PREV: float (initial previous X2 average)

    Returns:
      CUMDPT: float (mm)
      DSMID: list[float] (mm)
      TDL: float
      TMA: list[float] length 5
      NDays: int
      WetDay: list[int] length 30
      X2_PREV: float
      SRFTEMP: float
      ST: list[float] length NLAYR
    """
    import math

    SWI = SW[:]

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0
        CUMDPT = CUMDPT + DLAYR[L] * 10.0
        TBD += BD[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SWI[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    TMA = [float(_fortran_nint(TAVG * 10000.) / 10000.0) for _ in range(5)]
    X2_AVG = TMA[0] * 5.0

    ST = [TAVG for _ in range(NLAYR)]

    WFT = 0.1
    WetDay = [0] * 30
    NDays = 0

    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG
    for _ in range(8):
        SRFTEMP, ST, TMA, X2_AVG, X2_PREV = SOILT_EPIC(
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
    SRFTEMP: float,
    ST: List[float],
) -> Tuple[float, List[float], int, List[int], float, float, List[float], float]:
    """
    Daily EPIC soil temperature process.

    Inputs:
      ISWWAT: str ('Y' or 'N')
      BD: list[float]
      DLAYR: list[float] (cm)
      DS: list[float] (cm cumulative)
      DUL: list[float]
      LL: list[float]
      NLAYR: int
      SW: list[float]
      TAVG: float
      TMAX: float
      TMIN: float
      TAV: float
      RAIN: float (mm)
      DEPIR: float (mm)
      BIOMAS: float (kg/ha)
      MULCHMASS: float (kg/ha)
      SNOW: float (mm)
      CUMDPT: float (mm)
      DSMID: list[float] (mm)
      TDL: float (InOut)
      TMA: list[float] length 5 (InOut)
      NDays: int (InOut)
      WetDay: list[int] length 30 (InOut)
      X2_PREV: float (InOut)
      SRFTEMP: float (InOut)
      ST: list[float] length NLAYR (InOut)

    Returns:
      TDL: float
      TMA: list[float]
      NDays: int
      WetDay: list[int]
      X2_PREV: float
      SRFTEMP: float
      ST: list[float]
      X2_AVG: float
    """
    import math

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
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    if NDays == 30:
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
    else:
        NDays = NDays + 1

    WetDay_value = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    WetDay_index = min(NDays - 1, 29)
    WetDay[WetDay_index] = WetDay_value

    NWetDays = sum(WetDay[:NDays])
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP, ST, TMA, X2_AVG, X2_PREV = SOILT_EPIC(
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
        WetDay=WetDay_value,
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, X2_AVG


def test_STEMP_EPIC_smoke() -> Tuple[float, List[float], float, List[float]]:
    """
    Smoke test derived from ASKEE_DSSAT_EPIC.for example.
    Returns:
      SRFTEMP after daily step
      ST after daily step
      CUMDPT
      DSMID
    """
    ISWWAT = 'Y'
    NLAYR = 4
    BD = [1.6, 1.6, 1.6, 1.6]
    DLAYR = [10.0, 10.0, 10.0, 10.0]
    DS = [10.0, 20.0, 30.0, 40.0]
    DUL = [0.3, 0.3, 0.3, 0.3]
    LL = [0.2, 0.2, 0.2, 0.2]
    SW = [0.2, 0.2, 0.2, 0.2]

    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0

    MULCHMASS = 0.0
    SNOW = 0.0
    X2_PREV = 0.0

    CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST = STEMP_EPIC_SEASINIT(
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

    RAIN = 0.0
    DEPIR = 0.0
    BIOMAS = 0.0

    TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST, _ = STEMP_EPIC(
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

    return SRFTEMP, ST, CUMDPT, DSMID