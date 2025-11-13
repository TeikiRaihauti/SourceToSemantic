from typing import List, Tuple


def _truncate(value: float, decimals: int) -> float:
    # Truncate a float to a fixed number of decimal places (Java-style truncation toward zero).
    factor = 10 ** decimals
    return int(value * factor) / factor


def Calculate_SOILT(
    t_NL: int,
    ALBEDO: float,
    B: float,
    t_CUMDPT: float,
    t_DOY: int,
    DP: float,
    t_HDAY: float,
    t_NLAYR: int,
    PESW: float,
    t_SRAD: float,
    t_TAMP: float,
    t_TAV: float,
    t_TAVG: float,
    t_TMAX: float,
    WW: float,
    t_DSMID: List[float],
    t_ATOT: float,
    t_TMA: List[float],
    t_ST: List[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    Calculate soil temperature (SOILT).

    Inputs:
    - t_NL: int
    - ALBEDO: float
    - B: float
    - t_CUMDPT: float
    - t_DOY: int
    - DP: float
    - t_HDAY: float
    - t_NLAYR: int
    - PESW: float
    - t_SRAD: float
    - t_TAMP: float
    - t_TAV: float
    - t_TAVG: float
    - t_TMAX: float
    - WW: float
    - t_DSMID: List[float]
    - t_ATOT: float
    - t_TMA: List[float] (length 5)
    - t_ST: List[float] (length >= t_NLAYR)

    Returns:
    - t_ATOT: float
    - t_TMA: List[float]
    - t_SRFTEMP: float
    - t_ST: List[float]
    """
    ALX = (float(t_DOY) - t_HDAY) * 0.0174
    t_ATOT = t_ATOT - t_TMA[5 - 1]
    for K in range(5, 1, -1):
        t_TMA[K - 1] = t_TMA[K - 1 - 1]
    t_TMA[1 - 1] = (1.0 - ALBEDO) * (t_TAVG + ((t_TMAX - t_TAVG) * (t_SRAD * 0.03) ** 0.5)) + (ALBEDO * t_TMA[(1 - 1)])
    t_TMA[1 - 1] = _truncate(t_TMA[(1 - 1)], 4)
    t_ATOT = t_ATOT + t_TMA[1 - 1]
    WC = max(0.01, PESW) / (WW * t_CUMDPT) * 10.0
    FX = pow(2.718281828459045, B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP
    TA = t_TAV + (t_TAMP * __import__("math").cos(ALX) / 2.0)
    DT = t_ATOT / 5.0 - TA
    for L in range(1, t_NLAYR + 1):
        ZD = -(t_DSMID[(L - 1)] / DD)
        t_ST[L - 1] = t_TAV + (((t_TAMP / 2.0) * __import__("math").cos((ALX + ZD)) + DT) * __import__("math").exp(ZD))
        t_ST[L - 1] = _truncate(t_ST[(L - 1)], 3)
    t_SRFTEMP = t_TAV + (t_TAMP / 2.0 * __import__("math").cos(ALX) + DT)
    return t_ATOT, t_TMA, t_SRFTEMP, t_ST


def init(
    NL: int,
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: List[float],
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    DOY: int,
) -> Tuple[float, List[float], float, List[float], float, float, List[float], float]:
    """
    Initialization for STEMP.

    Inputs:
    - NL: int (Number of soil layers)
    - ISWWAT: str (Water simulation control switch, e.g., "Y"/"N")
    - BD: List[float] (Bulk density per layer)
    - DLAYR: List[float] (Layer thickness, cm)
    - DS: List[float] (Cumulative depth in soil layer L, cm)
    - DUL: List[float] (Volumetric water content at Drained Upper Limit per layer)
    - LL: List[float] (Volumetric water content at lower limit per layer)
    - NLAYR: int (Actual number of soil layers)
    - MSALB: float (Soil albedo with mulch and water effects)
    - SRAD: float (Solar radiation)
    - SW: List[float] (Volumetric soil water content per layer)
    - TAVG: float (Average daily temperature)
    - TMAX: float (Maximum daily temperature)
    - XLAT: float (Latitude)
    - TAV: float (Average annual soil temperature)
    - TAMP: float (Amplitude for soil temperature)
    - DOY: int (Day of year)

    Returns:
    - CUMDPT: float
    - DSMID: List[float]
    - TDL: float
    - TMA: List[float] (last 5 days soil temperature)
    - ATOT: float (sum of TMA)
    - SRFTEMP: float
    - ST: List[float] (soil temperature per layer)
    - HDAY: float
    """
    t_NL = NL
    t_ISWWAT = ISWWAT
    t_BD = BD
    t_DLAYR = DLAYR
    t_DS = DS
    t_DUL = DUL
    t_LL = LL
    t_NLAYR = NLAYR
    t_MSALB = MSALB
    t_SRAD = SRAD
    t_SW = SW
    t_TAVG = TAVG
    t_TMAX = TMAX
    t_XLAT = XLAT
    t_TAV = TAV
    t_TAMP = TAMP
    t_DOY = DOY

    t_CUMDPT = 0.0
    t_DSMID = [0.0] * t_NL
    t_TDL = 0.0
    t_TMA = [0.0] * 5
    t_ATOT = 0.0
    t_SRFTEMP = 0.0
    t_ST = [0.0] * t_NL
    t_HDAY = 0.0

    DLI = [0.0] * t_NL
    DSI = t_DS
    SWI = t_SW

    t_HDAY = 20.0 if t_XLAT < 0.0 else 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    t_TDL = 0.0
    t_CUMDPT = 0.0
    for L in range(1, t_NLAYR + 1):
        if L == 1:
            DLI[L - 1] = DSI[L - 1]
        else:
            DLI[L - 1] = DSI[L - 1] - DSI[L - 1 - 1]
        t_DSMID[L - 1] = t_CUMDPT + (DLI[(L - 1)] * 5.0)
        t_CUMDPT = t_CUMDPT + (DLI[(L - 1)] * 10.0)
        TBD = TBD + (t_BD[(L - 1)] * DLI[(L - 1)])
        TLL = TLL + (t_LL[(L - 1)] * DLI[(L - 1)])
        TSW = TSW + (SWI[(L - 1)] * DLI[(L - 1)])
        t_TDL = t_TDL + (t_DUL[(L - 1)] * DLI[(L - 1)])

    if t_ISWWAT.strip().upper() == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, t_TDL - TLL)

    ABD = TBD / DSI[(t_NLAYR - 1)]
    FX = ABD / (ABD + (686.0 * __import__("math").exp(-(5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = __import__("math").log(500.0 / DP)
    ALBEDO = t_MSALB

    for I in range(1, 5 + 1):
        t_TMA[I - 1] = _truncate(t_TAVG, 4)
    t_ATOT = t_TMA[(1 - 1)] * 5.0
    for L in range(1, t_NLAYR + 1):
        t_ST[L - 1] = t_TAVG

    for _ in range(8):
        t_ATOT, t_TMA, t_SRFTEMP, t_ST = Calculate_SOILT(
            t_NL, ALBEDO, B, t_CUMDPT, t_DOY, DP, t_HDAY, t_NLAYR, PESW, t_SRAD, t_TAMP, t_TAV, t_TAVG, t_TMAX, WW,
            t_DSMID, t_ATOT, t_TMA, t_ST
        )

    return t_CUMDPT, t_DSMID, t_TDL, t_TMA, t_ATOT, t_SRFTEMP, t_ST, t_HDAY


def process(
    NL: int,
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: List[float],
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CUMDPT: float,
    DSMID: List[float],
    TDL: float,
    TMA: List[float],
    ATOT: float,
    SRFTEMP: float,
    ST: List[float],
    DOY: int,
    HDAY: float,
) -> Tuple[float, List[float], float, List[float], float, float, List[float]]:
    """
    Main daily process for STEMP.

    Inputs:
    - NL: int
    - ISWWAT: str
    - BD: List[float]
    - DLAYR: List[float]
    - DS: List[float]
    - DUL: List[float]
    - LL: List[float]
    - NLAYR: int
    - MSALB: float
    - SRAD: float
    - SW: List[float]
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - CUMDPT: float
    - DSMID: List[float]
    - TDL: float
    - TMA: List[float]
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float]
    - DOY: int
    - HDAY: float

    Returns:
    - CUMDPT: float
    - DSMID: List[float]
    - TDL: float
    - TMA: List[float]
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float]
    """
    t_NL = NL
    t_ISWWAT = ISWWAT
    t_BD = BD
    t_DLAYR = DLAYR
    t_DS = DS
    t_DUL = DUL
    t_LL = LL
    t_NLAYR = NLAYR
    t_MSALB = MSALB
    t_SRAD = SRAD
    t_SW = SW
    t_TAVG = TAVG
    t_TMAX = TMAX
    t_XLAT = XLAT
    t_TAV = TAV
    t_TAMP = TAMP
    t_CUMDPT = CUMDPT
    t_DSMID = DSMID
    t_TDL = TDL
    t_TMA = TMA
    t_ATOT = ATOT
    t_SRFTEMP = SRFTEMP
    t_ST = ST
    t_DOY = DOY
    t_HDAY = HDAY

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(1, t_NLAYR + 1):
        TBD = TBD + (t_BD[(L - 1)] * t_DLAYR[(L - 1)])
        t_TDL = t_TDL + (t_DUL[(L - 1)] * t_DLAYR[(L - 1)])
        TLL = TLL + (t_LL[(L - 1)] * t_DLAYR[(L - 1)])
        TSW = TSW + (t_SW[(L - 1)] * t_DLAYR[(L - 1)])
    ABD = TBD / t_DS[(t_NLAYR - 1)]
    FX = ABD / (ABD + (686.0 * __import__("math").exp(-(5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = __import__("math").log(500.0 / DP)
    ALBEDO = t_MSALB
    if t_ISWWAT.strip().upper() == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, t_TDL - TLL)

    t_ATOT, t_TMA, t_SRFTEMP, t_ST = Calculate_SOILT(
        t_NL, ALBEDO, B, t_CUMDPT, t_DOY, DP, t_HDAY, t_NLAYR, PESW, t_SRAD, t_TAMP, t_TAV, t_TAVG, t_TMAX, WW,
        t_DSMID, t_ATOT, t_TMA, t_ST
    )

    return t_CUMDPT, t_DSMID, t_TDL, t_TMA, t_ATOT, t_SRFTEMP, t_ST