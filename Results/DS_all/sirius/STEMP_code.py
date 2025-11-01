from typing import List, Tuple
import math

def _truncate(value: float, decimals: int) -> float:
    """
    Truncate a float to a fixed number of decimal places (toward zero), like C# cast to int then divide.
    """
    factor = 10.0 ** decimals
    return int(value * factor) / factor


def SOILT(
    NL: int,
    ALBEDO: float,
    B: float,
    CUMDPT: float,
    DOY: int,
    DP: float,
    HDAY: float,
    NLAYR: int,
    PESW: float,
    SRAD: float,
    TAMP: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    WW: float,
    DSMID: List[float],
    ATOT: float,
    TMA: List[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    Compute soil temperatures across layers.

    Inputs:
    - NL: int, number of soil layers (array length)
    - ALBEDO: float, soil albedo (not used in function, kept for signature parity)
    - B: float, soil damping coefficient term (ln(500/DP))
    - CUMDPT: float, cumulative soil depth (cm)
    - DOY: int, day of year
    - DP: float, soil damping depth (mm? as per model)
    - HDAY: float, reference day for annual temp cycle (20 or 200 depending on hemisphere)
    - NLAYR: int, actual number of soil layers used
    - PESW: float, plant extractable soil water
    - SRAD: float, solar radiation (not used in function, kept for signature parity)
    - TAMP: float, amplitude of annual temperature function (degC)
    - TAV: float, average annual temperature (degC)
    - TAVG: float, average daily temperature (degC)
    - TMAX: float, max daily temperature (degC, not used here)
    - WW: float, water content term
    - DSMID: List[float], depth to midpoint of each soil layer (cm)
    - ATOT: float, sum of last 5 days average soil temperatures
    - TMA: List[float], last 5 days average soil temperatures

    Returns:
    - ATOT: float, updated sum of last 5 days average soil temperatures
    - TMA: List[float], updated last 5 days average soil temperatures
    - SRFTEMP: float, surface temperature (degC)
    - ST: List[float], updated soil temperature profile (degC)
    """
    ALX = (float(DOY) - HDAY) * 0.01740
    # Update last 5 days average soil temperatures
    ATOT = ATOT - TMA[5 - 1]
    for K in range(5 - 1, 0, -1):
        TMA[K] = TMA[K - 1]
    TMA[0] = TAVG
    TMA[0] = _truncate(TMA[0], 4)
    ATOT = ATOT + TMA[0]

    WC = max(0.010, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + (TAMP * math.cos(ALX) / 2.0)
    DT = ATOT / 5.0 - TA

    ST = [0.0] * NL
    for L in range(NLAYR):
        ZD = -(DSMID[L] / DD)
        ST[L] = TAV + ((TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD))
        ST[L] = _truncate(ST[L], 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)
    return ATOT, TMA, SRFTEMP, ST


def Init(
    NL: int,
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    MSALB: float,
    SW: List[float],
    XLAT: float,
    SRAD: float,
    TAVG: float,
    TMAX: float,
    TAV: float,
    TAMP: float,
    DOY: int,
) -> Tuple[
    float, List[float], float, List[float], float, float, List[float], float
]:
    """
    Initialize state variables for the STEMP soil temperature model.

    Inputs:
    - NL: int, number of soil layers (array length)
    - ISWWAT: str, water simulation control switch ("Y" or other)
    - BD: List[float], bulk density per layer (g/cm3)
    - DLAYR: List[float], thickness per layer (cm) [unused in Init; kept for parity]
    - DS: List[float], cumulative depth per layer (cm)
    - DUL: List[float], water content at drained upper limit per layer (cm3/cm3)
    - LL: List[float], lower limit water content per layer (cm3/cm3)
    - NLAYR: int, number of layers actually used
    - MSALB: float, soil albedo
    - SW: List[float], volumetric soil water per layer (cm3/cm3)
    - XLAT: float, latitude (deg) used to set HDAY reference
    - SRAD: float, solar radiation (MJ/m2-d)
    - TAVG: float, average daily temperature (degC)
    - TMAX: float, maximum daily temperature (degC)
    - TAV: float, annual average temperature (degC)
    - TAMP: float, amplitude of annual temperature (degC)
    - DOY: int, day of year

    Returns:
    - CUMDPT: float, cumulative profile depth (cm)
    - DSMID: List[float], depth to midpoint of each layer (cm)
    - TDL: float, total water content at drained upper limit (cm)
    - TMA: List[float], last 5 days average soil temperature (degC)
    - ATOT: float, sum of TMA (degC)
    - SRFTEMP: float, surface temperature (degC)
    - ST: List[float], soil temperature per layer (degC)
    - HDAY: float, reference day in the annual cycle (day)
    """
    CUMDPT = 0.0
    DSMID = [0.0] * NL
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NL
    HDAY = 0.0

    SWI = list(SW)
    DSI = list(DS)

    # Hemisphere-based reference day in the annual cycle
    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    DLI = [0.0] * NL
    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DSI[L]
        else:
            DLI[L] = DSI[L] - DSI[L - 1]
        DSMID[L] = CUMDPT + (DLI[L] * 5.0)
        CUMDPT = CUMDPT + (DLI[L] * 10.0)
        TBD = TBD + (BD[L] * DLI[L])
        TLL = TLL + (LL[L] * DLI[L])
        TSW = TSW + (SWI[L] * DLI[L])
        TDL = TDL + (DUL[L] * DLI[L])

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.630 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.3560 - (0.1440 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    for I in range(5):
        TMA[I] = _truncate(TAVG, 4)
    ATOT = TMA[0] * 5.0

    for L in range(NLAYR):
        ST[L] = TAVG

    # Spin-up iterations
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def CalculateModel(
    NL: int,
    ISWWAT: str,
    BD: List[float],
    DLAYR: List[float],
    DS: List[float],
    DUL: List[float],
    LL: List[float],
    NLAYR: int,
    MSALB: float,
    SW: List[float],
    SRAD: float,
    TAVG: float,
    TMAX: float,
    TAV: float,
    TAMP: float,
    DOY: int,
    CUMDPT: float,
    DSMID: List[float],
    TDL: float,
    TMA: List[float],
    ATOT: float,
    SRFTEMP: float,
    ST: List[float],
    HDAY: float,
) -> Tuple[
    float, List[float], float, List[float], float, float, List[float], float
]:
    """
    Main daily soil temperature update.

    Inputs:
    - NL: int, number of soil layers (array length)
    - ISWWAT: str, water simulation control switch ("Y" or other)
    - BD: List[float], bulk density per layer (g/cm3)
    - DLAYR: List[float], thickness per layer (cm)
    - DS: List[float], cumulative depth per layer (cm)
    - DUL: List[float], water content at drained upper limit per layer (cm3/cm3)
    - LL: List[float], lower limit water content per layer (cm3/cm3)
    - NLAYR: int, number of layers actually used
    - MSALB: float, soil albedo
    - SW: List[float], volumetric soil water per layer (cm3/cm3)
    - SRAD: float, solar radiation (MJ/m2-d)
    - TAVG: float, average daily temperature (degC)
    - TMAX: float, maximum daily temperature (degC)
    - TAV: float, annual average temperature (degC)
    - TAMP: float, amplitude of annual temperature (degC)
    - DOY: int, day of year
    - CUMDPT: float, cumulative profile depth (cm) [state]
    - DSMID: List[float], depth to midpoint per layer (cm) [state]
    - TDL: float, total water content at drained upper limit (cm) [state]
    - TMA: List[float], last 5 days average soil temperature (degC) [state]
    - ATOT: float, sum of TMA (degC) [state]
    - SRFTEMP: float, surface temperature (degC) [state]
    - ST: List[float], soil temperature per layer (degC) [state]
    - HDAY: float, reference day in the annual cycle (day) [state]

    Returns:
    - CUMDPT: float, unchanged cumulative depth (cm)
    - DSMID: List[float], unchanged midpoint depths (cm)
    - TDL: float, updated total water content at drained upper limit (cm)
    - TMA: List[float], updated last 5 days average soil temperature (degC)
    - ATOT: float, updated sum of TMA (degC)
    - SRFTEMP: float, updated surface temperature (degC)
    - ST: List[float], updated soil temperature profile (degC)
    - HDAY: float, unchanged reference day (day)
    """
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0

    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.630 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.3560 - (0.1440 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT, TMA, SRFTEMP, ST = SOILT(
        NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
    )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY