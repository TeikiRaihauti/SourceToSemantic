from typing import List, Tuple
import math

def _truncate(value: float, decimals: int) -> float:
    """
    Truncate a float to a given number of decimal places without rounding.
    """
    if decimals < 0:
        return value
    factor = 10 ** decimals
    return math.trunc(value * factor) / factor


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
    Compute soil temperature profile and surface temperature for one time step.

    Inputs:
    - NL: int, Number of soil layers
    - ALBEDO: float, Soil albedo (unused in algorithm, retained for API parity)
    - B: float, Empirical coefficient derived from damping depth
    - CUMDPT: float, Cumulative depth of soil profile
    - DOY: int, Day of year
    - DP: float, Damping depth (mm)
    - HDAY: float, Reference day based on latitude
    - NLAYR: int, Actual number of soil layers
    - PESW: float, Plant extractable soil water in profile (cm)
    - SRAD: float, Solar radiation (unused in algorithm, retained for API parity)
    - TAMP: float, Annual temperature amplitude (degC)
    - TAV: float, Annual average temperature (degC)
    - TAVG: float, Average daily air temperature (degC)
    - TMAX: float, Max daily air temperature (unused in algorithm, retained for API parity)
    - WW: float, Parameter dependent on bulk density
    - DSMID: List[float], Depth to midpoint of soil layers
    - ATOT: float, Sum of last 5 days soil temps
    - TMA: List[float], Last 5 days average soil temps

    Returns:
    - ATOT: float, Updated sum of last 5 days soil temps
    - TMA: List[float], Updated last 5 days average soil temps
    - SRFTEMP: float, Surface (litter) temperature (degC)
    - ST: List[float], Updated soil temperature by layer (degC)
    """
    ST = [0.0] * NL
    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT = ATOT - TMA[5 - 1]
    for K in range(5, 1, -1):
        TMA[K - 1] = TMA[K - 2]
    TMA[0] = TAVG
    TMA[0] = _truncate(TMA[0], 4)
    ATOT = ATOT + TMA[0]
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP
    TA = TAV + (TAMP * math.cos(ALX) / 2.0)
    DT = ATOT / 5.0 - TA
    for L in range(1, NLAYR + 1):
        ZD = -(DSMID[L - 1] / DD)
        ST[L - 1] = TAV + ((TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD))
        ST[L - 1] = _truncate(ST[L - 1], 3)
    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)
    return ATOT, TMA, SRFTEMP, ST


def Init(
    NL: int,
    ISWWAT: str,
    BD: List[float],
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
) -> Tuple[float, List[float], float, List[float], float, float, List[float], float]:
    """
    Initialization of the STEMP component.

    Inputs:
    - NL: int, Number of soil layers
    - ISWWAT: str, Water simulation control switch ("Y" or "N")
    - BD: List[float], Bulk density per layer (g/cm3)
    - DS: List[float], Cumulative depth per layer (cm)
    - DUL: List[float], Volumetric water at drained upper limit per layer
    - LL: List[float], Volumetric water at lower limit per layer
    - NLAYR: int, Actual number of soil layers
    - MSALB: float, Soil albedo with mulch and water effects
    - SW: List[float], Current volumetric soil water per layer
    - XLAT: float, Latitude (deg)
    - SRAD: float, Solar radiation (MJ/m2-d)
    - TAVG: float, Average daily air temperature (degC)
    - TMAX: float, Maximum daily air temperature (degC)
    - TAV: float, Average annual soil temperature (degC)
    - TAMP: float, Amplitude of annual soil temperature (degC)
    - DOY: int, Day of year

    Returns (state variables):
    - CUMDPT: float, Cumulative depth of soil profile
    - DSMID: List[float], Depth to midpoint of soil layer L
    - TDL: float, Total water at drained upper limit (cm)
    - TMA: List[float], Previous 5 days average soil temperatures
    - ATOT: float, Sum of TMA
    - SRFTEMP: float, Surface (litter) temperature (degC)
    - ST: List[float], Soil temperature in soil layer L (degC)
    - HDAY: float, Harvest/reference day used in temperature seasonal phase
    """
    CUMDPT = 0.0
    DSMID = [0.0] * NL
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NL
    HDAY = 0.0

    I = 0
    L = 0
    ABD = 0.0
    ALBEDO = 0.0
    B = 0.0
    DP = 0.0
    FX = 0.0
    PESW = 0.0
    TBD = 0.0
    WW = 0.0
    TLL = 0.0
    TSW = 0.0
    DLI = [0.0] * NL
    DSI = DS[:]
    SWI = SW[:]

    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    for L in range(1, NLAYR + 1):
        if L == 1:
            DLI[L - 1] = DSI[L - 1]
        else:
            DLI[L - 1] = DSI[L - 1] - DSI[L - 2]
        DSMID[L - 1] = CUMDPT + (DLI[L - 1] * 5.0)
        CUMDPT = CUMDPT + (DLI[L - 1] * 10.0)
        TBD = TBD + (BD[L - 1] * DLI[L - 1])
        TLL = TLL + (LL[L - 1] * DLI[L - 1])
        TSW = TSW + (SWI[L - 1] * DLI[L - 1])
        TDL = TDL + (DUL[L - 1] * DLI[L - 1])

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    for I in range(1, 5 + 1):
        TMA[I - 1] = _truncate(TAVG, 4)

    ATOT = TMA[0] * 5.0

    for L in range(1, NLAYR + 1):
        ST[L - 1] = TAVG

    for I in range(1, 8 + 1):
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
) -> Tuple[float, List[float], float, List[float], float, float, List[float]]:
    """
    Main biophysical process function of the STEMP component (one time step).

    Inputs:
    - NL: int, Number of soil layers
    - ISWWAT: str, Water simulation control switch ("Y" or "N")
    - BD: List[float], Bulk density per layer (g/cm3)
    - DLAYR: List[float], Thickness of soil layer L (cm)
    - DS: List[float], Cumulative depth per layer (cm)
    - DUL: List[float], Volumetric water at drained upper limit per layer
    - LL: List[float], Volumetric water at lower limit per layer
    - NLAYR: int, Actual number of soil layers
    - MSALB: float, Soil albedo with mulch and water effects
    - SW: List[float], Current volumetric soil water per layer
    - SRAD: float, Solar radiation (MJ/m2-d)
    - TAVG: float, Average daily air temperature (degC)
    - TMAX: float, Maximum daily air temperature (degC)
    - TAV: float, Average annual soil temperature (degC)
    - TAMP: float, Amplitude of annual soil temperature (degC)
    - DOY: int, Day of year
    - CUMDPT: float, State - Cumulative depth of soil profile
    - DSMID: List[float], State - Depth to midpoint of soil layers
    - TDL: float, State - Total water at drained upper limit (cm); updated cumulatively
    - TMA: List[float], State - Last 5 days average soil temps
    - ATOT: float, State - Sum of last 5 days soil temps
    - SRFTEMP: float, State - Surface (litter) temperature (degC)
    - ST: List[float], State - Soil temperature per layer (degC)
    - HDAY: float, State - Reference day based on latitude

    Returns (updated state variables):
    - CUMDPT: float
    - DSMID: List[float]
    - TDL: float
    - TMA: List[float]
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float]
    """
    L = 0
    ABD = 0.0
    ALBEDO = 0.0
    B = 0.0
    DP = 0.0
    FX = 0.0
    PESW = 0.0
    TBD = 0.0
    WW = 0.0
    TLL = 0.0
    TSW = 0.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0

    for L in range(1, NLAYR + 1):
        TBD = TBD + (BD[L - 1] * DLAYR[L - 1])
        TDL = TDL + (DUL[L - 1] * DLAYR[L - 1])
        TLL = TLL + (LL[L - 1] * DLAYR[L - 1])
        TSW = TSW + (SW[L - 1] * DLAYR[L - 1])

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT, TMA, SRFTEMP, ST = SOILT(
        NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
    )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST