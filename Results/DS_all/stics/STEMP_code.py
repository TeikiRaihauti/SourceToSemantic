from typing import List, Sequence, Tuple
import math


def _trunc(x: float, decimals: int) -> float:
    """Truncate toward zero to a fixed number of decimals (Fortran INT behavior)."""
    factor = 10.0 ** decimals
    return math.trunc(x * factor) / factor


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
    DSMID: Sequence[float],
    ATOT: float,
    TMA: Sequence[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    Compute soil temperature profile and surface temperature.

    Inputs:
    - NL: int
    - ALBEDO: float
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float
    - WW: float
    - DSMID: Sequence[float] (len >= NLAYR)
    - ATOT: float
    - TMA: Sequence[float] (len == 5)

    Returns:
    - ATOT: float (updated)
    - TMA: List[float] (len 5, updated)
    - SRFTEMP: float
    - ST: List[float] (len NL)
    """
    # Copy state arrays to avoid mutating inputs
    TMA_out = list(TMA)
    ST_out = [0.0] * NL

    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average of soil temperatures (last 5 days)
    ATOT_out = ATOT - TMA_out[4]
    for k in range(4, 0, -1):
        TMA_out[k] = TMA_out[k - 1]
    TMA_out[0] = _trunc(TAVG, 4)
    ATOT_out = ATOT_out + TMA_out[0]

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * (((1.0 - WC) / (1.0 + WC)) ** 2))
    DD = FX * DP

    TA = TAV + (TAMP * math.cos(ALX) / 2.0)
    DT = ATOT_out / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        st_val = TAV + ((TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD))
        ST_out[L] = _trunc(st_val, 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_out, TMA_out, SRFTEMP, ST_out


def init_stemp(
    NL: int,
    ISWWAT: str,
    BD: Sequence[float],
    DLAYR: Sequence[float],
    DS: Sequence[float],
    DUL: Sequence[float],
    LL: Sequence[float],
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    DOY: int,
) -> Tuple[float, List[float], float, List[float], float, float, List[float], float]:
    """
    Initialize STEMP state.

    Inputs:
    - NL: int
    - ISWWAT: str
    - BD: Sequence[float] (len >= NLAYR)
    - DLAYR: Sequence[float] (len >= NLAYR)
    - DS: Sequence[float] (len >= NLAYR)
    - DUL: Sequence[float] (len >= NLAYR)
    - LL: Sequence[float] (len >= NLAYR)
    - NLAYR: int
    - MSALB: float
    - SRAD: float
    - SW: Sequence[float] (len >= NLAYR)
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - DOY: int

    Returns:
    - CUMDPT: float
    - DSMID: List[float] (len NL)
    - TDL: float
    - TMA: List[float] (len 5)
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float] (len NL)
    - HDAY: float
    """
    DSMID = [0.0] * NL
    ST = [0.0] * NL
    TMA = [0.0] * 5

    # Hemisphere-based "harvest day" seed
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    # Layer thickness from cumulative depths and midpoints (units: mm)
    for L in range(NLAYR):
        if L == 0:
            DLI = DS[L]
        else:
            DLI = DS[L] - DS[L - 1]
        DSMID[L] = CUMDPT + (DLI * 5.0)  # midpoint depth (mm)
        CUMDPT = CUMDPT + (DLI * 10.0)   # cumulative depth (mm)

        TBD += BD[L] * DLI
        TLL += LL[L] * DLI
        TSW += SW[L] * DLI
        TDL += DUL[L] * DLI

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp((-5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # not directly used in SOILT but kept for signature consistency

    # Initialize last-5-day average soil temperatures and layer temps
    for i in range(5):
        TMA[i] = _trunc(TAVG, 4)
    ATOT = TMA[0] * 5.0
    for L in range(NLAYR):
        ST[L] = TAVG

    SRFTEMP = 0.0
    # Spin-up calls (8 cycles) to establish initial ST profile
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def model_stemp(
    NL: int,
    ISWWAT: str,
    BD: Sequence[float],
    DLAYR: Sequence[float],
    DS: Sequence[float],
    DUL: Sequence[float],
    LL: Sequence[float],
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CUMDPT: float,
    DSMID: Sequence[float],
    TDL: float,
    TMA: Sequence[float],
    ATOT: float,
    SRFTEMP: float,
    ST: Sequence[float],
    DOY: int,
    HDAY: float,
) -> Tuple[float, List[float], float, List[float], float, float, List[float]]:
    """
    Daily STEMP process step.

    Inputs:
    - NL: int
    - ISWWAT: str
    - BD: Sequence[float] (len >= NLAYR)
    - DLAYR: Sequence[float] (len >= NLAYR)
    - DS: Sequence[float] (len >= NLAYR)
    - DUL: Sequence[float] (len >= NLAYR)
    - LL: Sequence[float] (len >= NLAYR)
    - NLAYR: int
    - MSALB: float
    - SRAD: float
    - SW: Sequence[float] (len >= NLAYR)
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - CUMDPT: float
    - DSMID: Sequence[float] (len >= NLAYR)
    - TDL: float
    - TMA: Sequence[float] (len == 5)
    - ATOT: float
    - SRFTEMP: float
    - ST: Sequence[float] (len >= NLAYR)
    - DOY: int
    - HDAY: float

    Returns:
    - CUMDPT: float
    - DSMID: List[float] (len NL)
    - TDL: float (updated as in Fortran)
    - TMA: List[float] (len 5)
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float] (len NL)
    """
    # Copy state arrays to avoid mutating inputs
    DSMID_out = list(DSMID) + [0.0] * max(0, NL - len(DSMID))
    ST_in = list(ST) + [0.0] * max(0, NL - len(ST))
    TMA_in = list(TMA)

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL is intentionally not reset here to preserve original Fortran behavior.
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp((-5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
        NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_out, ATOT, TMA_in
    )

    return CUMDPT, DSMID_out, TDL, TMA_out, ATOT_out, SRFTEMP_out, ST_out