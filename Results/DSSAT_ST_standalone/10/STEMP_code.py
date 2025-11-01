from typing import List, Tuple
import math


def YR_DOY(YRDOY: int) -> Tuple[int, int]:
    """
    Convert YRDOY to YEAR and DOY.

    Inputs:
    - YRDOY: int

    Returns:
    - YEAR: int
    - DOY: int
    """
    YEAR = YRDOY // 1000
    DOY = YRDOY - YEAR * 1000
    return YEAR, DOY


def _nint(x: float) -> int:
    """
    Fortran-like NINT: nearest integer with halves away from zero.

    Inputs:
    - x: float

    Returns:
    - int
    """
    if x >= 0.0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


def _round_nint_decimals(x: float, decimals: int) -> float:
    """
    Round to given decimals using Fortran NINT behavior.

    Inputs:
    - x: float
    - decimals: int

    Returns:
    - float
    """
    factor = 10.0 ** decimals
    return _nint(x * factor) / factor


def SOILT(
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
    SRFTEMP: float,
    ST: List[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    Determines soil temperature by layer.

    Inputs:
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
    - DSMID: List[float] (length >= NLAYR)
    - ATOT: float
    - TMA: List[float] (length 5)
    - SRFTEMP: float
    - ST: List[float] (length >= NLAYR)

    Returns (updated):
    - ATOT: float
    - TMA: List[float] (length 5)
    - SRFTEMP: float
    - ST: List[float] (length NLAYR)
    """
    # Copy arrays to avoid mutating inputs
    TMA_out = list(TMA)
    ST_out = list(ST)

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT_out = ATOT - TMA_out[4]

    # shift moving average window
    for K in range(4, 0, -1):
        TMA_out[K] = TMA_out[K - 1]

    TMA_out[0] = TAVG
    TMA_out[0] = _round_nint_decimals(TMA_out[0], 4)
    ATOT_out = ATOT_out + TMA_out[0]

    # Corrected water content function (EPIC correction)
    denom = WW * CUMDPT
    if denom == 0.0:
        WC = 0.01  # avoid division by zero; consistent minimal value
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_out / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD if DD != 0.0 else 0.0
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_out[L] = _round_nint_decimals(ST_val, 3)

    SRFTEMP_out = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_out, TMA_out, SRFTEMP_out, ST_out


def STEMP(
    CONTROL_DYNAMIC: int,
    CONTROL_YRDOY: int,
    CONTROL_RUN: int,
    CONTROL_RNMODE: str,
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: List[float],
    SOILPROP_DLAYR: List[float],
    SOILPROP_DS: List[float],
    SOILPROP_DUL: List[float],
    SOILPROP_LL: List[float],
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
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
) -> Tuple[float, List[float], float, List[float], float, float, List[float]]:
    """
    Determines soil temperature by layer across dynamic phases.

    Inputs:
    - CONTROL_DYNAMIC: int (1=RUNINIT, 2=SEASINIT, 3=RATE, 5=OUTPUT, 6=SEASEND)
    - CONTROL_YRDOY: int
    - CONTROL_RUN: int
    - CONTROL_RNMODE: str
    - ISWITCH_ISWWAT: str ('Y'/'N')
    - SOILPROP_BD: List[float] (length >= NLAYR)
    - SOILPROP_DLAYR: List[float] (length >= NLAYR)
    - SOILPROP_DS: List[float] (length >= NLAYR)
    - SOILPROP_DUL: List[float] (length >= NLAYR)
    - SOILPROP_LL: List[float] (length >= NLAYR)
    - SOILPROP_NLAYR: int
    - SOILPROP_MSALB: float
    - SRAD: float
    - SW: List[float] (length >= NLAYR)
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - CUMDPT: float (InOut)
    - DSMID: List[float] (InOut, length >= NLAYR)
    - TDL: float (InOut)
    - TMA: List[float] (InOut, length 5)
    - ATOT: float (InOut)
    - SRFTEMP: float (InOut)
    - ST: List[float] (InOut, length >= NLAYR)

    Returns (updated):
    - CUMDPT: float
    - DSMID: List[float] (length NLAYR)
    - TDL: float
    - TMA: List[float] (length 5)
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float] (length NLAYR)
    """
    # Local copies for outputs
    CUMDPT_out = CUMDPT
    DSMID_out = list(DSMID)
    TDL_out = TDL
    TMA_out = list(TMA)
    ATOT_out = ATOT
    SRFTEMP_out = SRFTEMP
    ST_out = list(ST)

    # Dynamic phase codes
    RUNINIT = 1
    SEASINIT = 2
    RATE = 3
    OUTPUT = 5
    SEASEND = 6

    # Date conversion
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Hemispheric hottest day
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    if CONTROL_DYNAMIC == SEASINIT:
        # Initialize for season only when required
        should_initialize = (CONTROL_RUN == 1) or (CONTROL_RNMODE not in ("Q", "F"))
        if should_initialize:
            SWI = list(SW)
            DSI = list(SOILPROP_DS)

            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_out = 0.0
            CUMDPT_out = 0.0

            for L in range(SOILPROP_NLAYR):
                if L == 0:
                    DLI = DSI[L]
                else:
                    DLI = DSI[L] - DSI[L - 1]

                DSMID_out[L] = CUMDPT_out + DLI * 5.0
                CUMDPT_out = CUMDPT_out + DLI * 10.0
                TBD += SOILPROP_BD[L] * DLI
                TLL += SOILPROP_LL[L] * DLI
                TSW += SWI[L] * DLI
                TDL_out += SOILPROP_DUL[L] * DLI

            if ISWITCH_ISWWAT == "Y":
                PESW = max(0.0, TSW - TLL)
            else:
                PESW = max(0.0, TDL_out - TLL)

            ABD = TBD / DSI[SOILPROP_NLAYR - 1]
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP)
            ALBEDO = SOILPROP_MSALB

            # Initialize moving average temperatures and profile temps
            init_tavg = _round_nint_decimals(TAVG, 4)
            TMA_out = [init_tavg] * 5
            ATOT_out = TMA_out[0] * 5.0

            for L in range(SOILPROP_NLAYR):
                ST_out[L] = TAVG

            # Spin-up
            for _ in range(8):
                ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
                    ALBEDO,
                    B,
                    CUMDPT_out,
                    DOY,
                    DP,
                    HDAY,
                    SOILPROP_NLAYR,
                    PESW,
                    SRAD,
                    TAMP,
                    TAV,
                    TAVG,
                    TMAX,
                    WW,
                    DSMID_out,
                    ATOT_out,
                    TMA_out,
                    SRFTEMP_out,
                    ST_out,
                )

    elif CONTROL_DYNAMIC == RATE:
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        for L in range(SOILPROP_NLAYR):
            TBD += SOILPROP_BD[L] * SOILPROP_DLAYR[L]
            TDL_out += SOILPROP_DUL[L] * SOILPROP_DLAYR[L]  # Note: not reset (matches original behavior)
            TLL += SOILPROP_LL[L] * SOILPROP_DLAYR[L]
            TSW += SW[L] * SOILPROP_DLAYR[L]

        ABD = TBD / SOILPROP_DS[SOILPROP_NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = SOILPROP_MSALB

        if ISWITCH_ISWWAT == "Y":
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_out - TLL)

        ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
            ALBEDO,
            B,
            CUMDPT_out,
            DOY,
            DP,
            HDAY,
            SOILPROP_NLAYR,
            PESW,
            SRAD,
            TAMP,
            TAV,
            TAVG,
            TMAX,
            WW,
            DSMID_out,
            ATOT_out,
            TMA_out,
            SRFTEMP_out,
            ST_out,
        )

    elif CONTROL_DYNAMIC in (OUTPUT, SEASEND):
        # No computation; output was I/O in original and is omitted here.
        pass

    return CUMDPT_out, DSMID_out, TDL_out, TMA_out, ATOT_out, SRFTEMP_out, ST_out