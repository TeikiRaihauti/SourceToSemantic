from typing import List, Sequence, Tuple
import math

def _nint(value: float) -> int:
    """
    Fortran-like NINT: nearest integer, halves away from zero.
    """
    if value >= 0:
        return int(math.floor(value + 0.5))
    else:
        return int(math.ceil(value - 0.5))


def _round_nint(value: float, decimals: int) -> float:
    """
    Fortran-like rounding via NINT at specified decimals.
    """
    factor = 10.0 ** decimals
    return _nint(value * factor) / factor


def YR_DOY(YRDOY: int) -> Tuple[int, int]:
    """
    Compute (YEAR, DOY) from YRDOY (e.g., 2024123 -> (2024, 123)).

    Inputs:
    - YRDOY: int

    Returns:
    - YEAR: int
    - DOY: int
    """
    YEAR = YRDOY // 1000
    DOY = YRDOY - YEAR * 1000
    return YEAR, DOY


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
    DSMID: Sequence[float],
    ATOT: float,
    TMA: Sequence[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    SOILT: Determines soil temperature by layer.

    Inputs:
    - ALBEDO: float (not used in algorithm; kept for interface parity)
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float (not used in algorithm; kept for interface parity)
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float (not used in algorithm; kept for interface parity)
    - WW: float
    - DSMID: Sequence[float], length >= NLAYR
    - ATOT: float (InOut)
    - TMA: Sequence[float], length 5 (InOut)

    Returns:
    - ATOT: float (updated)
    - TMA: List[float], length 5 (updated)
    - SRFTEMP: float
    - ST: List[float], soil temperature by layer, length NLAYR
    """
    # Copy and enforce TMA length 5
    TMA_local = list(TMA[:5]) + [0.0] * max(0, 5 - len(TMA))
    # Seasonal phase angle
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average queue (5 days)
    ATOT = ATOT - TMA_local[4]
    for K in range(4, 0, -1):
        TMA_local[K] = TMA_local[K - 1]
    TMA_local[0] = TAVG
    TMA_local[0] = _round_nint(TMA_local[0] * 1.0, 4)
    ATOT = ATOT + TMA_local[0]

    # Water content function (corrected)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    ST = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[L] = _round_nint(ST_val * 1.0, 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA_local, SRFTEMP, ST


def STEMP(
    CONTROL_DYNAMIC: int,
    CONTROL_YRDOY: int,
    CONTROL_RUN: int,
    CONTROL_RNMODE: str,
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: Sequence[float],
    SOILPROP_DLAYR: Sequence[float],
    SOILPROP_DS: Sequence[float],
    SOILPROP_DUL: Sequence[float],
    SOILPROP_LL: Sequence[float],
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
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
) -> Tuple[float, List[float], float, List[float], float, List[float], List[float]]:
    """
    STEMP: Determines soil temperature by layer. Handles seasonal init and daily rate.

    Inputs:
    - CONTROL_DYNAMIC: int (2=SEASINIT, 3=RATE; other values ignored here)
    - CONTROL_YRDOY: int
    - CONTROL_RUN: int
    - CONTROL_RNMODE: str
    - ISWITCH_ISWWAT: str ('Y' or 'N')
    - SOILPROP_BD: Sequence[float], bulk density by layer
    - SOILPROP_DLAYR: Sequence[float], layer thickness (mm)
    - SOILPROP_DS: Sequence[float], depth to bottom of layer (cm)
    - SOILPROP_DUL: Sequence[float], drained upper limit (volumetric)
    - SOILPROP_LL: Sequence[float], lower limit (volumetric)
    - SOILPROP_NLAYR: int, number of soil layers
    - SOILPROP_MSALB: float, soil/mulch albedo
    - SRAD: float, daily solar radiation (MJ/m2-d)
    - SW: Sequence[float], soil water content (volumetric) by layer
    - TAVG: float, average air temperature (C)
    - TMAX: float, max air temperature (C)
    - XLAT: float, latitude (deg)
    - TAV: float, annual average air temperature (C)
    - TAMP: float, annual amplitude of air temperature (C)
    - CUMDPT: float (InOut), cumulative profile depth (mm)
    - DSMID: Sequence[float] (InOut), depth to mid-point of each layer (mm)
    - TDL: float (InOut), sum of DUL*thickness across profile
    - TMA: Sequence[float] (InOut), 5-day moving average queue [most recent first]
    - ATOT: float (InOut), sum of TMA elements
    - SRFTEMP: float (InOut), surface temperature (C)
    - ST: Sequence[float] (InOut), soil temperature by layer (C)

    Returns:
    - CUMDPT: float (possibly updated at SEASINIT)
    - DSMID: List[float] (possibly updated at SEASINIT)
    - TDL: float (updated in RATE for ISWWAT='N' behavior per original)
    - TMA: List[float] (updated)
    - ATOT: float (updated)
    - SRFTEMP: float (updated)
    - ST: List[float] (updated)
    """
    # Local copies to ensure purity
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SW_local = list(SW[:NLAYR])

    DSMID_local = list(DSMID[:NLAYR]) + [0.0] * max(0, NLAYR - len(DSMID))
    TMA_local = list(TMA[:5]) + [0.0] * max(0, 5 - len(TMA))
    ST_local = list(ST[:NLAYR]) + [0.0] * max(0, NLAYR - len(ST))
    CUMDPT_local = float(CUMDPT)
    TDL_local = float(TDL)
    ATOT_local = float(ATOT)
    SRFTEMP_local = float(SRFTEMP)

    # Compute DOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Hottest day (hemisphere)
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    if CONTROL_DYNAMIC == 2:  # SEASINIT
        RUN = int(CONTROL_RUN)
        RNMODE = str(CONTROL_RNMODE)

        if (RUN == 1) or (RNMODE not in ('Q', 'F')):
            # Initialize from initial soil water (SW) and DS
            SWI = SW_local
            DSI = DS

            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_local = 0.0
            CUMDPT_local = 0.0

            for L in range(NLAYR):
                if L == 0:
                    DLI = DSI[L]
                else:
                    DLI = DSI[L] - DSI[L - 1]
                DSMID_local[L] = CUMDPT_local + DLI * 5.0
                CUMDPT_local = CUMDPT_local + DLI * 10.0
                TBD += BD[L] * DLI
                TLL += LL[L] * DLI
                TSW += SWI[L] * DLI
                TDL_local += DUL[L] * DLI

            if ISWITCH_ISWWAT == 'Y':
                PESW = max(0.0, TSW - TLL)
            else:
                PESW = max(0.0, TDL_local - TLL)

            ABD = TBD / DSI[NLAYR - 1]
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP)
            ALBEDO = MSALB

            for i in range(5):
                TMA_local[i] = _round_nint(TAVG * 1.0, 4)
            ATOT_local = TMA_local[0] * 5.0

            for L in range(NLAYR):
                ST_local[L] = TAVG

            for _ in range(8):
                ATOT_local, TMA_local, SRFTEMP_local, ST_local = SOILT(
                    ALBEDO, B, CUMDPT_local, DOY, DP, HDAY, NLAYR,
                    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_local,
                    ATOT_local, TMA_local
                )

    elif CONTROL_DYNAMIC == 3:  # RATE
        TBD = 0.0
        TLL = 0.0
        TSW_sum = 0.0
        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_local = TDL_local + DUL[L] * DLAYR[L]
            TLL += LL[L] * DLAYR[L]
            TSW_sum += SW_local[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = MSALB

        if ISWITCH_ISWWAT == 'Y':
            PESW = max(0.0, TSW_sum - TLL)
        else:
            PESW = max(0.0, TDL_local - TLL)

        ATOT_local, TMA_local, SRFTEMP_local, ST_local = SOILT(
            ALBEDO, B, CUMDPT_local, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_local,
            ATOT_local, TMA_local
        )

    # For OUTPUT/SEASEND dynamic steps, original code only performs I/O which is intentionally omitted.

    return (
        CUMDPT_local,
        DSMID_local,
        TDL_local,
        TMA_local,
        ATOT_local,
        SRFTEMP_local,
        ST_local,
    )


def test_STEMP_basic() -> None:
    """
    Basic test derived from ASKEE.for-like setup. Does not assert values; ensures no errors.
    """
    # Control and switches
    CONTROL_DYNAMIC = 2  # SEASINIT
    CONTROL_YRDOY = 2021100
    CONTROL_RUN = 1
    CONTROL_RNMODE = 'B'
    ISWITCH_ISWWAT = 'Y'

    # Soil properties (4 layers)
    SOILPROP_BD = [1.6, 1.5, 1.4, 1.3]
    SOILPROP_DLAYR = [100.0, 100.0, 100.0, 100.0]
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.30, 0.30, 0.30, 0.30]
    SOILPROP_LL = [0.20, 0.20, 0.20, 0.20]
    SOILPROP_NLAYR = 4
    SOILPROP_MSALB = 0.13

    # Weather and water
    SRAD = 20.0
    SW = [0.20, 0.20, 0.20, 0.20]
    TAVG = 25.0
    TMAX = 30.0
    XLAT = 28.0
    TAV = 20.0
    TAMP = 10.0

    # State variables
    CUMDPT = 0.0
    DSMID = [0.0] * SOILPROP_NLAYR
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * SOILPROP_NLAYR

    # Seasonal init
    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        CONTROL_DYNAMIC, CONTROL_YRDOY, CONTROL_RUN, CONTROL_RNMODE,
        ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
        SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
        CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    )

    # Daily rate
    CONTROL_DYNAMIC = 3
    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        CONTROL_DYNAMIC, CONTROL_YRDOY, CONTROL_RUN, CONTROL_RNMODE,
        ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
        SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
        CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    )
    # No assertions; this function is for manual inspection or further extension.