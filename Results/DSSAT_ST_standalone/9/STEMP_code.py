from typing import List, Sequence, Tuple


def YR_DOY(YRDOY: int) -> Tuple[int, int]:
    """
    YR_DOY
    Converts YRDOY to YR and DOY.

    Inputs:
    - YRDOY: int

    Returns:
    - (YR, DOY): Tuple[int, int]
    """
    YR = YRDOY // 1000
    DOY = YRDOY - YR * 1000
    return YR, DOY


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
    SRFTEMP: float,
    ST: Sequence[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    SOILT
    Determines soil temperature by layer.

    Inputs:
    - ALBEDO: float
    - B: float
    - CUMDPT: float (mm)
    - DOY: int
    - DP: float (mm)
    - HDAY: float (DOY of hottest day; 200 N hemisphere, 20 S hemisphere)
    - NLAYR: int
    - PESW: float (cm)
    - SRAD: float (unused)
    - TAMP: float (deg C)
    - TAV: float (deg C)
    - TAVG: float (deg C)
    - TMAX: float (unused)
    - WW: float (dimensionless)
    - DSMID: Sequence[float] (mm), length >= NLAYR
    - ATOT: float (running total of TMA)
    - TMA: Sequence[float], length >= 5
    - SRFTEMP: float (deg C)
    - ST: Sequence[float], length >= NLAYR

    Returns:
    - ATOT: float
    - TMA: List[float] (length 5)
    - SRFTEMP: float
    - ST: List[float] (length NLAYR)
    """
    import math

    def nint(x: float) -> int:
        return int(math.floor(x + 0.5)) if x >= 0.0 else int(math.ceil(x - 0.5))

    # Defensive copies and sizing
    TMA_list = list(TMA[:5]) + [0.0] * max(0, 5 - len(TMA))
    ST_list = list(ST[:NLAYR]) + [0.0] * max(0, NLAYR - len(ST))
    DSMID_list = list(DSMID[:NLAYR]) + [0.0] * max(0, NLAYR - len(DSMID))

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT = ATOT - TMA_list[4]

    # Shift TMA
    for K in range(4, 0, -1):
        TMA_list[K] = TMA_list[K - 1]
    TMA_list[0] = TAVG

    # Keep only 4 decimals
    TMA_list[0] = nint(TMA_list[0] * 10000.0) / 10000.0
    ATOT = ATOT + TMA_list[0]

    # Corrected water content function
    denom = WW * CUMDPT
    WC = (max(0.01, PESW) / denom * 10.0) if denom != 0.0 else 0.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID_list[L] / DD if DD != 0.0 else 0.0
        val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        # Keep only 3 decimals
        ST_list[L] = nint(val * 1000.0) / 1000.0

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA_list, SRFTEMP, ST_list


def STEMP(
    CONTROL_DYNAMIC: int,
    CONTROL_YRDOY: int,
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
    CONTROL_RUN: int,
    CONTROL_RNMODE: str,
) -> Tuple[float, List[float], float, List[float], float, float, List[float]]:
    """
    STEMP
    Determines soil temperature by layer. Executes initialization or daily rate
    calculations depending on CONTROL_DYNAMIC value.

    Inputs:
    - CONTROL_DYNAMIC: int (2=SEASINIT, 3=RATE, 5=OUTPUT, 6=SEASEND)
    - CONTROL_YRDOY: int
    - ISWITCH_ISWWAT: str ('Y' or 'N')
    - SOILPROP_BD: Sequence[float], bulk density by layer
    - SOILPROP_DLAYR: Sequence[float], thickness by layer
    - SOILPROP_DS: Sequence[float], cumulative depth by layer (same units as input data; used as in original)
    - SOILPROP_DUL: Sequence[float], drained upper limit by layer
    - SOILPROP_LL: Sequence[float], lower limit by layer
    - SOILPROP_NLAYR: int, number of layers
    - SOILPROP_MSALB: float, mulch/soil albedo
    - SRAD: float
    - SW: Sequence[float], soil water by layer
    - TAVG: float, daily average soil surface temperature (deg C)
    - TMAX: float
    - XLAT: float, latitude (deg; sign determines hemisphere)
    - TAV: float, average annual temperature (deg C)
    - TAMP: float, annual temperature amplitude (deg C)
    - CUMDPT: float (InOut), cumulative profile depth (mm)
    - DSMID: Sequence[float] (InOut), depth to middle of layer (mm), length >= NLAYR
    - TDL: float (InOut), running sum (see original code behavior)
    - TMA: Sequence[float] (InOut), moving average memory, length >= 5
    - ATOT: float (InOut), running total of TMA entries
    - SRFTEMP: float (InOut), surface temperature (deg C)
    - ST: Sequence[float] (InOut), soil temperature by layer (deg C), length >= NLAYR
    - CONTROL_RUN: int
    - CONTROL_RNMODE: str (single character; typically 'B','Q','F', etc.)

    Returns (updated):
    - CUMDPT: float
    - DSMID: List[float] (length NLAYR)
    - TDL: float
    - TMA: List[float] (length 5)
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float] (length NLAYR)
    """
    import math

    # Defensive copies and sizing
    NLAYR = int(SOILPROP_NLAYR)
    BD = _expand_to_layers(SOILPROP_BD, NLAYR)
    DLAYR = _expand_to_layers(SOILPROP_DLAYR, NLAYR)
    DS = _expand_to_layers(SOILPROP_DS, NLAYR)
    DUL = _expand_to_layers(SOILPROP_DUL, NLAYR)
    LL = _expand_to_layers(SOILPROP_LL, NLAYR)
    SW_arr = _expand_to_layers(SW, NLAYR)

    DSMID_list = list(DSMID[:NLAYR]) + [0.0] * max(0, NLAYR - len(DSMID))
    ST_list = list(ST[:NLAYR]) + [0.0] * max(0, NLAYR - len(ST))
    TMA_list = list(TMA[:5]) + [0.0] * max(0, 5 - len(TMA))

    # Compute DOY from YRDOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Hemisphere hottest day
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    # Seasonal initialization
    if CONTROL_DYNAMIC == 2:
        RUN = CONTROL_RUN
        RNMODE = CONTROL_RNMODE

        # Only initialize for first run or when RNMODE not in {'Q','F'}
        if RUN == 1 or (RNMODE not in ('Q', 'F')):
            # Initial water and depth arrays
            SWI = list(SW_arr)
            DSI = list(DS)

            TBD = 0.0
            TLL = 0.0
            TSW_sum = 0.0
            TDL = 0.0
            CUMDPT = 0.0

            # Determine layer increments from cumulative depth and initialize profile metrics
            for L in range(NLAYR):
                if L == 0:
                    DLI = DSI[L]
                else:
                    DLI = DSI[L] - DSI[L - 1]

                DSMID_list[L] = CUMDPT + DLI * 5.0  # mm
                CUMDPT = CUMDPT + DLI * 10.0        # mm

                TBD += BD[L] * DLI
                TLL += LL[L] * DLI
                TSW_sum += SWI[L] * DLI
                TDL += DUL[L] * DLI

            if ISWITCH_ISWWAT == 'Y':
                PESW = max(0.0, TSW_sum - TLL)
            else:
                PESW = max(0.0, TDL - TLL)

            ABD = TBD / DSI[NLAYR - 1] if DSI[NLAYR - 1] != 0.0 else 0.0
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP) if DP != 0.0 else 0.0
            ALBEDO = SOILPROP_MSALB  # Not used in SOILT

            # Initialize TMA and ATOT
            TMA_list = [_nint4(TAVG)] * 5
            ATOT = TMA_list[0] * 5.0

            # Initialize ST by layer
            for L in range(NLAYR):
                ST_list[L] = TAVG

            # Warm-up iterations
            for _ in range(8):
                ATOT, TMA_list, SRFTEMP, ST_list = SOILT(
                    ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
                    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_list,
                    ATOT, TMA_list, SRFTEMP, ST_list
                )

    # Daily rate calculations
    elif CONTROL_DYNAMIC == 3:
        TBD = 0.0
        TLL = 0.0
        TSW_sum = 0.0
        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL += DUL[L] * DLAYR[L]
            TLL += LL[L] * DLAYR[L]
            TSW_sum += SW_arr[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1] if DS[NLAYR - 1] != 0.0 else 0.0
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP) if DP != 0.0 else 0.0
        ALBEDO = SOILPROP_MSALB  # Not used in SOILT

        if ISWITCH_ISWWAT == 'Y':
            PESW = max(0.0, TSW_sum - TLL)
        else:
            PESW = max(0.0, TDL - TLL)

        ATOT, TMA_list, SRFTEMP, ST_list = SOILT(
            ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_list,
            ATOT, TMA_list, SRFTEMP, ST_list
        )

    # OUTPUT or SEASEND: no state changes (I/O removed)
    elif CONTROL_DYNAMIC in (5, 6):
        pass

    return CUMDPT, DSMID_list, TDL, TMA_list, ATOT, SRFTEMP, ST_list


def _expand_to_layers(value: Sequence[float] or float, NLAYR: int) -> List[float]:
    """
    Expand scalar or sequence to list of length NLAYR (Fortran-style array assignment).
    """
    if isinstance(value, (int, float)):
        return [float(value)] * NLAYR
    lst = list(value)
    if len(lst) < NLAYR:
        lst = lst + [lst[-1] if lst else 0.0] * (NLAYR - len(lst))
    return lst[:NLAYR]


def _nint4(x: float) -> float:
    """
    Round to 4 decimals using Fortran-like NINT behavior.
    """
    import math

    def nint(val: float) -> int:
        return int(math.floor(val + 0.5)) if val >= 0.0 else int(math.ceil(val - 0.5))

    return nint(x * 10000.0) / 10000.0