from typing import List, Tuple


def _nint(x: float) -> int:
    """
    Fortran-like nearest integer with halves away from zero.
    """
    if x >= 0.0:
        return int((x + 0.5) // 1)
    else:
        return int((x - 0.5) // 1)


def YR_DOY(YRDOY: int) -> Tuple[int, int]:
    """
    Convert YYYYDDD date (e.g., 2021100) to (YEAR, DOY).
    Inputs:
      - YRDOY: int
    Returns:
      - YEAR: int
      - DOY: int
    """
    YEAR = int(YRDOY // 1000)
    DOY = int(YRDOY - YEAR * 1000)
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
    DSMID: List[float],
    ATOT: float,
    TMA: List[float],
    SRFTEMP: float,
    ST: List[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    Compute soil temperature profile by layer.
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
      - DSMID: List[float] length >= NLAYR
      - ATOT: float (InOut)
      - TMA: List[float] length 5 (InOut)
      - SRFTEMP: float (Output)
      - ST: List[float] length >= NLAYR (Output)
    Returns:
      - ATOT_new: float
      - TMA_new: List[float] length 5
      - SRFTEMP_new: float
      - ST_new: List[float] length NLAYR
    """
    import math

    # Local copies
    TMA_loc = list(TMA)
    ST_loc = list(ST)

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT_loc = ATOT - TMA_loc[4]

    # Shift TMA history
    for K in range(4, 0, -1):
        TMA_loc[K] = TMA_loc[K - 1]
    TMA_loc[0] = TAVG

    # Fortran NINT to 4 decimals
    TMA_loc[0] = _nint(TMA_loc[0] * 10000.0) / 10000.0
    ATOT_loc = ATOT_loc + TMA_loc[0]

    # Corrected water content equation
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_loc / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_loc[L] = _nint(ST_val * 1000.0) / 1000.0

    SRFTEMP_loc = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_loc, TMA_loc, SRFTEMP_loc, ST_loc


def STEMP(
    control_DYNAMIC: int,
    control_YRDOY: int,
    control_RUN: int,
    control_RNMODE: str,
    iswitch_ISWWAT: str,
    soilprop_BD: List[float],
    soilprop_DLAYR: List[float],
    soilprop_DS: List[float],
    soilprop_DUL: List[float],
    soilprop_LL: List[float],
    soilprop_NLAYR: int,
    soilprop_MSALB: float,
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
    Soil temperature main process function (daily timestep).
    Inputs:
      - control_DYNAMIC: int (1..7; uses 2=SEASINIT, 3=RATE, 5=OUTPUT, 6=SEASEND)
      - control_YRDOY: int
      - control_RUN: int
      - control_RNMODE: str
      - iswitch_ISWWAT: str ('Y'/'N')
      - soilprop_BD: List[float]
      - soilprop_DLAYR: List[float]
      - soilprop_DS: List[float] (cumulative depth, cm)
      - soilprop_DUL: List[float]
      - soilprop_LL: List[float]
      - soilprop_NLAYR: int
      - soilprop_MSALB: float
      - SRAD: float
      - SW: List[float]
      - TAVG: float
      - TMAX: float
      - XLAT: float
      - TAV: float
      - TAMP: float
      - CUMDPT: float (InOut)
      - DSMID: List[float] length >= NLAYR (InOut)
      - TDL: float (InOut)
      - TMA: List[float] length 5 (InOut)
      - ATOT: float (InOut)
      - SRFTEMP: float (InOut)
      - ST: List[float] length >= NLAYR (InOut)
    Returns:
      - CUMDPT_new: float
      - DSMID_new: List[float]
      - TDL_new: float
      - TMA_new: List[float]
      - ATOT_new: float
      - SRFTEMP_new: float
      - ST_new: List[float]
    """
    import math

    # Local copies for outputs
    CUMDPT_loc = CUMDPT
    DSMID_loc = list(DSMID)
    TDL_loc = TDL
    TMA_loc = list(TMA)
    ATOT_loc = ATOT
    SRFTEMP_loc = SRFTEMP
    ST_loc = list(ST)

    # Dynamic phase constants
    SEASINIT = 2
    RATE = 3
    OUTPUT = 5
    SEASEND = 6

    # Extract frequently used data
    NLAYR = soilprop_NLAYR
    BD = soilprop_BD
    DLAYR = soilprop_DLAYR
    DS = soilprop_DS
    DUL = soilprop_DUL
    LL = soilprop_LL
    MSALB = soilprop_MSALB
    ISWWAT = iswitch_ISWWAT

    # Date decomposition
    YEAR, DOY = YR_DOY(control_YRDOY)

    # Seasonal initialization
    if control_DYNAMIC == SEASINIT:
        RUN = control_RUN
        RNMODE = control_RNMODE

        if RUN == 1 or (RNMODE not in "QF"):
            # Initialize depths and water
            SWI = list(SW)
            DSI = list(DS)

            if XLAT < 0.0:
                HDAY = 20.0
            else:
                HDAY = 200.0

            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_loc = 0.0
            CUMDPT_loc = 0.0

            # Compute layer thicknesses from cumulative depths (cm)
            DLI = [0.0] * NLAYR
            for L in range(NLAYR):
                if L == 0:
                    DLI[L] = DSI[L]
                else:
                    DLI[L] = DSI[L] - DSI[L - 1]
                DSMID_loc[L] = CUMDPT_loc + DLI[L] * 5.0  # mm midpoint
                CUMDPT_loc = CUMDPT_loc + DLI[L] * 10.0   # mm profile depth
                TBD = TBD + BD[L] * DLI[L]
                TLL = TLL + LL[L] * DLI[L]
                TSW = TSW + SWI[L] * DLI[L]
                TDL_loc = TDL_loc + DUL[L] * DLI[L]

            if ISWWAT == "Y":
                PESW = max(0.0, TSW - TLL)
            else:
                PESW = max(0.0, TDL_loc - TLL)

            ABD = TBD / DSI[NLAYR - 1]
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP)
            ALBEDO = MSALB

            # Initialize TMA history and ST profile
            for i in range(5):
                TMA_loc[i] = _nint(TAVG * 10000.0) / 10000.0
            ATOT_loc = TMA_loc[0] * 5.0

            for L in range(NLAYR):
                ST_loc[L] = TAVG

            # Spin-up calls
            for _ in range(8):
                ATOT_loc, TMA_loc, SRFTEMP_loc, ST_loc = SOILT(
                    ALBEDO,
                    B,
                    CUMDPT_loc,
                    DOY,
                    DP,
                    HDAY,
                    NLAYR,
                    PESW,
                    SRAD,
                    TAMP,
                    TAV,
                    TAVG,
                    TMAX,
                    WW,
                    DSMID_loc,
                    ATOT_loc,
                    TMA_loc,
                    SRFTEMP_loc,
                    ST_loc,
                )

    # Daily rate calculations
    elif control_DYNAMIC == RATE:
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        for L in range(NLAYR):
            TBD = TBD + BD[L] * DLAYR[L]
            TDL_loc = TDL_loc + DUL[L] * DLAYR[L]
            TLL = TLL + LL[L] * DLAYR[L]
            TSW = TSW + SW[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = MSALB

        if ISWWAT == "Y":
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_loc - TLL)

        # Recompute HDAY from latitude each day (stateless)
        HDAY = 20.0 if XLAT < 0.0 else 200.0

        ATOT_loc, TMA_loc, SRFTEMP_loc, ST_loc = SOILT(
            ALBEDO,
            B,
            CUMDPT_loc,
            DOY,
            DP,
            HDAY,
            NLAYR,
            PESW,
            SRAD,
            TAMP,
            TAV,
            TAVG,
            TMAX,
            WW,
            DSMID_loc,
            ATOT_loc,
            TMA_loc,
            SRFTEMP_loc,
            ST_loc,
        )

    # OUTPUT and SEASEND branches do only I/O in original; omitted here.

    return CUMDPT_loc, DSMID_loc, TDL_loc, TMA_loc, ATOT_loc, SRFTEMP_loc, ST_loc