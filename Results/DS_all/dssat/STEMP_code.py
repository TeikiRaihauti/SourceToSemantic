def _trunc(value: float, decimals: int) -> float:
    """
    Truncate a floating-point number to a fixed number of decimal places (toward zero).
    Types:
    - value: float
    - decimals: int
    Returns:
    - float
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
    DSMID: list,
    ATOT: float,
    TMA: list,
):
    """
    Compute daily soil temperature profile and surface temperature.

    Inputs (types):
    - NL: int
    - ALBEDO: float (not used in calculations but preserved for interface compatibility)
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float (not used in calculations but preserved)
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float (not used in calculations but preserved)
    - WW: float
    - DSMID: list[float] of length NL
    - ATOT: float
    - TMA: list[float] of length 5

    Returns (types):
    - ATOT: float (updated)
    - TMA: list[float] of length 5 (updated)
    - SRFTEMP: float
    - ST: list[float] of length NL
    """
    from math import cos, exp  # imports inside function per constraints

    # Phase angle of annual temperature cycle
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update 5-day moving average array TMA and its sum ATOT
    ATOT = ATOT - TMA[4]
    # Shift right: 1..4 -> 2..5
    for K in range(4, 0, -1):
        TMA[K] = TMA[K - 1]
    # Insert today's average temperature (truncated to 4 decimals)
    TMA[0] = _trunc(TAVG, 4)
    ATOT = ATOT + TMA[0]

    # Soil water content effect
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    FX = exp(B * (((1.0 - WC) / (1.0 + WC)) ** 2))
    DD = FX * DP

    # Mean annual temperature at surface and daily deviation component
    TA = TAV + (TAMP * cos(ALX) / 2.0)
    DT = ATOT / 5.0 - TA

    # Temperature by layer (truncate to 3 decimals)
    ST = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST[L] = TAV + ((TAMP / 2.0 * cos(ALX + ZD) + DT) * exp(ZD))
        ST[L] = _trunc(ST[L], 3)

    # Surface temperature
    SRFTEMP = TAV + (TAMP / 2.0 * cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, ST


def init_stemp(
    NL: int,
    ISWWAT: str,
    BD: list,
    DLAYR: list,
    DS: list,
    DUL: list,
    LL: list,
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: list,
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    DOY: int,
):
    """
    Initialization for soil temperature model.

    Inputs (types):
    - NL: int
    - ISWWAT: str
    - BD: list[float] length NL
    - DLAYR: list[float] length NL
    - DS: list[float] length NL (cumulative depth)
    - DUL: list[float] length NL
    - LL: list[float] length NL
    - NLAYR: int
    - MSALB: float
    - SRAD: float
    - SW: list[float] length NL
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - DOY: int

    Returns (types):
    - CUMDPT: float
    - DSMID: list[float] length NL
    - TDL: float
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NL
    - HDAY: float
    """
    from math import exp, log  # imports inside function per constraints

    # Initialize arrays and state
    SWI = list(SW)
    DSI = list(DS)
    DLI = [0.0] * NLAYR
    DSMID = [0.0] * NLAYR
    ST = [0.0] * NLAYR
    TMA = [0.0] * 5

    # Harvest day based on latitude sign
    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    # Layer thickness from cumulative depths, layer midpoints, and cumulative depth
    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DSI[L]
        else:
            DLI[L] = DSI[L] - DSI[L - 1]
        DSMID[L] = CUMDPT + (DLI[L] * 5.0)
        CUMDPT = CUMDPT + (DLI[L] * 10.0)
        TBD += BD[L] * DLI[L]
        TLL += LL[L] * DLI[L]
        TSW += SWI[L] * DLI[L]
        TDL += DUL[L] * DLI[L]

    # Plant available/extractable soil water
    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Bulk density (profile), damping depth parameters
    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * exp((-5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = log(500.0 / DP)
    ALBEDO = MSALB  # preserved variable, not directly used in calculations below

    # Initialize moving average of soil temperatures
    for I in range(5):
        TMA[I] = _trunc(TAVG, 4)
    ATOT = TMA[0] * 5.0

    # Initialize layer temperatures
    for L in range(NLAYR):
        ST[L] = TAVG

    # Spin-up calls to stabilize temperature profile
    SRFTEMP = 0.0
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD,
            TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def model_stemp(
    NL: int,
    ISWWAT: str,
    BD: list,
    DLAYR: list,
    DS: list,
    DUL: list,
    LL: list,
    NLAYR: int,
    MSALB: float,
    SRAD: float,
    SW: list,
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CUMDPT: float,
    DSMID: list,
    TDL: float,
    TMA: list,
    ATOT: float,
    SRFTEMP: float,
    ST: list,
    DOY: int,
    HDAY: float,
):
    """
    Daily soil temperature update.

    Inputs (types):
    - NL: int
    - ISWWAT: str
    - BD: list[float] length NL
    - DLAYR: list[float] length NL
    - DS: list[float] length NL
    - DUL: list[float] length NL
    - LL: list[float] length NL
    - NLAYR: int
    - MSALB: float
    - SRAD: float
    - SW: list[float] length NL
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - CUMDPT: float
    - DSMID: list[float] length NL
    - TDL: float
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NL
    - DOY: int
    - HDAY: float

    Returns (types) updated state variables:
    - CUMDPT: float (unchanged)
    - DSMID: list[float] length NL (unchanged)
    - TDL: float (updated per original algorithm)
    - TMA: list[float] length 5 (updated)
    - ATOT: float (updated)
    - SRFTEMP: float (updated)
    - ST: list[float] length NL (updated)
    """
    from math import exp, log  # imports inside function per constraints

    # Aggregate profile properties
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL_updated = TDL
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL_updated += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * exp((-5.63 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.356 - (0.144 * ABD)
    B = log(500.0 / DP)
    ALBEDO = MSALB

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL_updated - TLL)

    ATOT, TMA, SRFTEMP_new, ST_new = SOILT(
        NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD,
        TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA
    )

    return CUMDPT, list(DSMID), TDL_updated, TMA, ATOT, SRFTEMP_new, ST_new