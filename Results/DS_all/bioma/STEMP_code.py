def SOILT(NL: int, ALBEDO: float, B: float, CUMDPT: float, DOY: int, DP: float, HDAY: float, NLAYR: int, PESW: float, SRAD: float, TAMP: float, TAV: float, TAVG: float, TMAX: float, WW: float, DSMID: list[float], ATOT: float, TMA: list[float]) -> tuple[float, list[float], float, list[float]]:
    """
    Supporting function computing soil temperatures.
    Inputs:
    - NL: int
    - ALBEDO: float (unused in algorithm but preserved for signature parity)
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float (unused in algorithm but preserved for signature parity)
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float (unused in algorithm but preserved for signature parity)
    - WW: float
    - DSMID: list[float]
    - ATOT: float
    - TMA: list[float] length 5

    Returns:
    - ATOT: float (updated)
    - TMA: list[float] length 5 (updated)
    - SRFTEMP: float
    - ST: list[float] length NL
    """
    import math

    def trunc(x: float, nd: int) -> float:
        f = 10 ** nd
        return int(x * f) / f

    # Angular position for the day relative to harvest day
    ALX = (float(DOY) - HDAY) * 0.01740

    # Update moving average array of last 5 days of average soil temperatures
    ATOT = ATOT - TMA[4]
    for K in range(4, 0, -1):
        TMA[K] = TMA[K - 1]
    TMA[0] = trunc(TAVG, 4)
    ATOT = ATOT + TMA[0]

    # Water content factor
    WC = max(0.010, PESW) / (WW * CUMDPT) * 10.0
    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + (TAMP * math.cos(ALX) / 2.0)
    DT = ATOT / 5.0 - TA

    ST = [0.0] * NL
    for L in range(NLAYR):
        ZD = -(DSMID[L] / DD)
        ST[L] = TAV + ((TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD))
        ST[L] = trunc(ST[L], 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)
    return ATOT, TMA, SRFTEMP, ST


def Init(NL: int,
         ISWWAT: str,
         BD: list[float],
         DS: list[float],
         DUL: list[float],
         LL: list[float],
         NLAYR: int,
         MSALB: float,
         SW: list[float],
         XLAT: float,
         SRAD: float,
         TAVG: float,
         TMAX: float,
         TAV: float,
         TAMP: float,
         DOY: int) -> tuple[float, list[float], float, list[float], float, float, list[float], float]:
    """
    Initialization function.
    Inputs:
    - NL: int (number of soil layers used to size arrays)
    - ISWWAT: str ("Y" for water-limited, otherwise not)
    - BD: list[float] (bulk density per layer)
    - DS: list[float] (cumulative depth per layer)
    - DUL: list[float] (drained upper limit volumetric water content per layer)
    - LL: list[float] (lower limit volumetric water content per layer)
    - NLAYR: int (actual number of soil layers used)
    - MSALB: float (soil albedo)
    - SW: list[float] (volumetric soil water content per layer)
    - XLAT: float (latitude)
    - SRAD: float (solar radiation)
    - TAVG: float (average daily air temperature)
    - TMAX: float (maximum daily air temperature)
    - TAV: float (average annual soil temperature)
    - TAMP: float (amplitude of temperature function)
    - DOY: int (day of year)

    Returns (state variables):
    - CUMDPT: float
    - DSMID: list[float] length NL (only first NLAYR filled)
    - TDL: float
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NL (only first NLAYR filled)
    - HDAY: float
    """
    import math

    def trunc(x: float, nd: int) -> float:
        f = 10 ** nd
        return int(x * f) / f

    # Initialize state variables
    CUMDPT = 0.0
    DSMID = [0.0] * NL
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NL
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    # Layer-wise intermediate values
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0

    # Compute layer thicknesses from cumulative depths and initialize DSMID, CUMDPT, and totals
    for L in range(NLAYR):
        if L == 0:
            DLI = DS[L]
        else:
            DLI = DS[L] - DS[L - 1]
        DSMID[L] = CUMDPT + (DLI * 5.0)
        CUMDPT = CUMDPT + (DLI * 10.0)
        TBD += BD[L] * DLI
        TLL += LL[L] * DLI
        TSW += SW[L] * DLI
        TDL += DUL[L] * DLI

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.630 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.3560 - (0.1440 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # preserved but unused downstream

    # Initialize TMA series and ATOT
    for I in range(5):
        TMA[I] = trunc(TAVG, 4)
    ATOT = TMA[0] * 5.0

    # Initialize soil temperature profile
    for L in range(NLAYR):
        ST[L] = TAVG

    # Spin-up soil temperature with 8 iterations
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA)

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def CalculateModel(NL: int,
                   ISWWAT: str,
                   BD: list[float],
                   DLAYR: list[float],
                   DS: list[float],
                   DUL: list[float],
                   LL: list[float],
                   NLAYR: int,
                   MSALB: float,
                   SW: list[float],
                   SRAD: float,
                   TAVG: float,
                   TMAX: float,
                   TAV: float,
                   TAMP: float,
                   DOY: int,
                   CUMDPT: float,
                   DSMID: list[float],
                   TDL: float,
                   TMA: list[float],
                   ATOT: float,
                   SRFTEMP: float,
                   ST: list[float],
                   HDAY: float) -> tuple[float, list[float], float, list[float], float, float, list[float]]:
    """
    Main biophysical process function (one time step update).
    Inputs:
    - NL: int
    - ISWWAT: str
    - BD: list[float]
    - DLAYR: list[float]
    - DS: list[float]
    - DUL: list[float]
    - LL: list[float]
    - NLAYR: int
    - MSALB: float
    - SW: list[float]
    - SRAD: float
    - TAVG: float
    - TMAX: float
    - TAV: float
    - TAMP: float
    - DOY: int
    - CUMDPT: float
    - DSMID: list[float]
    - TDL: float (carried-over state; accumulated as in original code)
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NL
    - HDAY: float

    Returns (updated state variables):
    - CUMDPT: float (unchanged here)
    - DSMID: list[float] (unchanged here)
    - TDL: float (updated/accumulated)
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NL
    """
    import math

    # Layer totals
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0

    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # accumulated per original code
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + (686.0 * math.exp(-(5.630 * ABD))))
    DP = 1000.0 + (2500.0 * FX)
    WW = 0.3560 - (0.1440 * ABD)
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # preserved but unused downstream

    if ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT, TMA, SRFTEMP, ST = SOILT(NL, ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA)

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST