def YR_DOY(YRDOY: int) -> tuple:
    """
    Convert YYYYDOY integer to (year, day-of-year).
    Inputs:
      - YRDOY: int
    Returns:
      - (YR, DOY): tuple of ints
    """
    YR = YRDOY // 1000
    DOY = YRDOY - YR * 1000
    return YR, DOY


def _nint(x: float) -> int:
    """
    Fortran-like NINT: nearest integer, halves away from zero.
    """
    if x >= 0.0:
        return int(x + 0.5)
    else:
        return -int(abs(x) + 0.5)


def _round_nint_decimals(x: float, decimals: int) -> float:
    """
    Fortran-like rounding using NINT at specified decimals.
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
    DSMID: list,
    ATOT: float,
    TMA: list,
    SRFTEMP: float,
    ST: list,
) -> tuple:
    """
    Soil temperature by layer (daily).
    Inputs:
      - ALBEDO: float (not used in algorithm, kept for interface fidelity)
      - B: float
      - CUMDPT: float
      - DOY: int
      - DP: float
      - HDAY: float
      - NLAYR: int
      - PESW: float
      - SRAD: float (not used in algorithm)
      - TAMP: float
      - TAV: float
      - TAVG: float
      - TMAX: float (not used in algorithm)
      - WW: float
      - DSMID: list[float]
      - ATOT: float
      - TMA: list[float] length 5
      - SRFTEMP: float
      - ST: list[float] length >= NLAYR
    Returns (updated):
      - ATOT: float
      - TMA: list[float] length 5
      - SRFTEMP: float
      - ST: list[float] length NLAYR
    """
    import math

    # Phase offset by day of year
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average accumulator and queue
    ATOT = ATOT - TMA[4]
    # Shift TMA (5->2)
    for k in range(4, 0, -1):
        TMA[k] = TMA[k - 1]
    # Insert today's mean temp with 4-decimal NINT rounding
    TMA[0] = _round_nint_decimals(TAVG, 4)
    ATOT = ATOT + TMA[0]

    # Corrected water content effect (EPIC corrected)
    # WC is dimensionless ratio
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # damping depth (mm)

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    # Layer temperatures
    st_out = list(ST[:NLAYR])
    for l in range(NLAYR):
        ZD = -DSMID[l] / DD
        st_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        st_out[l] = _round_nint_decimals(st_val, 3)

    # Surface temperature (litter/soil surface)
    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, st_out


def STEMP_Initialize(
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: list,
    SOILPROP_DS: list,
    SOILPROP_DUL: list,
    SOILPROP_LL: list,
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
    SRAD: float,
    SW: list,
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CONTROL_YRDOY: int,
) -> tuple:
    """
    Seasonal initialization for soil temperature state.
    Inputs:
      - ISWITCH_ISWWAT: str ('Y' if soil water simulated, else not)
      - SOILPROP_BD: list[float], bulk density by layer
      - SOILPROP_DS: list[float], cumulative depth to layer bottom (cm)
      - SOILPROP_DUL: list[float], drained upper limit (vol. frac.) by layer
      - SOILPROP_LL: list[float], lower limit (vol. frac.) by layer
      - SOILPROP_NLAYR: int, number of soil layers
      - SOILPROP_MSALB: float, albedo including mulch/soil water effects
      - SRAD: float, solar radiation (MJ/m2/d) [not used directly]
      - SW: list[float], initial soil water content (vol. frac.) by layer
      - TAVG: float, daily average air temperature (C)
      - TMAX: float, daily max air temperature (C) [not used directly]
      - XLAT: float, latitude (deg); sign used to set hottest day
      - TAV: float, annual average air temperature (C)
      - TAMP: float, annual temperature amplitude (C)
      - CONTROL_YRDOY: int, YYYYDOY date
    Returns (state):
      - CUMDPT: float, cumulative profile depth (mm)
      - DSMID: list[float], depth to mid-point of each layer (mm)
      - TDL: float, profile sum of DUL*layer_thickness (cm)
      - TMA: list[float] length 5, moving average queue of TAVG
      - ATOT: float, sum of TMA queue
      - SRFTEMP: float, surface temperature (C)
      - ST: list[float], soil temperature by layer (C)
      - HDAY: float, day-of-year of hottest day (phase anchor)
    """
    # Day-of-year for initialization/spin-up
    _, DOY = YR_DOY(CONTROL_YRDOY)

    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    SWI = list(SW[:NLAYR])  # initial soil water by layer
    MSALB = float(SOILPROP_MSALB)

    # Hemisphere-dependent hottest day
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    # Initialize accumulators
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    # Compute layer thickness (cm) from cumulative DS and DSMID (mm)
    DLI = [0.0] * NLAYR
    DSMID = [0.0] * NLAYR
    for l in range(NLAYR):
        if l == 0:
            DLI[l] = DS[l]
        else:
            DLI[l] = DS[l] - DS[l - 1]
        DSMID[l] = CUMDPT + DLI[l] * 5.0  # mm to layer midpoint
        CUMDPT = CUMDPT + DLI[l] * 10.0   # mm cumulative depth

        TBD += BD[l] * DLI[l]
        TLL += LL[l] * DLI[l]
        TSW += SWI[l] * DLI[l]
        TDL += DUL[l] * DLI[l]

    # Plant-available water in profile (cm)
    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        # If water not simulated, use DUL as water content
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]  # bulk density averaged over profile
    FX = ABD / (ABD + 686.0 * __import__("math").exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # damping depth (mm)
    WW = 0.356 - 0.144 * ABD   # water capacity (vol. fraction)
    B = __import__("math").log(500.0 / DP)
    ALBEDO = MSALB

    # Initialize moving average of TAVG (rounded to 4 decimals)
    TMA = [_round_nint_decimals(TAVG, 4) for _ in range(5)]
    ATOT = TMA[0] * 5.0

    # Initialize soil temperatures by layer
    ST = [float(TAVG) for _ in range(NLAYR)]
    SRFTEMP = float(TAVG)

    # Spin-up SOILT 8 times to stabilize temperatures
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            ALBEDO=ALBEDO,
            B=B,
            CUMDPT=CUMDPT,
            DOY=int(DOY),
            DP=DP,
            HDAY=HDAY,
            NLAYR=NLAYR,
            PESW=PESW,
            SRAD=SRAD,
            TAMP=TAMP,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            WW=WW,
            DSMID=DSMID,
            ATOT=ATOT,
            TMA=TMA,
            SRFTEMP=SRFTEMP,
            ST=ST,
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def STEMP(
    CONTROL_YRDOY: int,
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: list,
    SOILPROP_DLAYR: list,
    SOILPROP_DS: list,
    SOILPROP_DUL: list,
    SOILPROP_LL: list,
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
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
    HDAY: float,
) -> tuple:
    """
    Daily soil temperature process (RATE).
    Inputs:
      - CONTROL_YRDOY: int, YYYYDOY
      - ISWITCH_ISWWAT: str, 'Y' if simulating water, else other
      - SOILPROP_BD: list[float], bulk density by layer
      - SOILPROP_DLAYR: list[float], layer thickness (same length as layers)
      - SOILPROP_DS: list[float], cumulative depth to layer bottom (cm)
      - SOILPROP_DUL: list[float], DUL by layer (vol. frac.)
      - SOILPROP_LL: list[float], LL by layer (vol. frac.)
      - SOILPROP_NLAYR: int, number of layers
      - SOILPROP_MSALB: float, albedo incl. mulch/soil water
      - SRAD: float, solar radiation (not used by core equations)
      - SW: list[float], soil water content by layer (vol. frac.)
      - TAVG: float, average air temperature (C)
      - TMAX: float, max air temperature (C) [not used directly]
      - XLAT: float, latitude (deg) [not used directly here]
      - TAV: float, annual average temperature (C)
      - TAMP: float, annual temperature amplitude (C)
      - CUMDPT: float, cumulative profile depth (mm)
      - DSMID: list[float], depth to mid-point of each layer (mm)
      - TDL: float, previous day's TDL (cm); recalculated daily here
      - TMA: list[float] length 5, moving average queue of TAVG
      - ATOT: float, sum of TMA queue
      - SRFTEMP: float, previous surface temperature
      - ST: list[float], previous soil temperature by layer
      - HDAY: float, hottest-day phase anchor
    Returns (updated):
      - TDL: float, profile sum of DUL*layer_thickness (cm) for the day
      - TMA: list[float] length 5
      - ATOT: float
      - SRFTEMP: float
      - ST: list[float]
    """
    import math

    _, DOY = YR_DOY(CONTROL_YRDOY)

    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SW_in = list(SW[:NLAYR])

    # Aggregates over profile
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL_local = 0.0  # recompute daily

    for l in range(NLAYR):
        TBD += BD[l] * DLAYR[l]
        TDL_local += DUL[l] * DLAYR[l]
        TLL += LL[l] * DLAYR[l]
        TSW += SW_in[l] * DLAYR[l]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL_local - TLL)  # cm

    ATOT, TMA, SRFTEMP, ST_out = SOILT(
        ALBEDO=ALBEDO,
        B=B,
        CUMDPT=CUMDPT,
        DOY=int(DOY),
        DP=DP,
        HDAY=HDAY,
        NLAYR=NLAYR,
        PESW=PESW,
        SRAD=SRAD,
        TAMP=TAMP,
        TAV=TAV,
        TAVG=TAVG,
        TMAX=TMAX,
        WW=WW,
        DSMID=DSMID,
        ATOT=ATOT,
        TMA=TMA,
        SRFTEMP=SRFTEMP,
        ST=ST,
    )

    return TDL_local, TMA, ATOT, SRFTEMP, ST_out