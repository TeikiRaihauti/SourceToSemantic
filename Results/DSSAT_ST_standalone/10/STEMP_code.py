import math

def YR_DOY(CONTROL_YRDOY):
    """
    Converts YRDOY (YYYYDDD) to (YEAR, DOY).
    """
    YEAR = int(CONTROL_YRDOY // 1000)
    DOY = int(CONTROL_YRDOY - YEAR * 1000)
    return YEAR, DOY


def _round_half_away_from_zero(x, decimals=0):
    """
    Fortran-like rounding to nearest with ties away from zero.
    """
    factor = 10.0 ** decimals
    y = x * factor
    if y >= 0:
        return math.floor(y + 0.5) / factor
    else:
        return -math.floor(-y + 0.5) / factor


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA, SRFTEMP, ST):
    """
    Determines soil temperature by layer.

    Inputs
    - ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR
    - PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID
    InOut
    - ATOT, TMA
    Outputs
    - SRFTEMP, ST
    """
    # Ensure local copies
    TMA = list(TMA)
    ST = list(ST)
    DSMID = list(DSMID)

    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average array of past 5 days
    ATOT = ATOT - TMA[4]
    for k in range(4, 0, -1):
        TMA[k] = TMA[k - 1]
    TMA[0] = TAVG
    TMA[0] = _round_half_away_from_zero(TMA[0], 4)
    ATOT = ATOT + TMA[0]

    # Corrected water content equation (EPIC)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)

    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        st_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[L] = _round_half_away_from_zero(st_val, 3)

    # Surface litter temperature
    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, ST


def STEMP_Initialize(CONTROL_YRDOY, ISWITCH_ISWWAT, 
                     SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL, SOILPROP_NLAYR, SOILPROP_MSALB,
                     SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP):
    """
    Seasonal initialization for soil temperature (STEMP).

    Inputs:
    - CONTROL_YRDOY: integer YYYYDDD
    - ISWITCH_ISWWAT: 'Y' or 'N'
    - SOILPROP_BD: list of bulk density per layer (g/cm3)
    - SOILPROP_DS: list cumulative depth to bottom of layer (cm)
    - SOILPROP_DUL: list of volumetric soil water content at DUL per layer (cm3/cm3)
    - SOILPROP_LL: list of volumetric soil water content at lower limit per layer (cm3/cm3)
    - SOILPROP_NLAYR: number of soil layers
    - SOILPROP_MSALB: mulch/soil albedo
    - SRAD: solar radiation (MJ/m2-d)
    - SW: list of volumetric soil water content initial per layer (cm3/cm3)
    - TAVG: average daily temperature (degC)
    - TMAX: maximum daily temperature (degC)
    - XLAT: latitude (deg)
    - TAV: average annual soil temperature (degC)
    - TAMP: amplitude of soil temperature function (degC)

    Returns:
    - CUMDPT: cumulative depth of soil profile (mm)
    - DSMID: list of depth to midpoint of soil layer (mm)
    - TDL: total water content at DUL over profile (cm)
    - TMA: list of previous 5 days average soil temps (degC)
    - ATOT: sum of TMA (degC)
    - SRFTEMP: soil surface temperature (degC)
    - ST: list of soil temperature per layer (degC)
    """
    NLAYR = SOILPROP_NLAYR
    YEAR, DOY = YR_DOY(CONTROL_YRDOY)

    # Hemisphere hottest day-of-year
    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    # Initial values
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    DSI = list(SOILPROP_DS)
    DLI = [0.0] * NLAYR
    SWI = list(SW)

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0

    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DSI[L]
        else:
            DLI[L] = DSI[L] - DSI[L - 1]

        DSMID[L] = CUMDPT + DLI[L] * 5.0  # mm depth to midpt of layer
        CUMDPT = CUMDPT + DLI[L] * 10.0   # mm profile depth

        TBD += SOILPROP_BD[L] * DLI[L]
        TLL += SOILPROP_LL[L] * DLI[L]
        TSW += SWI[L] * DLI[L]
        TDL += SOILPROP_DUL[L] * DLI[L]

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = SOILPROP_MSALB

    # Initialize TMA and ATOT
    TMA = [0.0] * 5
    for i in range(5):
        TMA[i] = _round_half_away_from_zero(TAVG, 4)
    ATOT = TMA[0] * 5.0

    # Initialize soil temperature by layer
    ST = [TAVG for _ in range(NLAYR)]
    SRFTEMP = TAVG

    # Spin-up calls to SOILT to stabilize thermal state
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
            ATOT, TMA, SRFTEMP, ST
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST


def STEMP(CONTROL_YRDOY, ISWITCH_ISWWAT,
          SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL, SOILPROP_NLAYR, SOILPROP_MSALB,
          SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
          CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST):
    """
    Daily rate calculations for soil temperature (STEMP).

    Inputs:
    - CONTROL_YRDOY: integer YYYYDDD
    - ISWITCH_ISWWAT: 'Y' or 'N'
    - SOILPROP_BD: list of bulk density per layer (g/cm3)
    - SOILPROP_DLAYR: list of thickness of soil layer (cm)
    - SOILPROP_DS: list cumulative depth to bottom of layer (cm)
    - SOILPROP_DUL: list of volumetric soil water content at DUL per layer (cm3/cm3)
    - SOILPROP_LL: list of volumetric soil water content at lower limit per layer (cm3/cm3)
    - SOILPROP_NLAYR: number of soil layers
    - SOILPROP_MSALB: mulch/soil albedo
    - SRAD: solar radiation (MJ/m2-d)
    - SW: list of volumetric soil water content per layer (cm3/cm3)
    - TAVG: average daily temperature (degC)
    - TMAX: maximum daily temperature (degC)
    - XLAT: latitude (deg)
    - TAV: average annual soil temperature (degC)
    - TAMP: amplitude of soil temperature function (degC)
    - CUMDPT: cumulative depth of soil profile (mm)
    - DSMID: list of depth to midpoint of soil layer (mm)
    - TDL: total water content of soil at drained upper limit (cm)
    - TMA: list of previous 5 days average soil temps (degC)
    - ATOT: sum of TMA (degC)
    - SRFTEMP: soil surface temperature (degC)
    - ST: list of soil temperature per layer (degC)

    Returns:
    - TDL: updated total water content at DUL over profile (cm)
    - TMA: updated TMA
    - ATOT: updated ATOT
    - SRFTEMP: updated surface temperature
    - ST: updated soil layer temperatures
    """
    NLAYR = SOILPROP_NLAYR
    _, DOY = YR_DOY(CONTROL_YRDOY)

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0

    # Note: TDL is not reset here to preserve original behavior
    for L in range(NLAYR):
        TBD += SOILPROP_BD[L] * SOILPROP_DLAYR[L]
        TDL += SOILPROP_DUL[L] * SOILPROP_DLAYR[L]
        TLL += SOILPROP_LL[L] * SOILPROP_DLAYR[L]
        TSW += SW[L] * SOILPROP_DLAYR[L]

    ABD = TBD / SOILPROP_DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = SOILPROP_MSALB

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    # Hemisphere hottest day-of-year
    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    ATOT, TMA, SRFTEMP, ST = SOILT(
        ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
        PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
        ATOT, TMA, SRFTEMP, ST
    )

    return TDL, TMA, ATOT, SRFTEMP, ST


def test_STEMP_initialization_and_rate():
    """
    Derived test from ASKEE example values.
    This test does not assert exact numeric results but checks shapes and basic consistency.
    """
    # Inputs based on ASKEE.for example
    CONTROL_YRDOY = 2021100
    ISWITCH_ISWWAT = 'Y'
    NLAYR = 4
    SOILPROP_BD = [1.6] * NLAYR
    SOILPROP_DLAYR = [10.0] * NLAYR
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * NLAYR
    SOILPROP_LL = [0.2] * NLAYR
    SOILPROP_NLAYR = NLAYR
    SOILPROP_MSALB = 0.13

    SRAD = 20.0
    SW = [0.2] * NLAYR
    TAVG = 25.0
    TMAX = 30.0
    XLAT = 28.0
    TAV = 20.0
    TAMP = 10.0

    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP_Initialize(
        CONTROL_YRDOY, ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL, SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP
    )

    assert len(DSMID) == NLAYR
    assert DSMID == [50.0, 150.0, 250.0, 350.0]
    assert abs(CUMDPT - 400.0) < 1e-6
    assert len(TMA) == 5
    assert len(ST) == NLAYR
    assert isinstance(SRFTEMP, float)
    assert abs(ATOT - TMA[0] * 5.0) < 1e-6

    TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        CONTROL_YRDOY, ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL, SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
        CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    )

    assert len(TMA) == 5
    assert len(ST) == NLAYR
    assert isinstance(SRFTEMP, float)