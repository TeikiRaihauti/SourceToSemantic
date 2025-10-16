def YR_DOY(YRDOY):
    """
    Converts YRDOY (YYYYDDD) to (YEAR, DOY).
    """
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def _nint(x):
    """
    Fortran NINT: nearest integer, halves away from zero.
    """
    import math
    if x >= 0:
        return math.floor(x + 0.5)
    else:
        return -math.floor(-x + 0.5)


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA):
    """
    Soil temperature by layer. Supports STEMP.

    Inputs:
      ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID (list)

    InOut:
      ATOT (float), TMA (list length 5)

    Outputs:
      ATOT, TMA, SRFTEMP, ST (list length NLAYR)
    """
    import math

    # Phase angle for annual temperature wave
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average array TMA (length 5) and ATOT (sum of last 5 TAVG)
    ATOT = ATOT - TMA[4]
    # Shift TMA down
    for K in range(4, 0, -1):
        TMA[K] = TMA[K - 1]
    # Set newest average daily temperature
    TMA[0] = TAVG
    # Keep only 4 decimals (debug vs release parity)
    TMA[0] = _nint(TMA[0] * 10000.0) / 10000.0
    ATOT = ATOT + TMA[0]

    # Corrected water content function (EPIC method)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)

    DD = FX * DP  # Damping depth in mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    ST = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST[L] = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        # Keep only 3 decimals (debug vs release parity)
        ST[L] = _nint(ST[L] * 1000.0) / 1000.0

    # Surface temperature (temporary approach)
    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, ST


def STEMP_Initialize(CONTROL_YRDOY,
                     ISWITCH_ISWWAT,
                     SOILPROP_BD,
                     SOILPROP_DS,
                     SOILPROP_DUL,
                     SOILPROP_LL,
                     SOILPROP_NLAYR,
                     SOILPROP_MSALB,
                     SRAD,
                     SW,
                     TAVG,
                     TMAX,
                     XLAT,
                     TAV,
                     TAMP):
    """
    Seasonal initialization for STEMP component.

    Inputs:
      CONTROL_YRDOY           int YYYYDDD
      ISWITCH_ISWWAT          'Y' or 'N'
      SOILPROP_BD             list of bulk density (g/cm3) per layer
      SOILPROP_DS             list cumulative depth (cm) per layer
      SOILPROP_DUL            list drained upper limit (cm3/cm3) per layer
      SOILPROP_LL             list lower limit (cm3/cm3) per layer
      SOILPROP_NLAYR          int number of soil layers
      SOILPROP_MSALB          float mulch/soil albedo
      SRAD                    float daily solar radiation (MJ/m2-d) [passed to spin-up]
      SW                      list initial soil water (cm3/cm3) per layer
      TAVG                    float average daily air temperature (C)
      TMAX                    float max daily air temperature (C) [not used directly]
      XLAT                    float latitude (deg, sign determines hemisphere)
      TAV                     float average annual soil temperature (C)
      TAMP                    float amplitude of annual soil temperature (C)

    Returns (state variables):
      CUMDPT                  float cumulative depth of soil profile (mm)
      DSMID                   list depth to midpoint of each soil layer (mm)
      TDL                     float total water content at DUL (cm) [state, cumulative as in original]
      TMA                     list length 5 of last 5 days average soil temperatures (C)
      ATOT                    float sum of TMA array (C)
      SRFTEMP                 float surface soil temperature (C)
      ST                      list soil temperature by layer (C)
      HDAY                    float day-of-year of climatological hottest day (20 SH, 200 NH)
    """
    import math

    # Unpack
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DSI = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SWI = list(SW[:NLAYR])

    # Determine hottest day offset by hemisphere
    if XLAT < 0.0:
        HDAY = 20.0   # Southern hemisphere
    else:
        HDAY = 200.0  # Northern hemisphere

    # Compute layer thickness (cm) from cumulative depths and midpoints, CUMDPT (mm)
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    DLI = [0.0] * NLAYR  # layer thickness (cm)

    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DSI[L]
        else:
            DLI[L] = DSI[L] - DSI[L - 1]
        DSMID[L] = CUMDPT + DLI[L] * 5.0   # midpoint depth (mm)
        CUMDPT = CUMDPT + DLI[L] * 10.0    # cumulative profile depth (mm)

        TBD += BD[L] * DLI[L]
        TLL += LL[L] * DLI[L]
        TSW += SWI[L] * DLI[L]
        TDL += DUL[L] * DLI[L]

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm (use DUL if water not simulated)

    ABD = TBD / DSI[NLAYR - 1]  # average bulk density (g/cm3)
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # not used in SOILT but preserved

    # Initialize TMA (5-day mean) and ATOT, layer temps ST
    TMA = [0.0] * 5
    for I in range(5):
        TMA[I] = _nint(TAVG * 10000.0) / 10000.0
    ATOT = TMA[0] * 5.0

    ST = [TAVG for _ in range(NLAYR)]

    # Compute DOY from YRDOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Spin-up SOILT 8 times to stabilize temperatures
    SRFTEMP = 0.0
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
            ATOT, TMA
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def STEMP(CONTROL_YRDOY,
          ISWITCH_ISWWAT,
          SOILPROP_BD,
          SOILPROP_DLAYR,
          SOILPROP_DS,
          SOILPROP_DUL,
          SOILPROP_LL,
          SOILPROP_NLAYR,
          SOILPROP_MSALB,
          SRAD,
          SW,
          TAVG,
          TMAX,
          TAV,
          TAMP,
          CUMDPT,
          DSMID,
          TDL,
          TMA,
          ATOT,
          SRFTEMP,
          ST,
          HDAY):
    """
    Main daily soil temperature process (RATE step).

    Inputs:
      CONTROL_YRDOY           int YYYYDDD
      ISWITCH_ISWWAT          'Y' or 'N'
      SOILPROP_BD             list bulk density (g/cm3) per layer
      SOILPROP_DLAYR          list layer thickness (cm) per layer
      SOILPROP_DS             list cumulative depth (cm) per layer
      SOILPROP_DUL            list drained upper limit (cm3/cm3) per layer
      SOILPROP_LL             list lower limit (cm3/cm3) per layer
      SOILPROP_NLAYR          int number of layers
      SOILPROP_MSALB          float mulch/soil albedo
      SRAD                    float solar radiation (MJ/m2-d)
      SW                      list volumetric soil water (cm3/cm3) per layer
      TAVG                    float average daily temperature (C)
      TMAX                    float maximum daily temperature (C)
      TAV                     float average annual soil temperature (C)
      TAMP                    float amplitude of annual soil temperature (C)

    State In:
      CUMDPT                  float cumulative profile depth (mm)
      DSMID                   list midpoint depths (mm)
      TDL                     float total water at DUL (cm) [carried as in original code]
      TMA                     list length 5 of last 5 days' TAVG (C)
      ATOT                    float sum of TMA (C)
      SRFTEMP                 float previous surface temperature (C) [not used as input]
      ST                      list previous soil temperatures (C)
      HDAY                    float hottest day-of-year used for phase

    Returns (updated states):
      CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY
    """
    import math

    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SW_in = list(SW[:NLAYR])

    # Accumulate over layers
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]  # Note: not reset daily, matches original code
        TLL += LL[L] * DLAYR[L]
        TSW += SW_in[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm (use DUL if water not simulated)

    # Compute DOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    ATOT, TMA, SRFTEMP, ST = SOILT(
        ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
        PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
        ATOT, TMA
    )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY