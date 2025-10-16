def YR_DOY(YRDOY):
    """
    Converts YRDOY (YYYYDDD) to (YEAR, DOY).
    """
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def _is_sequence(x):
    try:
        len(x)  # noqa
        return True
    except Exception:
        return False


def _to_layer_array(val, n):
    """
    Ensure val is a list of length n. If scalar, replicate n times.
    If sequence, copy and truncate/pad as needed (padding with last value).
    """
    if _is_sequence(val):
        arr = [float(v) for v in val[:n]]
        if len(arr) < n:
            if len(arr) == 0:
                arr = [0.0] * n
            else:
                last = arr[-1]
                arr.extend([last] * (n - len(arr)))
        return arr
    else:
        return [float(val)] * n


def _fortran_nint(x):
    """
    Fortran NINT behavior: nearest integer; ties away from zero.
    """
    if x >= 0.0:
        return int(float(int(x)) + (1.0 if (x - float(int(x))) >= 0.5 else 0.0))
    else:
        ax = -x
        return -int(float(int(ax)) + (1.0 if (ax - float(int(ax))) >= 0.5 else 0.0))


def _round_via_nint(x, decimals):
    """
    Round x to 'decimals' places using Fortran NINT logic on scaled value.
    """
    scale = 10.0 ** decimals
    return _fortran_nint(x * scale) / scale


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA, SRFTEMP, ST):
    """
    Soil temperature by layer.
    Inputs:
      ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
      PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID
    InOut:
      ATOT, TMA
    Outputs:
      SRFTEMP, ST
    Returns: (ATOT, TMA, SRFTEMP, ST)
    """
    # Ensure copies to keep function pure
    TMA_local = [float(x) for x in _to_layer_array(TMA, 5)]
    DSMID_local = [float(x) for x in _to_layer_array(DSMID, NLAYR)]
    ST_local = [float(x) for x in _to_layer_array(ST, NLAYR)]
    ATOT_local = float(ATOT)

    # Angular day offset
    ALX = (float(DOY) - float(HDAY)) * 0.0174

    # Update running sum of previous 5 days average air temps (TMA)
    ATOT_local = ATOT_local - TMA_local[4]
    for K in range(4, 0, -1):
        TMA_local[K] = TMA_local[K - 1]
    TMA_local[0] = float(TAVG)
    # Prevent release/debug differences: keep only 4 decimals
    TMA_local[0] = _round_via_nint(TMA_local[0], 4)
    ATOT_local = ATOT_local + TMA_local[0]

    # Corrected water content function (EPIC corrected)
    WC = max(0.01, float(PESW)) / (float(WW) * float(CUMDPT)) * 10.0

    FX = __import__("math").exp(float(B) * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * float(DP)  # Damping depth in mm

    TA = float(TAV) + float(TAMP) * __import__("math").cos(ALX) / 2.0
    DT = ATOT_local / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID_local[L] / DD
        ST_val = float(TAV) + (float(TAMP) / 2.0 * __import__("math").cos(ALX + ZD) + DT) * __import__("math").exp(ZD)
        # Prevent debug vs release differences: keep only 3 decimals
        ST_local[L] = _round_via_nint(ST_val, 3)

    # Surface temperature
    SRFTEMP_local = float(TAV) + (float(TAMP) / 2.0 * __import__("math").cos(ALX) + DT)

    return ATOT_local, TMA_local, SRFTEMP_local, ST_local


def STEMP_Initialize(CONTROL_RUN, CONTROL_RNMODE, CONTROL_YRDOY,
                     ISWITCH_ISWWAT,
                     SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
                     SOILPROP_NLAYR, SOILPROP_MSALB,
                     SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
                     CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY):
    """
    Seasonal initialization for STEMP.
    Returns updated state variables:
      CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY
    """

    # Determine whether to (re)initialize
    # Fortran: IF (RUN .EQ. 1 .OR. INDEX('QF',RNMODE) .LE. 0) THEN
    rnmode_str = '' if CONTROL_RNMODE is None else str(CONTROL_RNMODE)
    do_init = (int(CONTROL_RUN) == 1) or (rnmode_str.find('Q') == -1 and rnmode_str.find('F') == -1)

    # If not initializing, return unchanged states
    if not do_init:
        return (CUMDPT, list(DSMID), TDL, list(_to_layer_array(TMA, 5)),
                ATOT, SRFTEMP, list(ST), HDAY)

    NLAYR = int(SOILPROP_NLAYR)
    BD = _to_layer_array(SOILPROP_BD, NLAYR)
    DS = _to_layer_array(SOILPROP_DS, NLAYR)
    DUL = _to_layer_array(SOILPROP_DUL, NLAYR)
    LL = _to_layer_array(SOILPROP_LL, NLAYR)
    SWI = _to_layer_array(SW, NLAYR)  # Initial soil water
    DSMID_local = [0.0] * NLAYR

    # HDAY: day of hottest day depends on hemisphere
    if float(XLAT) < 0.0:
        HDAY_local = 20.0
    else:
        HDAY_local = 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL_local = 0.0
    CUMDPT_local = 0.0

    # Compute layer thickness from cumulative depths and midpoints
    for L in range(NLAYR):
        if L == 0:
            DLI = DS[L]
        else:
            DLI = DS[L] - DS[L - 1]  # cm
        DSMID_local[L] = CUMDPT_local + DLI * 5.0       # mm to midpoint
        CUMDPT_local = CUMDPT_local + DLI * 10.0        # mm cumulative depth

        TBD += BD[L] * DLI
        TLL += LL[L] * DLI
        TSW += SWI[L] * DLI
        TDL_local += DUL[L] * DLI

    # Potential extractable soil water (cm)
    if str(ISWITCH_ISWWAT).upper() == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL_local - TLL)

    ABD = TBD / DS[NLAYR - 1]  # Average bulk density
    FX = ABD / (ABD + 686.0 * __import__("math").exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX               # mm
    WW = 0.356 - 0.144 * ABD                # volumetric fraction
    B = __import__("math").log(500.0 / DP)
    ALBEDO = float(SOILPROP_MSALB)

    # Initialize TMA (last 5 days average temperature) to TAVG with 4 decimals
    TMA_local = [_round_via_nint(float(TAVG), 4) for _ in range(5)]
    ATOT_local = TMA_local[0] * 5.0

    # Initialize soil temperatures to TAVG
    ST_local = [float(TAVG)] * NLAYR
    SRFTEMP_local = float(TAVG)

    # Compute DOY from YRDOY
    _, DOY = YR_DOY(int(CONTROL_YRDOY))

    # Spin-up: call SOILT 8 times
    for _ in range(8):
        ATOT_local, TMA_local, SRFTEMP_local, ST_local = SOILT(
            ALBEDO, B, CUMDPT_local, DOY, DP, HDAY_local, NLAYR,
            PESW, float(SRAD), float(TAMP), float(TAV), float(TAVG), float(TMAX), WW, DSMID_local,
            ATOT_local, TMA_local, SRFTEMP_local, ST_local
        )

    return (CUMDPT_local, DSMID_local, TDL_local, TMA_local,
            ATOT_local, SRFTEMP_local, ST_local, HDAY_local)


def STEMP(ISWITCH_ISWWAT,
          SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
          SOILPROP_NLAYR, SOILPROP_MSALB,
          SRAD, SW, TAVG, TMAX, TAV, TAMP,
          CONTROL_YRDOY,
          CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY):
    """
    Main daily soil temperature process.
    Inputs:
      ISWITCH_ISWWAT
      SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL
      SOILPROP_NLAYR, SOILPROP_MSALB
      SRAD, SW, TAVG, TMAX, TAV, TAMP
      CONTROL_YRDOY
      CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY (state variables)
    Returns updated state variables:
      CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY
    """

    NLAYR = int(SOILPROP_NLAYR)
    BD = _to_layer_array(SOILPROP_BD, NLAYR)
    DLAYR = _to_layer_array(SOILPROP_DLAYR, NLAYR)
    DS = _to_layer_array(SOILPROP_DS, NLAYR)
    DUL = _to_layer_array(SOILPROP_DUL, NLAYR)
    LL = _to_layer_array(SOILPROP_LL, NLAYR)
    SW_arr = _to_layer_array(SW, NLAYR)
    DSMID_local = [float(x) for x in _to_layer_array(DSMID, NLAYR)]
    TMA_local = [float(x) for x in _to_layer_array(TMA, 5)]
    ST_local = [float(x) for x in _to_layer_array(ST, NLAYR)]
    ATOT_local = float(ATOT)
    TDL_local = float(TDL)
    SRFTEMP_local = float(SRFTEMP)

    # Compute yearly day
    _, DOY = YR_DOY(int(CONTROL_YRDOY))

    # Daily sums
    TBD = 0.0
    TLL = 0.0
    TSW_sum = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        # Note: Original code does not reset TDL before accumulation
        TDL_local += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW_sum += SW_arr[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]            # Average bulk density
    FX = ABD / (ABD + 686.0 * __import__("math").exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX            # mm
    WW = 0.356 - 0.144 * ABD             # volumetric fraction
    B = __import__("math").log(500.0 / DP)
    ALBEDO = float(SOILPROP_MSALB)

    # Potential extractable soil water (cm)
    if str(ISWITCH_ISWWAT).upper() == 'Y':
        PESW = max(0.0, TSW_sum - TLL)
    else:
        PESW = max(0.0, TDL_local - TLL)

    # Update soil temperatures by layer
    ATOT_local, TMA_local, SRFTEMP_local, ST_local = SOILT(
        ALBEDO, B, float(CUMDPT), DOY, DP, float(HDAY), NLAYR,
        PESW, float(SRAD), float(TAMP), float(TAV), float(TAVG), float(TMAX), WW, DSMID_local,
        ATOT_local, TMA_local, SRFTEMP_local, ST_local
    )

    return (float(CUMDPT), DSMID_local, TDL_local, TMA_local,
            ATOT_local, SRFTEMP_local, ST_local, float(HDAY))