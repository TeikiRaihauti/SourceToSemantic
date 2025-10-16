def _nint(x):
    if x >= 0:
        return float(int(x + 0.5))
    else:
        return float(int(x - 0.5))


def YR_DOY(YRDOY):
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA):
    # Copy inputs to avoid mutation
    TMA_local = [float(v) for v in TMA]
    ATOT_local = float(ATOT)
    ST = [0.0] * NLAYR

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT_local = ATOT_local - TMA_local[4]

    # shift TMA
    for k in range(4, 0, -1):
        TMA_local[k] = TMA_local[k - 1]

    TMA_local[0] = float(TAVG)

    # Keep only 4 decimals (debug vs release equivalence)
    TMA_local[0] = _nint(TMA_local[0] * 10000.0) / 10000.0
    ATOT_local = ATOT_local + TMA_local[0]

    # Corrected EPIC water content function
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = pow((1.0 - WC) / (1.0 + WC), 2.0)
    FX = pow(2.718281828459045, B * FX)

    DD = FX * DP  # mm

    TA = TAV + TAMP * __import__('math').cos(ALX) / 2.0
    DT = ATOT_local / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        STL = TAV + (TAMP / 2.0 * __import__('math').cos(ALX + ZD) + DT) * __import__('math').exp(ZD)
        ST[L] = _nint(STL * 1000.0) / 1000.0

    SRFTEMP = TAV + (TAMP / 2.0 * __import__('math').cos(ALX) + DT)

    return ATOT_local, TMA_local, SRFTEMP, ST


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
    NLAYR = int(SOILPROP_NLAYR)
    BD = [float(v) for v in SOILPROP_BD[:NLAYR]]
    DS = [float(v) for v in SOILPROP_DS[:NLAYR]]
    DUL = [float(v) for v in SOILPROP_DUL[:NLAYR]]
    LL = [float(v) for v in SOILPROP_LL[:NLAYR]]
    MSALB = float(SOILPROP_MSALB)
    # SW may be scalar or list
    if hasattr(SW, "__len__") and not isinstance(SW, (str, bytes)):
        SWI = [float(v) for v in SW[:NLAYR]]
        if len(SWI) < NLAYR:
            SWI = SWI + [SWI[-1]] * (NLAYR - len(SWI))
    else:
        SWI = [float(SW)] * NLAYR

    _, DOY = YR_DOY(CONTROL_YRDOY)

    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    DLI_prev = 0.0
    for L in range(NLAYR):
        if L == 0:
            DLI = DS[L]
        else:
            DLI = DS[L] - DS[L - 1]
        DSMID[L] = CUMDPT + DLI * 5.0  # mm to mid-point
        CUMDPT = CUMDPT + DLI * 10.0   # mm profile depth
        TBD += BD[L] * DLI
        TLL += LL[L] * DLI
        TSW += SWI[L] * DLI
        TDL += DUL[L] * DLI
        DLI_prev = DLI

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * __import__('math').exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = __import__('math').log(500.0 / DP)
    ALBEDO = MSALB  # not used in SOILT but retained for completeness

    TMA = [0.0] * 5
    for i in range(5):
        TMA[i] = _nint(TAVG * 10000.0) / 10000.0
    ATOT = TMA[0] * 5.0

    ST = [float(TAVG)] * NLAYR
    SRFTEMP = float(TAVG)

    # Spin-up calls
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
          HDAY,
          CUMDPT,
          DSMID,
          TMA,
          ATOT):
    NLAYR = int(SOILPROP_NLAYR)
    BD = [float(v) for v in SOILPROP_BD[:NLAYR]]
    DLAYR = [float(v) for v in SOILPROP_DLAYR[:NLAYR]]
    DS = [float(v) for v in SOILPROP_DS[:NLAYR]]
    DUL = [float(v) for v in SOILPROP_DUL[:NLAYR]]
    LL = [float(v) for v in SOILPROP_LL[:NLAYR]]
    MSALB = float(SOILPROP_MSALB)
    if hasattr(SW, "__len__") and not isinstance(SW, (str, bytes)):
        SW_arr = [float(v) for v in SW[:NLAYR]]
        if len(SW_arr) < NLAYR:
            SW_arr = SW_arr + [SW_arr[-1]] * (NLAYR - len(SW_arr))
    else:
        SW_arr = [float(SW)] * NLAYR

    _, DOY = YR_DOY(CONTROL_YRDOY)

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL_local = 0.0
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL_local += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW_arr[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * __import__('math').exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = __import__('math').log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL_local - TLL)  # cm

    ATOT_new, TMA_new, SRFTEMP, ST = SOILT(
        ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
        PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
        ATOT, TMA
    )

    return CUMDPT, DSMID, TDL_local, TMA_new, ATOT_new, SRFTEMP, ST


def test_STEMP_example():
    # Example derived from ASKEE.for
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

    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY = STEMP_Initialize(
        CONTROL_YRDOY, ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
        SOILPROP_NLAYR, SOILPROP_MSALB, SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP
    )

    assert len(DSMID) == NLAYR
    assert len(ST) == NLAYR
    assert len(TMA) == 5

    CUMDPT2, DSMID2, TDL2, TMA2, ATOT2, SRFTEMP2, ST2 = STEMP(
        CONTROL_YRDOY, ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL,
        SOILPROP_LL, SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, TAV, TAMP, HDAY,
        CUMDPT, DSMID, TMA, ATOT
    )

    assert CUMDPT2 == CUMDPT
    assert DSMID2 == DSMID
    assert len(ST2) == NLAYR
    assert len(TMA2) == 5