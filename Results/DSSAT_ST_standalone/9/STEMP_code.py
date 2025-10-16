def _nint_round(x, decimals):
    import math
    factor = 10.0 ** decimals
    if x >= 0:
        return math.floor(x * factor + 0.5) / factor
    else:
        return math.ceil(x * factor - 0.5) / factor


def YR_DOY(YRDOY):
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA, SRFTEMP, ST):
    import math

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT = ATOT - TMA[4]

    # Shift TMA: K=5..2 becomes 4..1
    TMA = [TMA[0]] + TMA[:4]
    # Insert today's average temperature
    TMA[0] = TAVG
    TMA[0] = _nint_round(TMA[0], 4)
    ATOT = ATOT + TMA[0]

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    ST_out = list(ST)
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_out[L] = _nint_round(ST_val, 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, ST_out


def STEMP_Initialize(
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
    TAMP,
    CONTROL_YRDOY,
):
    import math

    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD[:NLAYR])
    DSI = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    SWI = list(SW[:NLAYR])
    MSALB = SOILPROP_MSALB

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
    DLI = [0.0] * NLAYR

    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DSI[L]
        else:
            DLI[L] = DSI[L] - DSI[L - 1]
        DSMID[L] = CUMDPT + DLI[L] * 5.0
        CUMDPT = CUMDPT + DLI[L] * 10.0
        TBD = TBD + BD[L] * DLI[L]
        TLL = TLL + LL[L] * DLI[L]
        TSW = TSW + SWI[L] * DLI[L]
        TDL = TDL + DUL[L] * DLI[L]

    if ISWITCH_ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    TMA = [_nint_round(TAVG, 4) for _ in range(5)]
    ATOT = TMA[0] * 5.0
    ST = [TAVG for _ in range(NLAYR)]
    SRFTEMP = TAVG

    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            ALBEDO,
            B,
            CUMDPT,
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
            DSMID,
            ATOT,
            TMA,
            SRFTEMP,
            ST,
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def STEMP(
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
    CONTROL_YRDOY,
    HDAY,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    ATOT,
    ST,
):
    import math

    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    SW_arr = list(SW[:NLAYR])
    MSALB = SOILPROP_MSALB

    _, DOY = YR_DOY(CONTROL_YRDOY)

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD = TBD + BD[L] * DLAYR[L]
        TDL = TDL + DUL[L] * DLAYR[L]
        TLL = TLL + LL[L] * DLAYR[L]
        TSW = TSW + SW_arr[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT == "Y":
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT, TMA, SRFTEMP, ST = SOILT(
        ALBEDO,
        B,
        CUMDPT,
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
        DSMID,
        ATOT,
        TMA,
        0.0,
        ST,
    )

    return TDL, TMA, ATOT, SRFTEMP, ST