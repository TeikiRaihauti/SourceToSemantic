import math


def YR_DOY(YRDOY):
    # Converts YRDOY (YYYYDDD) to (YEAR, DOY)
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def _fortran_nint(x):
    # Fortran NINT: round to nearest integer, half away from zero
    if x >= 0.0:
        return math.floor(x + 0.5)
    else:
        return math.ceil(x - 0.5)


def _fortran_round(x, decimals):
    scale = 10.0 ** decimals
    return _fortran_nint(x * scale) / scale


def SOILT(
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
):
    # Copy arrays to avoid mutating inputs
    DSMID_loc = list(DSMID[:NLAYR])
    TMA_loc = list(TMA[:5])
    ST_loc = list(ST[:NLAYR])

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT_loc = ATOT - TMA_loc[4]

    # Shift TMA
    for K in range(4, 0, -1):
        TMA_loc[K] = TMA_loc[K - 1]

    TMA_loc[0] = TAVG
    TMA_loc[0] = _fortran_round(TMA_loc[0], 4)
    ATOT_loc = ATOT_loc + TMA_loc[0]

    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_loc / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID_loc[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_loc[L] = _fortran_round(ST_val, 3)

    SRFTEMP_loc = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_loc, TMA_loc, SRFTEMP_loc, ST_loc


def STEMP(
    CONTROL_DYNAMIC,
    CONTROL_YRDOY,
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
    XLAT,
    TAV,
    TAMP,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    ATOT,
    SRFTEMP,
    ST,
):
    # Local copies and sizing
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SW_arr = list(SW[:NLAYR])

    # State arrays
    DSMID_state = list(DSMID[:NLAYR]) if len(DSMID) >= NLAYR else [0.0] * NLAYR
    ST_state = list(ST[:NLAYR]) if len(ST) >= NLAYR else [0.0] * NLAYR
    TMA_state = list(TMA[:5]) if len(TMA) >= 5 else [TAVG] * 5

    _, DOY = YR_DOY(CONTROL_YRDOY)

    if XLAT < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    if CONTROL_DYNAMIC == 2:
        # SEASINIT
        SWI = list(SW_arr)
        DSI = list(DS)

        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_state = 0.0
        CUMDPT_state = 0.0
        DLI = [0.0] * NLAYR

        for L in range(NLAYR):
            if L == 0:
                DLI[L] = DSI[L]
            else:
                DLI[L] = DSI[L] - DSI[L - 1]
            DSMID_state[L] = CUMDPT_state + DLI[L] * 5.0
            CUMDPT_state = CUMDPT_state + DLI[L] * 10.0
            TBD = TBD + BD[L] * DLI[L]
            TLL = TLL + LL[L] * DLI[L]
            TSW = TSW + SWI[L] * DLI[L]
            TDL_state = TDL_state + DUL[L] * DLI[L]

        if ISWITCH_ISWWAT == "Y":
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_state - TLL)

        ABD = TBD / DSI[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = MSALB

        TMA_state = [ _fortran_round(TAVG, 4) for _ in range(5) ]
        ATOT_state = TMA_state[0] * 5.0

        for L in range(NLAYR):
            ST_state[L] = TAVG

        for _ in range(8):
            ATOT_state, TMA_state, SRFTEMP_state, ST_state = SOILT(
                ALBEDO,
                B,
                CUMDPT_state,
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
                DSMID_state,
                ATOT_state,
                TMA_state,
                SRFTEMP,
                ST_state,
            )

        return (
            CUMDPT_state,
            DSMID_state,
            TDL_state,
            TMA_state,
            ATOT_state,
            SRFTEMP_state,
            ST_state,
        )

    elif CONTROL_DYNAMIC == 3:
        # RATE
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_state = 0.0

        for L in range(NLAYR):
            TBD = TBD + BD[L] * DLAYR[L]
            TDL_state = TDL_state + DUL[L] * DLAYR[L]
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
            PESW = max(0.0, TDL_state - TLL)

        ATOT_state, TMA_state, SRFTEMP_state, ST_state = SOILT(
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
            DSMID_state,
            ATOT,
            TMA_state,
            SRFTEMP,
            ST_state,
        )

        return (
            CUMDPT,
            DSMID_state,
            TDL_state,
            TMA_state,
            ATOT_state,
            SRFTEMP_state,
            ST_state,
        )

    else:
        # For other dynamics stages (e.g., OUTPUT, SEASEND), no biophysical updates
        return (CUMDPT, DSMID_state, TDL, TMA_state, ATOT, SRFTEMP, ST_state)


def test_STEMP_basic_from_ASKEE():
    # Construct inputs similar to ASKEE.for
    CONTROL_DYNAMIC = 2  # SEASINIT
    CONTROL_YRDOY = 2021100

    ISWITCH_ISWWAT = "Y"

    SOILPROP_NLAYR = 4
    SOILPROP_BD = [1.6] * SOILPROP_NLAYR
    SOILPROP_DLAYR = [10.0] * SOILPROP_NLAYR
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * SOILPROP_NLAYR
    SOILPROP_LL = [0.2] * SOILPROP_NLAYR
    SOILPROP_MSALB = 0.13

    SRAD = 20.0
    SW = [0.2] * SOILPROP_NLAYR
    TAVG = 25.0
    TMAX = 30.0
    XLAT = 28.0
    TAV = 20.0
    TAMP = 10.0

    CUMDPT = 0.0
    DSMID = [0.0] * SOILPROP_NLAYR
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * SOILPROP_NLAYR

    # Initialization
    (
        CUMDPT,
        DSMID,
        TDL,
        TMA,
        ATOT,
        SRFTEMP,
        ST,
    ) = STEMP(
        CONTROL_DYNAMIC,
        CONTROL_YRDOY,
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
        XLAT,
        TAV,
        TAMP,
        CUMDPT,
        DSMID,
        TDL,
        TMA,
        ATOT,
        SRFTEMP,
        ST,
    )

    # One daily RATE step
    CONTROL_DYNAMIC = 3
    (
        CUMDPT_rate,
        DSMID_rate,
        TDL_rate,
        TMA_rate,
        ATOT_rate,
        SRFTEMP_rate,
        ST_rate,
    ) = STEMP(
        CONTROL_DYNAMIC,
        CONTROL_YRDOY,
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
        XLAT,
        TAV,
        TAMP,
        CUMDPT,
        DSMID,
        TDL,
        TMA,
        ATOT,
        SRFTEMP,
        ST,
    )

    # Return values for inspection by a testing framework
    return {
        "init": {
            "CUMDPT": CUMDPT,
            "DSMID": DSMID,
            "TDL": TDL,
            "TMA": TMA,
            "ATOT": ATOT,
            "SRFTEMP": SRFTEMP,
            "ST": ST,
        },
        "rate": {
            "CUMDPT": CUMDPT_rate,
            "DSMID": DSMID_rate,
            "TDL": TDL_rate,
            "TMA": TMA_rate,
            "ATOT": ATOT_rate,
            "SRFTEMP": SRFTEMP_rate,
            "ST": ST_rate,
        },
    }