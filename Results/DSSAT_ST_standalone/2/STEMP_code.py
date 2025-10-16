import math

def YR_DOY(YRDOY):
    """
    Converts YRDOY (YYYYDDD) to (YR, DOY).
    """
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def _nint(x):
    """
    Fortran-like NINT: nearest integer, halves away from zero.
    """
    if x >= 0.0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


def SOILT(
    ALBEDO,  # Input (unused in core equations, kept for interface consistency)
    B,       # Input
    CUMDPT,  # Input (mm)
    DOY,     # Input (day of year)
    DP,      # Input (mm)
    HDAY,    # Input
    NLAYR,   # Input
    PESW,    # Input (cm)
    SRAD,    # Input (unused in core equations, kept for interface consistency)
    TAMP,    # Input (degC)
    TAV,     # Input (degC)
    TAVG,    # Input (degC)
    TMAX,    # Input (degC, unused here)
    WW,      # Input (dimensionless)
    DSMID,   # Input array (mm), length NLAYR
    ATOT,    # InOut (sum of last 5 days of average soil temperatures)
    TMA,     # InOut array (5) of previous 5 days of average soil temperatures
    ST       # Output array (NLAYR) soil temperature by layer (degC)
):
    """
    Determines soil temperature by layer.
    Returns updated (ATOT, TMA, SRFTEMP, ST).
    """
    # Copy inputs to avoid side-effects
    TMA_out = list(TMA) if TMA is not None and len(TMA) >= 5 else [0.0] * 5
    ST_out = list(ST) if ST is not None and len(ST) >= NLAYR else [0.0] * NLAYR

    # Seasonal phase angle
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update moving average array and sum
    ATOT_out = ATOT - TMA_out[4]
    for k in range(4, 0, -1):
        TMA_out[k] = TMA_out[k - 1]
    TMA_out[0] = TAVG
    # Keep only 4 decimals (debug vs release fix)
    TMA_out[0] = _nint(TMA_out[0] * 10000.0) / 10000.0
    ATOT_out = ATOT_out + TMA_out[0]

    # Water content function (corrected equation)
    # WC (ratio) = max(0.01, PESW) / (WW * CUMDPT) * 10
    denom = WW * CUMDPT
    if denom <= 0.0:
        WC = 0.01  # fallback to minimum to avoid divide-by-zero
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)

    DD = FX * DP  # damping depth in mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_out / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD if DD != 0.0 else 0.0
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        # Keep only 3 decimals (debug vs release fix)
        ST_out[L] = _nint(ST_val * 1000.0) / 1000.0

    # Surface litter layer temperature
    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_out, TMA_out, SRFTEMP, ST_out


def STEMP(
    CONTROL_DYNAMIC,  # integer dynamic flag (e.g., 2=SEASINIT, 3=RATE)
    CONTROL_YRDOY,    # YYYYDDD
    CONTROL_RUN,      # integer RUN number
    CONTROL_RNMODE,   # character RNMODE
    ISWITCH_ISWWAT,   # 'Y' or 'N'
    SOILPROP_BD,      # array, bulk density (g/cm3), length NLAYR
    SOILPROP_DLAYR,   # array, layer thickness (cm), length NLAYR
    SOILPROP_DS,      # array, cumulative depth (cm), length NLAYR
    SOILPROP_DUL,     # array, DUL (cm3/cm3), length NLAYR
    SOILPROP_LL,      # array, LL (cm3/cm3), length NLAYR
    SOILPROP_NLAYR,   # integer, number of layers
    SOILPROP_MSALB,   # float, soil/mulch albedo
    SRAD,             # float, solar radiation (MJ/m2-d)
    SW,               # array, volumetric soil water content (cm3/cm3), length NLAYR
    TAVG,             # float, average daily air temperature (degC)
    TMAX,             # float, maximum daily air temperature (degC)
    XLAT,             # float, latitude (deg)
    TAV,              # float, average annual soil temperature (degC)
    TAMP,             # float, amplitude of temperature function (degC)
    CUMDPT,           # InOut float (mm), cumulative depth of soil profile
    DSMID,            # InOut array (mm), depth to midpoint of soil layer, length NLAYR
    TDL,              # InOut float (cm), total water at DUL
    TMA,              # InOut array (5), previous 5 days avg soil temps
    ATOT,             # InOut float, sum of TMA
    SRFTEMP,          # InOut float, soil surface temperature
    ST                # InOut array (degC), soil temperature by layer (NLAYR)
):
    """
    Determines soil temperature by layer across dynamic phases.
    Returns updated (CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST).
    """
    # Copy inputs to avoid mutating caller state
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = float(SOILPROP_MSALB)
    SW_arr = list(SW[:NLAYR])

    DSMID_out = list(DSMID[:NLAYR]) if DSMID is not None and len(DSMID) >= NLAYR else [0.0] * NLAYR
    TMA_out = list(TMA[:5]) if TMA is not None and len(TMA) >= 5 else [0.0] * 5
    ST_out = list(ST[:NLAYR]) if ST is not None and len(ST) >= NLAYR else [0.0] * NLAYR

    ATOT_out = float(ATOT)
    CUMDPT_out = float(CUMDPT)
    TDL_out = float(TDL)
    SRFTEMP_out = float(SRFTEMP)

    # Compute DOY from YRDOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Latitude based hottest day (HDAY)
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    if CONTROL_DYNAMIC == 2:
        # SEASINIT
        # Gating to avoid re-initialization for certain run modes
        do_init = (CONTROL_RUN == 1) or (CONTROL_RNMODE not in ('Q', 'F'))
        if do_init:
            SWI = SW_arr[:]  # Initial soil water content per layer
            DSI = DS[:]      # Cumulative depth (cm)

            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_out = 0.0
            CUMDPT_out = 0.0

            for L in range(NLAYR):
                if L == 0:
                    DLI = DSI[L]
                else:
                    DLI = DSI[L] - DSI[L - 1]
                # Depth to midpoint of layer (mm)
                DSMID_out[L] = CUMDPT_out + DLI * 5.0
                # Update total profile depth (mm)
                CUMDPT_out = CUMDPT_out + DLI * 10.0

                TBD += BD[L] * DLI
                TLL += LL[L] * DLI
                TSW += SWI[L] * DLI
                TDL_out += DUL[L] * DLI

            if ISWITCH_ISWWAT == 'Y':
                PESW = max(0.0, TSW - TLL)  # cm
            else:
                PESW = max(0.0, TDL_out - TLL)  # cm

            ABD = TBD / DSI[NLAYR - 1] if DSI[NLAYR - 1] != 0.0 else 0.0
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP) if DP != 0.0 else 0.0
            ALBEDO = MSALB

            # Initialize TMA with TAVG, four-decimal rounding
            for i in range(5):
                TMA_out[i] = _nint(TAVG * 10000.0) / 10000.0
            ATOT_out = TMA_out[0] * 5.0

            # Initialize ST with TAVG
            for L in range(NLAYR):
                ST_out[L] = TAVG

            # Spin-up SOILT 8 times
            for _ in range(8):
                ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
                    ALBEDO, B, CUMDPT_out, DOY, DP, HDAY, NLAYR,
                    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_out,
                    ATOT_out, TMA_out, ST_out
                )

    elif CONTROL_DYNAMIC == 3:
        # RATE (daily calculation)
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_out += DUL[L] * DLAYR[L]  # Note: No reset per original code
            TLL += LL[L] * DLAYR[L]
            TSW += SW_arr[L] * DLAYR[L]

        total_depth_cm = DS[NLAYR - 1] if NLAYR > 0 else 0.0
        ABD = TBD / total_depth_cm if total_depth_cm != 0.0 else 0.0
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
        DP = 1000.0 + 2500.0 * FX  # mm
        WW = 0.356 - 0.144 * ABD   # vol. fraction
        B = math.log(500.0 / DP) if DP != 0.0 else 0.0
        ALBEDO = MSALB

        if ISWITCH_ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)  # cm
        else:
            PESW = max(0.0, TDL_out - TLL)  # cm

        ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
            ALBEDO, B, CUMDPT_out, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_out,
            ATOT_out, TMA_out, ST_out
        )

    # OUTPUT and SEASEND branches (and any printing) are intentionally omitted.

    return CUMDPT_out, DSMID_out, TDL_out, TMA_out, ATOT_out, SRFTEMP_out, ST_out