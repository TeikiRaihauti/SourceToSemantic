def YR_DOY(YRDOY):
    """
    Converts YRDOY (YYYYDDD) to YEAR and DOY.
    Returns (YEAR, DOY).
    """
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA, SRFTEMP, ST):
    """
    Soil temperature by layer.

    Inputs:
    - ALBEDO (not used in calculations, kept for signature compatibility)
    - B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID
    InOut:
    - ATOT, TMA
    Outputs:
    - SRFTEMP, ST

    Returns updated (ATOT, TMA_new, SRFTEMP_new, ST_new)
    """
    import math

    # Make copies to preserve function purity
    TMA_new = list(TMA)
    ST_new = list(ST)

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT_new = ATOT - TMA_new[4]

    # Shift TMA history
    for K in range(4, 0, -1):
        TMA_new[K] = TMA_new[K - 1]

    TMA_new[0] = TAVG
    # Keep only 4 decimals
    TMA_new[0] = round(TMA_new[0], 4)
    ATOT_new = ATOT_new + TMA_new[0]

    # Corrected water content equation (EPIC)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # Damping depth (mm)

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_new / 5.0 - TA

    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_new[L] = round(ST_val, 3)

    SRFTEMP_new = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_new, TMA_new, SRFTEMP_new, ST_new


def STEMP(CONTROL_DYNAMIC, CONTROL_YRDOY,
          ISWITCH_ISWWAT,
          SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
          SOILPROP_NLAYR, SOILPROP_MSALB,
          SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
          CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY):
    """
    Determines soil temperature by layer.

    This function encapsulates both seasonal initialization (CONTROL_DYNAMIC == 2)
    and daily rate calculations (CONTROL_DYNAMIC == 3). State variables must be
    provided and updated values are returned.

    Inputs:
    - CONTROL_DYNAMIC: integer dynamic flag (2 = SEASINIT, 3 = RATE)
    - CONTROL_YRDOY: integer YYYYDDD
    - ISWITCH_ISWWAT: 'Y' if water simulated, 'N' otherwise
    - SOILPROP_BD: list of bulk density per layer (g/cm3), length NLAYR
    - SOILPROP_DLAYR: list of layer thickness (cm), length NLAYR
    - SOILPROP_DS: list of cumulative depth (cm) to bottom of each layer, length NLAYR
    - SOILPROP_DUL: list of drained upper limit (cm3/cm3), length NLAYR
    - SOILPROP_LL: list of lower limit (cm3/cm3), length NLAYR
    - SOILPROP_NLAYR: number of layers
    - SOILPROP_MSALB: soil/mulch albedo (dimensionless)
    - SRAD: solar radiation (MJ/m2-d) [not used in calculations, kept for compatibility]
    - SW: list of volumetric water content by layer (cm3/cm3), length NLAYR
    - TAVG: average daily temperature (C)
    - TMAX: maximum daily temperature (C) [not used directly]
    - XLAT: latitude (deg)
    - TAV: average annual soil temperature (C)
    - TAMP: amplitude for soil temperature (C)
    - CUMDPT: cumulative soil profile depth (mm) [state]
    - DSMID: list of depth to midpoint of layers (mm), length NLAYR [state]
    - TDL: total water content at DUL (cm) [state]
    - TMA: list of last 5 days soil temperatures (C), length 5 [state]
    - ATOT: sum of TMA array (C) [state]
    - SRFTEMP: soil surface temperature (C) [state]
    - ST: list soil temperature by layer (C), length NLAYR [state]
    - HDAY: hottest day of year used in seasonal phase (DOY, float) [state]

    Returns updated tuple:
    (CUMDPT, DSMID_new, TDL_new, TMA_new, ATOT_new, SRFTEMP_new, ST_new, HDAY_new)
    """
    import math

    # Copy inputs to avoid mutation (purity)
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    SW_arr = list(SW)
    DSMID_new = list(DSMID)
    TMA_new = list(TMA)
    ST_new = list(ST)

    # Determine DOY from YRDOY
    _, DOY = YR_DOY(CONTROL_YRDOY)

    # Seasonal initialization
    if CONTROL_DYNAMIC == 2:
        SWI = list(SW_arr)
        DSI = list(DS)

        # Hottest day of year based on latitude
        if XLAT < 0.0:
            HDAY_new = 20.0   # Southern hemisphere
        else:
            HDAY_new = 200.0  # Northern hemisphere

        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_new = 0.0
        CUMDPT_new = 0.0

        # Compute DSMID (mm) and cumulative depth (mm); sums weighted by layer thickness (cm)
        for L in range(NLAYR):
            if L == 0:
                DLI = DSI[L]
            else:
                DLI = DSI[L] - DSI[L - 1]

            DSMID_new[L] = CUMDPT_new + DLI * 5.0
            CUMDPT_new = CUMDPT_new + DLI * 10.0

            TBD += BD[L] * DLI
            TLL += LL[L] * DLI
            TSW += SWI[L] * DLI
            TDL_new += DUL[L] * DLI

        if ISWITCH_ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)  # cm
        else:
            PESW = max(0.0, TDL_new - TLL)  # cm

        ABD = TBD / DSI[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX  # mm
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = SOILPROP_MSALB

        # Initialize TMA and ATOT
        TMA_val = round(TAVG, 4)
        TMA_new = [TMA_val for _ in range(5)]
        ATOT_new = TMA_new[0] * 5.0

        # Initialize ST to TAVG then spin-up
        for L in range(NLAYR):
            ST_new[L] = TAVG

        for _ in range(8):
            ATOT_new, TMA_new, SRFTEMP_new, ST_new = SOILT(
                ALBEDO, B, CUMDPT_new, DOY, DP, HDAY_new, NLAYR,
                PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_new,
                ATOT_new, TMA_new, SRFTEMP, ST_new
            )

        return (CUMDPT_new, DSMID_new, TDL_new, TMA_new, ATOT_new,
                SRFTEMP_new, ST_new, HDAY_new)

    # Daily rate calculations
    elif CONTROL_DYNAMIC == 3:
        # Use provided HDAY
        HDAY_new = HDAY

        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_new = TDL  # Note: not reset per original code

        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_new += DUL[L] * DLAYR[L]
            TLL += LL[L] * DLAYR[L]
            TSW += SW_arr[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX  # mm
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)
        ALBEDO = SOILPROP_MSALB

        if ISWITCH_ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)  # cm
        else:
            PESW = max(0.0, TDL_new - TLL)  # cm

        ATOT_new, TMA_new, SRFTEMP_new, ST_new = SOILT(
            ALBEDO, B, CUMDPT, DOY, DP, HDAY_new, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
            ATOT, TMA, SRFTEMP, ST
        )

        return (CUMDPT, DSMID, TDL_new, TMA_new, ATOT_new,
                SRFTEMP_new, ST_new, HDAY_new)

    else:
        # If dynamic phase not handled, return inputs unchanged
        return (CUMDPT, list(DSMID), TDL, list(TMA), ATOT, SRFTEMP, list(ST), HDAY)