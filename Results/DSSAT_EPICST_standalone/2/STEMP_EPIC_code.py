def STEMP_EPIC_Initialize(
    iswitch_ISWWAT,
    soilprop_BD,
    soilprop_DLAYR,
    soilprop_DS,
    soilprop_DUL,
    soilprop_LL,
    soilprop_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    orgc_MULCHMASS,
    water_SNOW,
    X2_PREV
):
    """
    Seasonal initialization for EPIC soil temperature routine.
    Inputs:
    - iswitch_ISWWAT: 'Y' or 'N'
    - soilprop_BD: list of bulk density per layer (g/cm3)
    - soilprop_DLAYR: list of layer thickness (cm)
    - soilprop_DS: list of cumulative depth by layer (cm)
    - soilprop_DUL: list of drained upper limit per layer (cm3/cm3)
    - soilprop_LL: list of lower limit per layer (cm3/cm3)
    - soilprop_NLAYR: number of soil layers (int)
    - SW: list of initial soil water content per layer (cm3/cm3)
    - TAVG: average daily temperature (degC)
    - TMAX: maximum daily temperature (degC)
    - TMIN: minimum daily temperature (degC)
    - TAV: average annual soil temperature (degC)
    - orgc_MULCHMASS: surface mulch mass (kg/ha)
    - water_SNOW: snow cover (mm)
    - X2_PREV: previous surface temperature memory (degC)
    Returns:
    - CUMDPT: profile depth (mm)
    - DSMID: list of depth to midpoint of each layer (mm)
    - TDL: total water content at DUL integrated across layers (cm)
    - TMA: list (len 5) of previous 5 days X2 values (degC)
    - NDays: number of days tracked for wet day memory (int)
    - WetDay: list (len 30) with 0/1 flags of wet days
    - X2_PREV: updated previous surface temperature memory (degC)
    - SRFTEMP: surface temperature (degC)
    - ST: list of soil temperature by layer (degC)
    """
    NLAYR = soilprop_NLAYR
    # Initialize accumulators
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR

    # Set initial soil water (SWI) = provided SW
    SWI = SW

    for L in range(NLAYR):
        DSMID[L] = CUMDPT + soilprop_DLAYR[L] * 5.0  # mm depth to midpt
        CUMDPT = CUMDPT + soilprop_DLAYR[L] * 10.0   # mm profile depth
        TBD += soilprop_BD[L] * soilprop_DLAYR[L]
        TLL += soilprop_LL[L] * soilprop_DLAYR[L]
        TSW += SWI[L] * soilprop_DLAYR[L]
        TDL += soilprop_DUL[L] * soilprop_DLAYR[L]

    if iswitch_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / soilprop_DS[NLAYR - 1]  # g/cm3
    FX = ABD / (ABD + 686.0 * (2.718281828459045 ** (-5.63 * ABD)))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    import math
    B = math.log(500.0 / DP)

    # Initialize 5-day memory TMA with TAVG (4 decimal places in Fortran NINT emulation not critical)
    TMA = [round(TAVG, 4) for _ in range(5)]

    # Initialize soil temperature by layer to TAVG
    ST = [TAVG for _ in range(NLAYR)]

    # Initialize wet day memory
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Soil cover function
    CV = (orgc_MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    SNOW = water_SNOW
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG  # initialize; will be updated by SOILT_EPIC
    X2_AVG = 0.0  # placeholder

    # Spin-up calls
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
            B=B,
            BCV=BCV,
            CUMDPT=CUMDPT,
            DP=DP,
            DSMID=DSMID,
            NLAYR=NLAYR,
            PESW=PESW,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            TMIN=TMIN,
            WetDay=0,
            WFT=WFT,
            WW=WW,
            TMA=TMA,
            ST=ST,
            X2_PREV=X2_PREV
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    iswitch_ISWWAT,
    soilprop_BD,
    soilprop_DLAYR,
    soilprop_DS,
    soilprop_DUL,
    soilprop_LL,
    soilprop_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    weather_RAIN,
    mgmt_DEPIR,
    plant_BIOMAS,
    orgc_MULCHMASS,
    water_SNOW,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    NDays,
    WetDay,
    X2_PREV,
    SRFTEMP,
    ST
):
    """
    Daily main biophysical process for EPIC soil temperature routine.
    Inputs:
    - iswitch_ISWWAT: 'Y' or 'N'
    - soilprop_BD, soilprop_DLAYR, soilprop_DS, soilprop_DUL, soilprop_LL: layer properties
    - soilprop_NLAYR: number of layers
    - SW: current soil water content by layer (cm3/cm3)
    - TAVG, TMAX, TMIN: daily temperatures (degC)
    - TAV: annual average air/soil temperature (degC)
    - weather_RAIN: daily rainfall (mm)
    - mgmt_DEPIR: daily irrigation depth (mm)
    - plant_BIOMAS: above-ground biomass (kg/ha)
    - orgc_MULCHMASS: mulch mass (kg/ha)
    - water_SNOW: snow cover (mm)
    State inputs (updated and returned):
    - CUMDPT, DSMID: geometry (mm)
    - TDL: total DUL water across profile (cm) [accumulates per original code]
    - TMA: previous 5-day X2 memory (degC)
    - NDays: number of days stored in WetDay (int, up to 30)
    - WetDay: list of last 30 day wet flags (0/1)
    - X2_PREV: previous surface temp memory (degC)
    - SRFTEMP: soil surface temperature (degC)
    - ST: soil temperature by layer (degC)
    Returns updated:
    - TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    """
    NLAYR = soilprop_NLAYR
    # Profile properties
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL is not reset here per original code; it accumulates across days if ISWWAT='N'
    for L in range(NLAYR):
        TBD += soilprop_BD[L] * soilprop_DLAYR[L]
        TDL += soilprop_DUL[L] * soilprop_DLAYR[L]
        TLL += soilprop_LL[L] * soilprop_DLAYR[L]
        TSW += SW[L] * soilprop_DLAYR[L]

    ABD = TBD / soilprop_DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * (2.718281828459045 ** (-5.63 * ABD)))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD   # vol. fraction
    import math
    B = math.log(500.0 / DP)

    if iswitch_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    # 30-day wet day memory
    RAIN = weather_RAIN
    DEPIR = mgmt_DEPIR
    if NDays == 30:
        # shift left
        for i in range(29):
            WetDay[i] = WetDay[i + 1]
        # NDays remains 30
    else:
        NDays = NDays + 1
    idx = NDays - 1 if NDays > 0 else 0
    WetDay[idx] = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    NWetDays = sum(WetDay[:NDays]) if NDays <= 30 else sum(WetDay)
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function
    BIOMAS = plant_BIOMAS
    MULCHMASS = orgc_MULCHMASS
    SNOW = water_SNOW
    CV = (BIOMAS + MULCHMASS) / 1000.0  # t/ha
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
    BCV = max(BCV1, BCV2)

    # Update temperatures via EPIC soil temperature subroutine
    TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
        B=B,
        BCV=BCV,
        CUMDPT=CUMDPT,
        DP=DP,
        DSMID=DSMID,
        NLAYR=NLAYR,
        PESW=PESW,
        TAV=TAV,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        WetDay=WetDay[idx],
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV
    )

    return TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def SOILT_EPIC(
    B,
    BCV,
    CUMDPT,
    DP,
    DSMID,
    NLAYR,
    PESW,
    TAV,
    TAVG,
    TMAX,
    TMIN,
    WetDay,
    WFT,
    WW,
    TMA,
    ST,
    X2_PREV
):
    """
    EPIC soil temperature calculation by layer.
    Inputs:
    - B: exponential decay factor
    - BCV: soil cover function (0-1)
    - CUMDPT: cumulative profile depth (mm)
    - DP: damping depth (mm)
    - DSMID: list of depth to midpoint of each layer (mm)
    - NLAYR: number of layers
    - PESW: potential extractable soil water over profile (cm)
    - TAV: average annual temperature (degC)
    - TAVG, TMAX, TMIN: daily temperatures (degC)
    - WetDay: wet day flag (0/1)
    - WFT: fraction wet days over last 30 days
    - WW: soil water content at saturation (volumetric fraction function of BD)
    - TMA: previous 5-day X2 array (degC) [InOut]
    - ST: soil temperature by layer (degC) [InOut]
    - X2_PREV: previous surface temp memory (degC) [InOut]
    Returns:
    - TMA: updated
    - SRFTEMP: surface temperature (degC)
    - ST: updated per layer temps
    - X2_AVG: average of last 5 X2 (degC)
    - X2_PREV: updated previous surface temp memory (degC) = X2_AVG
    """
    import math

    # Water content factor (dimensionless)
    denom = WW * CUMDPT
    if denom <= 0.0:
        WC = 0.01  # guard
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # EPIC surface driver X2
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update TMA memory: insert at front and shift
    TMA = [X2] + TMA[:4]
    X2_AVG = sum(TMA) / 5.0

    # Blend with cover effect
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    SRFTEMP = min(X2_AVG, X3)
    X1 = TAV - X3

    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        ST[L] = LAG * ST[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return TMA, SRFTEMP, ST, X2_AVG, X2_PREV