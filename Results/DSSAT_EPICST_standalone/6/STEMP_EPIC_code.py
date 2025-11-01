def initialize_STEMP_EPIC(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    ORGC_MULCHMASS,
    WATER_SNOW,
    X2_PREV
):
    """
    Seasonal initialization for EPIC soil temperature routine.
    Returns initialized state variables and spun-up soil temperatures.

    Inputs:
      ISWWAT: 'Y' or 'N' to indicate water simulation switch
      SOILPROP_BD: list of bulk density per layer (g/cm3)
      SOILPROP_DLAYR: list of layer thickness (cm)
      SOILPROP_DS: list of cumulative depth (cm)
      SOILPROP_DUL: list of drained upper limit (cm3/cm3)
      SOILPROP_LL: list of lower limit (cm3/cm3)
      SOILPROP_NLAYR: integer number of layers
      SW: list of soil water content (cm3/cm3) initial
      TAVG, TMAX, TMIN: daily air temperatures (degC)
      TAV: average annual temperature (degC)
      ORGC_MULCHMASS: surface mulch mass (kg/ha)
      WATER_SNOW: snow water equivalent (mm)
      X2_PREV: previous day surface temperature driver (degC) state

    Returns tuple:
      CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    """
    import math

    # Local copy/aliases for readability consistent with original code
    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    NLAYR = int(SOILPROP_NLAYR)
    SWI = list(SW)

    # Helper to mimic Fortran NINT for positive/negative values
    def nint(x):
        return int(x + 0.5) if x >= 0.0 else -int(abs(x) + 0.5)

    # Initialize cumulative depths and midpoints
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    for L in range(NLAYR):
        DSMID[L] = CUMDPT + DLAYR[L] * 5.0  # mm depth to midpoint of layer
        CUMDPT = CUMDPT + DLAYR[L] * 10.0   # mm cumulative profile depth
        TBD = TBD + BD[L] * DLAYR[L]
        TLL = TLL + LL[L] * DLAYR[L]
        TSW = TSW + SWI[L] * DLAYR[L]
        TDL = TDL + DUL[L] * DLAYR[L]

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    # Initialize TMA with rounded TAVG as in Fortran (to 4 decimals)
    TMA = [nint(TAVG * 10000.0) / 10000.0 for _ in range(5)]
    X2_AVG = TMA[0] * 5.0

    # Initialize soil temperature by layer
    ST = [TAVG for _ in range(NLAYR)]

    # Wet day memory
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Soil cover function from mulch and snow
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Spin-up surface and soil temperatures for 8 iterations
    SRFTEMP = TAVG
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_AVG, X2_PREV = SOILT_EPIC(
            B=B, BCV=BCV, CUMDPT=CUMDPT, DP=DP, DSMID=DSMID, NLAYR=NLAYR, PESW=PESW, TAV=TAV,
            TAVG=TAVG, TMAX=TMAX, TMIN=TMIN, WetDay=0, WFT=WFT, WW=WW,
            TMA=TMA, SRFTEMP=SRFTEMP, ST=ST, X2_AVG=X2_AVG, X2_PREV=X2_PREV
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def STEMP_EPIC(
    ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_RAIN,
    MGMT_DEPIR,
    PLANT_BIOMAS,
    ORGC_MULCHMASS,
    WATER_SNOW,
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
    Daily main biophysical process for EPIC soil temperature routine (RATE step).

    Inputs:
      ISWWAT: 'Y' or 'N'
      SOILPROP_*: soil property arrays and number of layers
      SW: current soil water content by layer (cm3/cm3)
      TAVG, TMAX, TMIN: daily air temperatures (degC)
      TAV: average annual air temperature (degC)
      WEATHER_RAIN: rainfall (mm)
      MGMT_DEPIR: irrigation (mm)
      PLANT_BIOMAS: plant biomass (kg/ha)
      ORGC_MULCHMASS: mulch mass (kg/ha)
      WATER_SNOW: snow (mm)
      CUMDPT: cumulative depth (mm)
      DSMID: list of mid-layer depths (mm)
      TDL: cumulative DUL thickness state (cm) [accumulated as in original code]
      TMA: last 5 days driver temperatures (degC)
      NDays: wet-day memory current length (days)
      WetDay: list of last 30 wet-day flags
      X2_PREV: previous day surface driver temperature (degC)
      SRFTEMP: previous surface temperature (degC)
      ST: soil temperature by layer (degC)

    Returns tuple of updated states:
      TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    """
    import math

    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    NLAYR = int(SOILPROP_NLAYR)
    SW_loc = list(SW)
    TMA_loc = list(TMA)
    ST_loc = list(ST)
    WetDay_loc = list(WetDay)

    # Profile sums
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for L in range(NLAYR):
        TBD = TBD + BD[L] * DLAYR[L]
        TDL = TDL + DUL[L] * DLAYR[L]  # Note: not reset daily in original code
        TLL = TLL + LL[L] * DLAYR[L]
        TSW = TSW + SW_loc[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)

    if ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    # Wet day memory update
    RAIN = WEATHER_RAIN
    DEPIR = MGMT_DEPIR
    if NDays == 30:
        for I in range(29):
            WetDay_loc[I] = WetDay_loc[I + 1]
    else:
        NDays = NDays + 1
    WetDay_loc_length = NDays  # last filled index is NDays-1
    WetDay_loc[NDays - 1] = 1 if (RAIN + DEPIR) > 1.0e-6 else 0
    NWetDays = sum(WetDay_loc[:NDays])
    WFT = float(NWetDays) / float(NDays) if NDays > 0 else 0.0

    # Soil cover function from plant + mulch and snow
    BIOMAS = PLANT_BIOMAS
    MULCHMASS = ORGC_MULCHMASS
    SNOW = WATER_SNOW
    CV = (BIOMAS + MULCHMASS) / 1000.0
    BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
    BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Advance soil temperatures
    X2_AVG = sum(TMA_loc) / 5.0 if len(TMA_loc) == 5 else 0.0
    TMA_loc, SRFTEMP, ST_loc, X2_AVG, X2_PREV = SOILT_EPIC(
        B=B, BCV=BCV, CUMDPT=CUMDPT, DP=DP, DSMID=DSMID, NLAYR=NLAYR, PESW=PESW, TAV=TAV,
        TAVG=TAVG, TMAX=TMAX, TMIN=TMIN, WetDay=WetDay_loc[NDays - 1], WFT=WFT, WW=WW,
        TMA=TMA_loc, SRFTEMP=SRFTEMP, ST=ST_loc, X2_AVG=X2_AVG, X2_PREV=X2_PREV
    )

    return TDL, TMA_loc, NDays, WetDay_loc, X2_PREV, SRFTEMP, ST_loc


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
    SRFTEMP,
    ST,
    X2_AVG,
    X2_PREV
):
    """
    Supporting function: Determines soil temperature by layer using EPIC method.

    Inputs:
      B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV, TAVG, TMAX, TMIN, WetDay, WFT, WW
      TMA: list length 5
      SRFTEMP: previous surface temperature
      ST: list of soil temperatures by layer
      X2_AVG: previous X2 average
      X2_PREV: previous X2 average stored for cover blending

    Returns updated:
      TMA, SRFTEMP, ST, X2_AVG, X2_PREV
    """
    import math

    # Local copies
    TMA_loc = list(TMA)
    ST_loc = list(ST)

    # Water content factor (ratio), guard against divide-by-zero
    denom = WW * CUMDPT
    if denom <= 0.0:
        denom = 1.0e-6
    WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Surface driver temperature (Potter & Williams 1994)
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory TMA: shift and insert
    # Fortran: TMA(1)=X2; then K=5..2: TMA(K)=TMA(K-1)
    TMA_loc = [X2, TMA_loc[0], TMA_loc[1], TMA_loc[2], TMA_loc[3]]

    X2_AVG = sum(TMA_loc) / 5.0

    # Blend with cover and previous value
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    # Surface temperature is min(X2_AVG, X3)
    SRFTEMP = min(X2_AVG, X3)

    # Temperature difference driver
    X1 = TAV - X3

    # Layer temperature update with lag
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD if DD != 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) if (ZD + math.exp(-0.8669 - 2.0775 * ZD)) != 0.0 else 0.0
        ST_loc[L] = LAG * ST_loc[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV = X2_AVG

    return TMA_loc, SRFTEMP, ST_loc, X2_AVG, X2_PREV