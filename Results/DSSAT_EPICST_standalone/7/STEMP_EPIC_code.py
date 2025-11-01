from math import exp, log
from typing import List, Sequence, Tuple


def SOILT_EPIC(
    B: float,
    BCV: float,
    CUMDPT: float,
    DP: float,
    DSMID: Sequence[float],
    NLAYR: int,
    PESW: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    TMIN: float,
    WetDay: int,
    WFT: float,
    WW: float,
    TMA: Sequence[float],
    ST: Sequence[float],
    X2_PREV: float,
) -> Tuple[List[float], float, List[float], float]:
    """
    EPIC soil temperature by layer.

    Inputs
    - B: float, log(500/DP).
    - BCV: float, soil cover factor (0-1) from biomass/mulch/snow.
    - CUMDPT: float, cumulative profile depth (mm).
    - DP: float, damping depth (mm).
    - DSMID: list[float], depth to mid-point of each layer (mm), length >= NLAYR.
    - NLAYR: int, number of soil layers.
    - PESW: float, plant extractable soil water in profile (cm).
    - TAV: float, long-term mean air temperature (C).
    - TAVG: float, daily average air temperature (C).
    - TMAX: float, daily maximum air temperature (C).
    - TMIN: float, daily minimum air temperature (C).
    - WetDay: int, 1 for wet day (rain+irrigation), else 0.
    - WFT: float, fraction of wet days over memory window (0-1).
    - WW: float, soil water content at saturation (volumetric fraction).
    - TMA: list[float], 5-day memory state (size 5).
    - ST: list[float], soil temperature by layer (C), size >= NLAYR.
    - X2_PREV: float, previous day's X2_AVG state.

    Returns
    - TMA_new: list[float], updated 5-day memory state.
    - SRFTEMP: float, surface soil temperature (C).
    - ST_new: list[float], updated soil temperature by layer (C).
    - X2_PREV_new: float, updated X2_PREV (= current X2_AVG).
    """
    # Guard against division by zero
    denom = WW * CUMDPT
    if denom <= 0.0:
        denom = 1e-6

    WC = max(0.01, PESW) / denom * 10.0  # ratio

    FX = exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    # Compute X2 based on wetness of the day
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory TMA: shift and insert X2
    TMA_new = [0.0] * 5
    TMA_new[0] = X2
    for k in range(4, 0, -1):
        TMA_new[k] = TMA[k - 1]
    X2_AVG = sum(TMA_new) / 5.0

    # Blend with cover
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV

    # Surface temperature
    SRFTEMP = min(X2_AVG, X3)

    # Temperature difference toward mean
    X1 = TAV - X3

    # Layer temperatures
    LAG = 0.5
    ST_new = list(ST)
    for l in range(NLAYR):
        ZD = DSMID[l] / max(DD, 1e-6)
        F = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        ST_new[l] = LAG * ST[l] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_new = X2_AVG
    return TMA_new, SRFTEMP, ST_new, X2_PREV_new


def stemp_epic_initialize(
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: Sequence[float],
    SOILPROP_DLAYR: Sequence[float],
    SOILPROP_DS: Sequence[float],
    SOILPROP_DUL: Sequence[float],
    SOILPROP_LL: Sequence[float],
    SOILPROP_NLAYR: int,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    ORGC_MULCHMASS: float,
    WATER_SNOW: float,
    X2_PREV: float = 0.0,
) -> Tuple[
    float, List[float], float, List[float], int, List[int], float, float, List[float]
]:
    """
    Seasonal initialization for EPIC soil temperature.

    Inputs
    - ISWITCH_ISWWAT: str, 'Y' if water balance simulated, else 'N'.
    - SOILPROP_BD: list[float], bulk density by layer (g/cm3).
    - SOILPROP_DLAYR: list[float], layer thickness (cm).
    - SOILPROP_DS: list[float], cumulative depth at bottom of layers (cm).
    - SOILPROP_DUL: list[float], drained upper limit (cm3/cm3).
    - SOILPROP_LL: list[float], lower limit (cm3/cm3).
    - SOILPROP_NLAYR: int, number of layers.
    - SW: list[float], current soil water content by layer (cm3/cm3).
    - TAVG: float, daily average air temperature (C) at initialization.
    - TMAX: float, daily maximum air temperature (C).
    - TMIN: float, daily minimum air temperature (C).
    - TAV: float, long-term mean air temperature (C).
    - ORGC_MULCHMASS: float, surface mulch mass (kg/ha).
    - WATER_SNOW: float, snow water equivalent (mm).
    - X2_PREV: float, previous X2_AVG state (default 0).

    Returns
    - CUMDPT: float, cumulative profile depth (mm).
    - DSMID: list[float], depth to mid-point of each layer (mm).
    - TDL: float, sum(DUL*DLAYR) (cm).
    - TMA: list[float], initialized 5-day memory state.
    - NDays: int, initialized number of days in wet-day memory (0).
    - WetDay: list[int], wet-day memory (length 30) initialized to zeros.
    - X2_PREV: float, updated previous X2_AVG after spin-up.
    - SRFTEMP: float, initialized surface temperature (C).
    - ST: list[float], initialized soil temperature by layer (C).
    """
    NLAYR = SOILPROP_NLAYR
    # Initialize cumulative depth and midpoints (mm)
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0

    for l in range(NLAYR):
        DSMID[l] = CUMDPT + SOILPROP_DLAYR[l] * 5.0  # mm
        CUMDPT += SOILPROP_DLAYR[l] * 10.0  # mm
        TBD += SOILPROP_BD[l] * SOILPROP_DLAYR[l]
        TLL += SOILPROP_LL[l] * SOILPROP_DLAYR[l]
        TSW += SW[l] * SOILPROP_DLAYR[l]
        TDL += SOILPROP_DUL[l] * SOILPROP_DLAYR[l]

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / max(SOILPROP_DS[NLAYR - 1], 1e-6)
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / max(DP, 1e-6))

    # Initialize memory TMA and layer temps
    init_t = round(TAVG, 4)
    TMA = [init_t for _ in range(5)]
    ST = [TAVG for _ in range(NLAYR)]

    # Wet-day memory
    WFT = 0.1
    WetDay = [0 for _ in range(30)]
    NDays = 0

    # Soil cover function (init uses mulch only as in original)
    CV = ORGC_MULCHMASS / 1000.0
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if CV > 0.0 else 0.0
    SNOW = WATER_SNOW
    BCV2 = SNOW / (SNOW + exp(2.303 - 0.2197 * SNOW)) if SNOW > 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    SRFTEMP = TAVG
    # Spin-up calls
    for _ in range(8):
        TMA, SRFTEMP, ST, X2_PREV = SOILT_EPIC(
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
            X2_PREV=X2_PREV,
        )

    return CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST


def stemp_epic_process(
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: Sequence[float],
    SOILPROP_DLAYR: Sequence[float],
    SOILPROP_DS: Sequence[float],
    SOILPROP_DUL: Sequence[float],
    SOILPROP_LL: Sequence[float],
    SOILPROP_NLAYR: int,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    WEATHER_RAIN: float,
    MGMT_DEPIR: float,
    PLANT_BIOMAS: float,
    ORGC_MULCHMASS: float,
    WATER_SNOW: float,
    CUMDPT: float,
    DSMID: Sequence[float],
    TDL: float,
    TMA: Sequence[float],
    NDays: int,
    WetDay: Sequence[int],
    X2_PREV: float,
    ST: Sequence[float],
) -> Tuple[float, List[float], int, List[int], float, float, List[float]]:
    """
    Daily EPIC soil temperature process.

    Inputs
    - ISWITCH_ISWWAT: str, 'Y' if water balance simulated, else 'N'.
    - SOILPROP_BD: list[float], bulk density by layer (g/cm3).
    - SOILPROP_DLAYR: list[float], layer thickness (cm).
    - SOILPROP_DS: list[float], cumulative depth at bottom of layers (cm).
    - SOILPROP_DUL: list[float], drained upper limit (cm3/cm3).
    - SOILPROP_LL: list[float], lower limit (cm3/cm3).
    - SOILPROP_NLAYR: int, number of layers.
    - SW: list[float], current soil water content by layer (cm3/cm3).
    - TAVG: float, daily average air temperature (C).
    - TMAX: float, daily maximum air temperature (C).
    - TMIN: float, daily minimum air temperature (C).
    - TAV: float, long-term mean air temperature (C).
    - WEATHER_RAIN: float, rainfall (mm).
    - MGMT_DEPIR: float, irrigation depth (mm).
    - PLANT_BIOMAS: float, standing biomass (kg/ha).
    - ORGC_MULCHMASS: float, surface mulch mass (kg/ha).
    - WATER_SNOW: float, snow water equivalent (mm).
    - CUMDPT: float, cumulative profile depth (mm), from initialization.
    - DSMID: list[float], depth to mid-point of each layer (mm).
    - TDL: float, prior sum(DUL*DLAYR) (cm). Note: accumulated per original code.
    - TMA: list[float], 5-day memory state (size 5).
    - NDays: int, number of days recorded in wet-day memory (<=30).
    - WetDay: list[int], wet-day memory (length 30).
    - X2_PREV: float, previous X2_AVG state.
    - ST: list[float], prior soil temperature by layer (C).

    Returns
    - TDL: float, updated (accumulated) sum(DUL*DLAYR) per original logic.
    - TMA: list[float], updated 5-day memory state.
    - NDays: int, updated number of days in wet-day memory.
    - WetDay: list[int], updated wet-day memory.
    - X2_PREV: float, updated previous X2_AVG state.
    - SRFTEMP: float, surface soil temperature (C).
    - ST: list[float], updated soil temperature by layer (C).
    """
    NLAYR = SOILPROP_NLAYR

    # Bulk properties and water sums
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL is intentionally NOT reset per original code
    for l in range(NLAYR):
        TBD += SOILPROP_BD[l] * SOILPROP_DLAYR[l]
        TDL += SOILPROP_DUL[l] * SOILPROP_DLAYR[l]
        TLL += SOILPROP_LL[l] * SOILPROP_DLAYR[l]
        TSW += SW[l] * SOILPROP_DLAYR[l]

    ABD = TBD / max(SOILPROP_DS[NLAYR - 1], 1e-6)
    FX = ABD / (ABD + 686.0 * exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = log(500.0 / max(DP, 1e-6))

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    # Update 30-day wet-day memory
    WetDay_list = list(WetDay)
    if NDays == 30:
        for i in range(29):
            WetDay_list[i] = WetDay_list[i + 1]
    else:
        NDays += 1
    WetDay_list[NDays - 1] = 1 if (WEATHER_RAIN + MGMT_DEPIR) > 1.0e-6 else 0

    NWetDays = sum(WetDay_list[:NDays])
    WFT = float(NWetDays) / float(max(NDays, 1))

    # Soil cover function (biomass + mulch, or snow)
    CV = (PLANT_BIOMAS + ORGC_MULCHMASS) / 1000.0
    BCV1 = CV / (CV + exp(5.3396 - 2.3951 * CV)) if CV > 0.0 else 0.0
    SNOW = WATER_SNOW
    BCV2 = SNOW / (SNOW + exp(2.303 - 0.2197 * SNOW)) if SNOW > 0.0 else 0.0
    BCV = max(BCV1, BCV2)

    # Temperature evolution
    TMA_new, SRFTEMP, ST_new, X2_PREV_new = SOILT_EPIC(
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
        WetDay=WetDay_list[NDays - 1],
        WFT=WFT,
        WW=WW,
        TMA=TMA,
        ST=ST,
        X2_PREV=X2_PREV,
    )

    return TDL, TMA_new, NDays, WetDay_list, X2_PREV_new, SRFTEMP, ST_new


def test_stemp_epic_askee_like() -> Tuple[
    float, List[float], float, List[float], int, List[int], float, float, List[float]
]:
    """
    Test derived from ASKEE_DSSAT_EPIC.for program.
    Does not execute at import; returns initialization outputs.

    Returns (from initialize):
    - CUMDPT, DSMID, TDL, TMA, NDays, WetDay, X2_PREV, SRFTEMP, ST
    """
    # Soil setup
    SOILPROP_NLAYR = 4
    SOILPROP_BD = [1.6] * SOILPROP_NLAYR
    SOILPROP_DLAYR = [10.0] * SOILPROP_NLAYR
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * SOILPROP_NLAYR
    SOILPROP_LL = [0.2] * SOILPROP_NLAYR

    # Weather and states
    ISWITCH_ISWWAT = 'Y'
    SW = [0.2] * SOILPROP_NLAYR
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0
    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0
    X2_PREV = 0.0

    return stemp_epic_initialize(
        ISWITCH_ISWWAT=ISWITCH_ISWWAT,
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        SOILPROP_NLAYR=SOILPROP_NLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW,
        X2_PREV=X2_PREV,
    )