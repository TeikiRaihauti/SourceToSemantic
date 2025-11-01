from typing import List, Sequence, Tuple
import math


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
    SRFTEMP: float,
    ST: Sequence[float],
    X2_AVG: float,
    X2_PREV: float
) -> Tuple[List[float], float, List[float], float, float]:
    """
    Soil temperature by layer (EPIC method).

    Inputs:
    - B: float
    - BCV: float
    - CUMDPT: float
    - DP: float
    - DSMID: Sequence[float], length >= NLAYR
    - NLAYR: int
    - PESW: float
    - TAV: float
    - TAVG: float
    - TMAX: float
    - TMIN: float
    - WetDay: int (0 or 1)
    - WFT: float
    - WW: float
    - TMA: Sequence[float], length 5 (in/out)
    - SRFTEMP: float (in/out)
    - ST: Sequence[float], length >= NLAYR (in/out)
    - X2_AVG: float (output in Fortran; here computed and returned)
    - X2_PREV: float (in/out)

    Returns:
    - TMA: List[float], updated length 5
    - SRFTEMP: float
    - ST: List[float], updated length NLAYR
    - X2_AVG: float
    - X2_PREV: float, updated
    """
    # Ensure lists for mutation
    TMA_list = list(TMA[:5])
    ST_list = list(ST[:NLAYR])

    # Water content term (WC): ratio
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX_local = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX_local * DP

    # Surface driver X2 based on wet/dry day
    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # EPIC memory (note: preserves original Fortran order)
    TMA_list[0] = X2
    for K in range(4, 0, -1):  # 5..2 in Fortran -> 4..1 in Python 0-based
        TMA_list[K] = TMA_list[K - 1]

    X2_AVG_out = sum(TMA_list) / 5.0

    # Soil cover influence
    X3 = (1.0 - BCV) * X2_AVG_out + BCV * X2_PREV
    SRFTEMP_out = min(X2_AVG_out, X3)

    # Layer temperatures
    X1 = TAV - X3
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID[L] / DD
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST_list[L] = LAG * ST_list[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_out = X2_AVG_out

    return TMA_list, SRFTEMP_out, ST_list, X2_AVG_out, X2_PREV_out


def STEMP_EPIC(
    DYNAMIC: int,
    ISWWAT: str,
    BD: Sequence[float],
    DLAYR: Sequence[float],
    DS: Sequence[float],
    DUL: Sequence[float],
    LL: Sequence[float],
    NLAYR: int,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    TMIN: float,
    TAV: float,
    WEATHER_RAIN: float,
    ORGC_MULCHMASS: float,
    WATER_SNOW: float,
    MGMT_DEPIR: float,
    PLANT_BIOMAS: float,
    CUMDPT: float,
    DSMID: Sequence[float],
    TDL: float,
    TMA: Sequence[float],
    NDays: int,
    WetDay: Sequence[int],
    X2_PREV: float,
    SRFTEMP: float,
    ST: Sequence[float]
) -> Tuple[float, List[float], float, List[float], int, List[int], float, float, List[float]]:
    """
    Soil temperature driver (EPIC method). Handles seasonal initialization and daily rate.

    Inputs:
    - DYNAMIC: int (2=SEASINIT, 3=RATE, 5=OUTPUT, 6=SEASEND)
    - ISWWAT: str ('Y' if water simulated)
    - BD: Sequence[float], bulk density by layer
    - DLAYR: Sequence[float], layer thickness (cm)
    - DS: Sequence[float], cumulative depth (cm) by layer
    - DUL: Sequence[float], drained upper limit by layer
    - LL: Sequence[float], lower limit by layer
    - NLAYR: int, number of layers
    - SW: Sequence[float], volumetric water content by layer
    - TAVG: float, daily average air temperature
    - TMAX: float, daily maximum air temperature
    - TMIN: float, daily minimum air temperature
    - TAV: float, annual average temperature
    - WEATHER_RAIN: float, daily rainfall (mm)
    - ORGC_MULCHMASS: float, mulch mass (kg/ha)
    - WATER_SNOW: float, snow water equivalent (mm)
    - MGMT_DEPIR: float, daily irrigation depth (mm)
    - PLANT_BIOMAS: float, plant biomass (kg/ha)
    - CUMDPT: float (in/out), cumulative profile depth (mm)
    - DSMID: Sequence[float] (in/out), depth to midpoint of layer (mm)
    - TDL: float (in/out), sum DUL*DLAYR (cm*cm)
    - TMA: Sequence[float] (in/out), 5-element memory array
    - NDays: int (in/out), number of days tracked in wet-day memory (<=30)
    - WetDay: Sequence[int] (in/out), wet-day memory (0/1), last NDays entries
    - X2_PREV: float (in/out), previous X2 average for cover effect
    - SRFTEMP: float (in/out), surface temperature
    - ST: Sequence[float] (in/out), soil temperature by layer

    Returns (in this exact order):
    - CUMDPT: float
    - DSMID: List[float]
    - TDL: float
    - TMA: List[float]
    - NDays: int
    - WetDay: List[int]
    - X2_PREV: float
    - SRFTEMP: float
    - ST: List[float]
    """
    # Ensure lists for mutation
    BD_list = list(BD[:NLAYR])
    DLAYR_list = list(DLAYR[:NLAYR])
    DS_list = list(DS[:NLAYR])
    DUL_list = list(DUL[:NLAYR])
    LL_list = list(LL[:NLAYR])
    SW_list = list(SW[:NLAYR])
    DSMID_list = list(DSMID[:NLAYR]) if len(DSMID) >= NLAYR else [0.0] * NLAYR
    TMA_list = list(TMA[:5]) if len(TMA) >= 5 else [0.0] * 5
    ST_list = list(ST[:NLAYR]) if len(ST) >= NLAYR else [TAVG] * NLAYR
    WetDay_list = list(WetDay[:NDays]) if NDays > 0 else []

    if DYNAMIC == 2:  # SEASINIT
        SWI = SW_list[:]

        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_out = 0.0
        CUMDPT_out = 0.0
        for L in range(NLAYR):
            DSMID_list[L] = CUMDPT_out + DLAYR_list[L] * 5.0
            CUMDPT_out += DLAYR_list[L] * 10.0
            TBD += BD_list[L] * DLAYR_list[L]
            TLL += LL_list[L] * DLAYR_list[L]
            TSW += SWI[L] * DLAYR_list[L]
            TDL_out += DUL_list[L] * DLAYR_list[L]

        if ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_out - TLL)

        ABD = TBD / DS_list[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)

        TMA_list = [round(TAVG, 4)] * 5
        X2_AVG = TMA_list[0] * 5.0

        ST_list = [TAVG] * NLAYR

        WFT = 0.1
        WetDay_list = []
        NDays_out = 0

        CV = ORGC_MULCHMASS / 1000.0
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV + math.exp(5.3396 - 2.3951 * CV) != 0 else 0.0
        BCV2 = WATER_SNOW / (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) if WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW) != 0 else 0.0
        BCV = max(BCV1, BCV2)

        SRFTEMP_out = SRFTEMP
        X2_PREV_out = X2_PREV
        for _ in range(8):
            TMA_list, SRFTEMP_out, ST_list, X2_AVG, X2_PREV_out = SOILT_EPIC(
                B, BCV, CUMDPT_out, DP, DSMID_list, NLAYR, PESW, TAV,
                TAVG, TMAX, TMIN, 0, WFT, WW, TMA_list, SRFTEMP_out, ST_list, X2_AVG, X2_PREV_out
            )

        return (
            CUMDPT_out,
            DSMID_list,
            TDL_out,
            TMA_list,
            NDays_out,
            WetDay_list,
            X2_PREV_out,
            SRFTEMP_out,
            ST_list
        )

    elif DYNAMIC == 3:  # RATE
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_out = TDL
        for L in range(NLAYR):
            TBD += BD_list[L] * DLAYR_list[L]
            TDL_out += DUL_list[L] * DLAYR_list[L]
            TLL += LL_list[L] * DLAYR_list[L]
            TSW += SW_list[L] * DLAYR_list[L]

        ABD = TBD / DS_list[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)

        if ISWWAT == 'Y':
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_out - TLL)

        RAIN = WEATHER_RAIN
        if NDays == 30:
            WetDay_list = WetDay_list[1:] if len(WetDay_list) > 0 else []
            NDays_out = 30
        else:
            NDays_out = NDays + 1

        wet_flag = 1 if (RAIN + MGMT_DEPIR) > 1.0e-6 else 0
        WetDay_list.append(wet_flag)
        if NDays_out > 30:
            NDays_out = 30
            WetDay_list = WetDay_list[-30:]

        NWetDays = sum(WetDay_list)
        WFT = float(NWetDays) / float(NDays_out) if NDays_out > 0 else 0.0

        CV = (PLANT_BIOMAS + ORGC_MULCHMASS) / 1000.0
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if CV + math.exp(5.3396 - 2.3951 * CV) != 0 else 0.0
        BCV2 = WATER_SNOW / (WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW)) if WATER_SNOW + math.exp(2.303 - 0.2197 * WATER_SNOW) != 0 else 0.0
        BCV = max(BCV1, BCV2)

        last_wet = WetDay_list[-1] if NDays_out > 0 else 0
        X2_AVG = 0.0  # will be computed inside SOILT_EPIC
        TMA_list, SRFTEMP_out, ST_list, X2_AVG, X2_PREV_out = SOILT_EPIC(
            B, BCV, CUMDPT, DP, DSMID_list, NLAYR, PESW, TAV,
            TAVG, TMAX, TMIN, last_wet, WFT, WW, TMA_list, SRFTEMP, ST_list, X2_AVG, X2_PREV
        )

        return (
            CUMDPT,
            DSMID_list,
            TDL_out,
            TMA_list,
            NDays_out,
            WetDay_list,
            X2_PREV_out,
            SRFTEMP_out,
            ST_list
        )

    else:
        # OUTPUT or SEASEND or other states: no changes to biophysical state
        return (
            CUMDPT,
            DSMID_list,
            TDL,
            TMA_list,
            NDays,
            WetDay_list,
            X2_PREV,
            SRFTEMP,
            ST_list
        )