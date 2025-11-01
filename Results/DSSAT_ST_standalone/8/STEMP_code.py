from typing import List, Tuple
import math


def _fortran_nint(x: float) -> int:
    # Fortran NINT: nearest integer, halves away from zero
    return int(math.floor(x + 0.5)) if x >= 0.0 else int(math.ceil(x - 0.5))


def _round_fortran(x: float, decimals: int) -> float:
    # Fortran-style rounding to given decimals using NINT
    scale = 10 ** decimals
    return _fortran_nint(x * scale) / scale


def SOILT(
    ALBEDO: float,  # dimensionless
    B: float,  # dimensionless
    CUMDPT: float,  # mm profile depth
    DOY: int,  # day of year (1-366)
    DP: float,  # mm damping depth
    HDAY: float,  # hottest day of year (20 SH, 200 NH)
    NLAYR: int,  # number of soil layers used
    PESW: float,  # cm plant-extractable soil water in profile
    SRAD: float,  # MJ/m2/d (not used here)
    TAMP: float,  # C, annual temperature amplitude
    TAV: float,  # C, annual average temperature
    TAVG: float,  # C, current day's average air temperature
    TMAX: float,  # C, current day's maximum air temperature (not used)
    WW: float,  # dimensionless, fraction
    DSMID: List[float],  # mm, depth to midpoint of layers (size NLAYR)
    ATOT: float,  # running sum of last 5 TAVG values
    TMA: List[float],  # last 5 TAVG values (length 5)
) -> Tuple[float, List[float], float, List[float]]:
    """
    Compute soil temperature by layer.

    Inputs:
      ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR, PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, ATOT, TMA

    Returns:
      ATOT_out: float
      TMA_out: List[float] length 5
      SRFTEMP: float
      ST: List[float] length NLAYR
    """
    # Copy and prepare local state
    TMA_work = [float(v) for v in TMA[:5]]
    ATOT_work = float(ATOT)

    # Phase angle for annual cycle
    ALX = (float(DOY) - HDAY) * 0.0174

    # Update running 5-day average storage
    ATOT_work = ATOT_work - TMA_work[4]
    for K in range(4, 0, -1):
        TMA_work[K] = TMA_work[K - 1]
    TMA_work[0] = _round_fortran(TAVG, 4)
    ATOT_work = ATOT_work + TMA_work[0]

    # Corrected water content function (ratio)
    # WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0
    # PESW in cm, CUMDPT in mm, multiply by 10 mm/cm
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # damping depth (mm)

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_work / 5.0 - TA

    ST = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD
        ST_L = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[L] = _round_fortran(ST_L, 3)

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_work, TMA_work, SRFTEMP, ST


def STEMP_seasonal_init(
    CONTROL_YRDOY: int,  # YYYY*1000 + DOY
    ISWITCH_ISWWAT: str,  # 'Y' or 'N'
    SOILPROP_BD: List[float],  # g/cm3, length >= NLAYR
    SOILPROP_DLAYR: List[float],  # cm, not used in init thickness derivation
    SOILPROP_DS: List[float],  # cm, cumulative depth at bottom of layer, length >= NLAYR
    SOILPROP_DUL: List[float],  # cm3/cm3, length >= NLAYR
    SOILPROP_LL: List[float],  # cm3/cm3, length >= NLAYR
    SOILPROP_NLAYR: int,  # number of layers
    SOILPROP_MSALB: float,  # albedo
    SRAD: float,  # MJ/m2/d
    SW: List[float],  # cm3/cm3, initial soil water by layer (length >= NLAYR)
    TAVG: float,  # C, today's average air temp
    TMAX: float,  # C, today's max air temp
    XLAT: float,  # degrees, latitude
    TAV: float,  # C, annual average air temperature
    TAMP: float,  # C, annual temperature amplitude
) -> Tuple[
    float,  # CUMDPT (mm)
    List[float],  # DSMID (mm), length NLAYR
    float,  # TDL (cm water)
    List[float],  # TMA (len 5)
    float,  # ATOT
    float,  # SRFTEMP (C)
    List[float],  # ST (C), length NLAYR
    float,  # HDAY
]:
    """
    Seasonal initialization for soil temperature.

    Returns:
      CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY
    """
    # Determine DOY
    YEAR = CONTROL_YRDOY // 1000
    DOY = CONTROL_YRDOY - YEAR * 1000

    NLAYR = SOILPROP_NLAYR
    BD = SOILPROP_BD
    DUL = SOILPROP_DUL
    LL = SOILPROP_LL
    DS = SOILPROP_DS
    MSALB = SOILPROP_MSALB

    # Southern/Northern hemisphere hottest day
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    # Initialize aggregates
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    # Layer thicknesses from cumulative depth (cm)
    DLI = [0.0] * NLAYR
    DSMID = [0.0] * NLAYR  # mm to midpoint of layer
    for L in range(NLAYR):
        if L == 0:
            DLI[L] = DS[L]
        else:
            DLI[L] = DS[L] - DS[L - 1]
        DSMID[L] = CUMDPT + DLI[L] * 5.0  # mm
        CUMDPT = CUMDPT + DLI[L] * 10.0  # mm
        TBD += BD[L] * DLI[L]
        TLL += LL[L] * DLI[L]
        TSW += SW[L] * DLI[L]
        TDL += DUL[L] * DLI[L]

    if ISWITCH_ISWWAT.upper() == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # not used directly in SOILT

    # Initialize running temp memory
    TMA = [_round_fortran(TAVG, 4) for _ in range(5)]
    ATOT = TMA[0] * 5.0

    # Initialize soil temperatures
    ST = [TAVG for _ in range(NLAYR)]
    SRFTEMP = TAVG

    # Spin-up SOILT 8 times
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(
            ALBEDO=ALBEDO,
            B=B,
            CUMDPT=CUMDPT,
            DOY=DOY,
            DP=DP,
            HDAY=HDAY,
            NLAYR=NLAYR,
            PESW=PESW,
            SRAD=SRAD,
            TAMP=TAMP,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            WW=WW,
            DSMID=DSMID,
            ATOT=ATOT,
            TMA=TMA,
        )

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def STEMP(
    ISWITCH_ISWWAT: str,  # 'Y' or 'N'
    SRAD: float,  # MJ/m2/d
    SW: List[float],  # cm3/cm3 current soil water by layer
    TAVG: float,  # C, daily average air temp
    TMAX: float,  # C, max air temp
    TAV: float,  # C, annual average air temp
    TAMP: float,  # C, annual temperature amplitude
    SOILPROP_BD: List[float],  # g/cm3
    SOILPROP_DLAYR: List[float],  # cm, layer thickness
    SOILPROP_DS: List[float],  # cm, cumulative depth at bottom of layer
    SOILPROP_DUL: List[float],  # cm3/cm3
    SOILPROP_LL: List[float],  # cm3/cm3
    SOILPROP_NLAYR: int,  # layers
    SOILPROP_MSALB: float,  # albedo
    DOY: int,  # day of year
    CUMDPT: float,  # mm, from init
    DSMID: List[float],  # mm, from init
    TDL: float,  # cm water, running (note: accumulates as in original code)
    TMA: List[float],  # last 5 TAVG values
    ATOT: float,  # running sum of last 5 TAVG values
    SRFTEMP: float,  # C, previous day value (not used in calc)
    ST: List[float],  # C, previous day soil temps (not used in calc)
    HDAY: float,  # hottest day (from init)
) -> Tuple[
    float,  # ATOT
    List[float],  # TMA
    float,  # SRFTEMP
    List[float],  # ST
    float,  # TDL
]:
    """
    Daily soil temperature rate calculations.

    Inputs include environment (SRAD, TAVG, TMAX, TAV, TAMP, DOY),
      soil properties (BD, DLAYR, DS, DUL, LL, NLAYR, MSALB),
      daily soil water (SW),
      and state variables from initialization (CUMDPT, DSMID, TDL, TMA, ATOT, HDAY).

    Returns:
      ATOT, TMA, SRFTEMP, ST, TDL (updated)
    """
    NLAYR = SOILPROP_NLAYR
    BD = SOILPROP_BD
    DLAYR = SOILPROP_DLAYR
    DUL = SOILPROP_DUL
    LL = SOILPROP_LL
    DS = SOILPROP_DS
    MSALB = SOILPROP_MSALB

    # Aggregate properties across profile
    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: Original code does not reinitialize TDL here; it accumulates
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX  # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT.upper() == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ATOT, TMA, SRFTEMP, ST = SOILT(
        ALBEDO=ALBEDO,
        B=B,
        CUMDPT=CUMDPT,
        DOY=DOY,
        DP=DP,
        HDAY=HDAY,
        NLAYR=NLAYR,
        PESW=PESW,
        SRAD=SRAD,
        TAMP=TAMP,
        TAV=TAV,
        TAVG=TAVG,
        TMAX=TMAX,
        WW=WW,
        DSMID=DSMID,
        ATOT=ATOT,
        TMA=TMA,
    )

    return ATOT, TMA, SRFTEMP, ST, TDL