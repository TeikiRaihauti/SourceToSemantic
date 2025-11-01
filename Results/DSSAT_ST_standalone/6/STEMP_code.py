from typing import Sequence, List, Tuple
import math


def fortran_nint(x: float) -> int:
    """
    Fortran NINT equivalent: nearest integer with halves rounded away from zero.
    """
    if x >= 0.0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


def YR_DOY(YRDOY: int) -> Tuple[int, int]:
    """
    Convert YRDOY to YEAR and DOY.

    Parameters
    - YRDOY: int

    Returns
    - (YEAR: int, DOY: int)
    """
    YR = int(YRDOY // 1000)
    DOY = YRDOY - YR * 1000
    return YR, DOY


def SOILT(
    ALBEDO: float,  # not used in algorithm but preserved
    B: float,
    CUMDPT: float,
    DOY: int,
    DP: float,
    HDAY: float,
    NLAYR: int,
    PESW: float,
    SRAD: float,  # not used in algorithm but preserved
    TAMP: float,
    TAV: float,
    TAVG: float,
    TMAX: float,  # not used in algorithm but preserved
    WW: float,
    DSMID: Sequence[float],
    ATOT: float,
    TMA: Sequence[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    SOILT - Determines soil temperature by layer.

    Inputs
    - ALBEDO: float
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float
    - WW: float
    - DSMID: Sequence[float] length >= NLAYR
    - ATOT: float
    - TMA: Sequence[float] length 5

    Returns (updated)
    - ATOT: float
    - TMA: List[float] length 5
    - SRFTEMP: float
    - ST: List[float] length NLAYR
    """
    # Copy to mutable structures
    TMA_list = [float(v) for v in TMA]
    DSMID_list = [float(v) for v in DSMID[:NLAYR]]

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT = ATOT - TMA_list[4]

    # Shift TMA history
    for K in range(4, 0, -1):
        TMA_list[K] = TMA_list[K - 1]
    TMA_list[0] = TAVG

    # Keep only 4 decimals (debug vs release fix)
    TMA_list[0] = fortran_nint(TMA_list[0] * 10000.0) / 10000.0
    ATOT = ATOT + TMA_list[0]

    # Corrected water content function (EPIC)
    # WC (ratio), PESW (cm), WW (dimensionless), CUMDPT (mm)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    ST: List[float] = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID_list[L] / DD
        ST_val = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[L] = fortran_nint(ST_val * 1000.0) / 1000.0  # 3 decimals

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA_list, SRFTEMP, ST


def STEMP_initialize(
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: Sequence[float],
    SOILPROP_DS: Sequence[float],
    SOILPROP_DUL: Sequence[float],
    SOILPROP_LL: Sequence[float],
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
    SRAD: float,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CONTROL_YRDOY: int,
) -> Tuple[float, List[float], float, List[float], float, float, List[float], float]:
    """
    Seasonal initialization for soil temperature (STEMP).

    Inputs
    - ISWITCH_ISWWAT: str ('Y' or 'N')
    - SOILPROP_BD: Sequence[float] bulk density by layer
    - SOILPROP_DS: Sequence[float] cumulative depth to bottom of layer (cm)
    - SOILPROP_DUL: Sequence[float] drained upper limit by layer (vol frac)
    - SOILPROP_LL: Sequence[float] lower limit by layer (vol frac)
    - SOILPROP_NLAYR: int number of layers
    - SOILPROP_MSALB: float mulch/soil albedo
    - SRAD: float solar radiation (MJ/m2/d) [not used here but passed through]
    - SW: Sequence[float] soil water by layer (vol frac)
    - TAVG: float average temperature (degC)
    - TMAX: float maximum temperature (degC) [not used by algorithm]
    - XLAT: float latitude (deg)
    - TAV: float annual average air temperature (degC)
    - TAMP: float annual temperature amplitude (degC)
    - CONTROL_YRDOY: int YYYYDDD

    Returns (initial state)
    - CUMDPT: float profile depth (mm)
    - DSMID: List[float] depth to mid-point of each layer (mm)
    - TDL: float sum of DUL*thickness (cm at init)
    - TMA: List[float] moving average temperatures (5-day history)
    - ATOT: float sum of TMA over 5 days
    - SRFTEMP: float surface temperature (degC)
    - ST: List[float] soil temperature by layer (degC)
    - HDAY: float hottest day of year (DOY) (20 for south, 200 for north)
    """
    NLAYR = int(SOILPROP_NLAYR)
    BD = [float(v) for v in SOILPROP_BD[:NLAYR]]
    DSI = [float(v) for v in SOILPROP_DS[:NLAYR]]
    DUL = [float(v) for v in SOILPROP_DUL[:NLAYR]]
    LL = [float(v) for v in SOILPROP_LL[:NLAYR]]
    SWI = [float(v) for v in SW[:NLAYR]]
    MSALB = float(SOILPROP_MSALB)

    # Compute DOY (used within SOILT)
    _, DOY = YR_DOY(int(CONTROL_YRDOY))

    # HDAY by hemisphere
    HDAY = 20.0 if XLAT < 0.0 else 200.0

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0
    DSMID: List[float] = [0.0] * NLAYR

    for L in range(NLAYR):
        if L == 0:
            DLI = DSI[L]
        else:
            DLI = DSI[L] - DSI[L - 1]
        DSMID[L] = CUMDPT + DLI * 5.0  # mm to middle of layer
        CUMDPT = CUMDPT + DLI * 10.0   # mm profile depth
        TBD += BD[L] * DLI
        TLL += LL[L] * DLI
        TSW += SWI[L] * DLI
        TDL += DUL[L] * DLI

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # not used in SOILT

    # Initialize moving average and layer temps
    TMA: List[float] = [fortran_nint(TAVG * 10000.0) / 10000.0 for _ in range(5)]
    ATOT = TMA[0] * 5.0
    ST: List[float] = [float(TAVG) for _ in range(NLAYR)]
    SRFTEMP = float(TAVG)

    # Spin-up
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
    ISWITCH_ISWWAT: str,
    SOILPROP_BD: Sequence[float],
    SOILPROP_DLAYR: Sequence[float],
    SOILPROP_DS: Sequence[float],
    SOILPROP_DUL: Sequence[float],
    SOILPROP_LL: Sequence[float],
    SOILPROP_NLAYR: int,
    SOILPROP_MSALB: float,
    SRAD: float,
    SW: Sequence[float],
    TAVG: float,
    TMAX: float,
    TAV: float,
    TAMP: float,
    CONTROL_YRDOY: int,
    HDAY: float,
    CUMDPT: float,
    DSMID: Sequence[float],
    TDL: float,
    TMA: Sequence[float],
    ATOT: float,
    SRFTEMP: float,
    ST: Sequence[float],
) -> Tuple[float, List[float], float, List[float]]:
    """
    STEMP - Daily soil temperature update (main biophysical process).

    Inputs
    - ISWITCH_ISWWAT: str ('Y' or 'N')
    - SOILPROP_BD: Sequence[float] bulk density by layer
    - SOILPROP_DLAYR: Sequence[float] layer thickness (mm) by layer
    - SOILPROP_DS: Sequence[float] cumulative depth to bottom of layer (cm)
    - SOILPROP_DUL: Sequence[float] drained upper limit by layer (vol frac)
    - SOILPROP_LL: Sequence[float] lower limit by layer (vol frac)
    - SOILPROP_NLAYR: int number of layers
    - SOILPROP_MSALB: float mulch/soil albedo
    - SRAD: float solar radiation (MJ/m2/d) [not used by algorithm]
    - SW: Sequence[float] soil water by layer (vol frac)
    - TAVG: float average temperature (degC)
    - TMAX: float maximum temperature (degC) [not used by algorithm]
    - TAV: float annual average air temperature (degC)
    - TAMP: float annual temperature amplitude (degC)
    - CONTROL_YRDOY: int YYYYDDD
    - HDAY: float hottest day of year (DOY), from initialization
    - CUMDPT: float profile depth (mm), from initialization
    - DSMID: Sequence[float] depth to mid-point of each layer (mm), from initialization
    - TDL: float state variable (accumulated; will be incremented)
    - TMA: Sequence[float] moving average temperatures (5-day history)
    - ATOT: float sum of TMA over 5 days
    - SRFTEMP: float surface temperature (degC)
    - ST: Sequence[float] soil temperature by layer (degC)

    Returns (updated)
    - TDL: float
    - TMA: List[float]
    - ATOT: float
    - SRFTEMP: float
    - ST: List[float]
    """
    NLAYR = int(SOILPROP_NLAYR)
    BD = [float(v) for v in SOILPROP_BD[:NLAYR]]
    DLAYR = [float(v) for v in SOILPROP_DLAYR[:NLAYR]]
    DS = [float(v) for v in SOILPROP_DS[:NLAYR]]
    DUL = [float(v) for v in SOILPROP_DUL[:NLAYR]]
    LL = [float(v) for v in SOILPROP_LL[:NLAYR]]
    SW_list = [float(v) for v in SW[:NLAYR]]
    MSALB = float(SOILPROP_MSALB)
    DSMID_list = [float(v) for v in DSMID[:NLAYR]]
    TMA_list = [float(v) for v in TMA]
    ST_list = [float(v) for v in ST[:NLAYR]]

    _, DOY = YR_DOY(int(CONTROL_YRDOY))

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    # Note: TDL is not reset daily (preserved original behavior)
    for L in range(NLAYR):
        TBD += BD[L] * DLAYR[L]
        TDL += DUL[L] * DLAYR[L]
        TLL += LL[L] * DLAYR[L]
        TSW += SW_list[L] * DLAYR[L]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB  # not used in SOILT

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)
    else:
        PESW = max(0.0, TDL - TLL)

    ATOT, TMA_list, SRFTEMP, ST_list = SOILT(
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
        DSMID=DSMID_list,
        ATOT=ATOT,
        TMA=TMA_list,
    )

    return TDL, TMA_list, ATOT, SRFTEMP, ST_list