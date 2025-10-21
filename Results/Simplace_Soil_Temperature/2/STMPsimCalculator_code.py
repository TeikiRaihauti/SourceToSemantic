from math import floor, ceil, fabs, exp, log
from typing import List, Tuple


def init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    cAVT: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize soil temperature profile and extended layer depths.

    Inputs:
    - cSoilLayerDepth: List[float], depth to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: float, mean air temperature on first day (°C)
    - cDampingDepth: float, initial damping depth of soil (m)
    - cAVT: float, long-term average annual air temperature (°C)

    Returns:
    - SoilTempArray: List[float], initialized soil temperature profile (°C)
    - rSoilTempArrayRate: List[float], initialized daily temperature change rates (°C d-1)
    - pSoilLayerDepth: List[float], extended depth to bottom of each layer including added layers to reach cDampingDepth (m)
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp

    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(fabs(ceil(additionalDepth))) + len(Z)

    tStmp = [0.0] * layers
    tStmpRate = [0.0] * layers
    tz = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            # first additional layer might be smaller than 1 m
            depth = tProfileDepth + firstAdditionalLayerHight + (i - len(Z))
        tz[i] = depth
        # Linear approximation between surface (first day mean) and damping depth (AVT)
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz
    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def process(
    iDoInitialize: bool,
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cAVT: float,
    cABD: float,
    cDampingDepth: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Advance soil temperature profile by one time step.

    Inputs:
    - iDoInitialize: bool, if True reinitialize state using init() before processing
    - SoilTempArray: List[float], current soil temperature profile (°C)
    - rSoilTempArrayRate: List[float], current daily temperature change rates (°C d-1)
    - pSoilLayerDepth: List[float], extended depth to bottom of each layer (m)
    - cSoilLayerDepth: List[float], original profile depths to bottom of each soil layer (m)
    - cFirstDayMeanTemp: float, mean air temperature on first day (°C) used if reinitializing
    - cAVT: float, long-term average annual air temperature (°C)
    - cABD: float, average bulk density of the soil profile (t m-3)
    - cDampingDepth: float, initial damping depth of soil (m) used if reinitializing
    - iSoilWaterContent: float, water stored in the soil profile (mm)
    - iSoilSurfaceTemperature: float, soil surface temperature (°C)

    Returns:
    - SoilTempArray: List[float], updated soil temperature profile (°C)
    - rSoilTempArrayRate: List[float], updated daily temperature change rates (°C d-1)
    - pSoilLayerDepth: List[float], (possibly reinitialized) extended layer depths (m)
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)

    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    n = len(SoilTempArray)
    # Ensure rate array length matches temperature array length
    if len(rSoilTempArrayRate) != n:
        rSoilTempArrayRate = [0.0] * n

    for i in range(n):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = tZp[i]
        rSoilTempArrayRate[i] = RATE
        SoilTempArray[i] = SoilTempArray[i] + rSoilTempArrayRate[i]

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def test_STMPsimCalculator_case0() -> Tuple[List[float], List[float], bool]:
    """
    Test derived from STMPsimCalculator.fillTestVariables test 0.

    Returns:
    - result: List[float], computed SoilTempArray after one process step
    - expected: List[float], expected SoilTempArray
    - passed: bool, True if within tolerance
    """
    # Define inputs (before process)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    # Process one step
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        iDoInitialize=False,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cABD=cABD,
        cDampingDepth=cDampingDepth,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
    )

    expected = [
        13.624360856350041,
        13.399968504634286,
        12.599999999999845,
        12.2,
        11.4,
        10.6,
        9.799999999999999,
        9.0,
    ]

    tol = 1e-9
    passed = len(SoilTempArray) == len(expected) and all(
        abs(a - b) <= tol for a, b in zip(SoilTempArray, expected)
    )
    return SoilTempArray, expected, passed