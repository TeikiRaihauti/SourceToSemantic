from math import ceil, floor, exp, log
from typing import List, Tuple, Optional


def init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    cAVT: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize soil temperature profile and layer depths.

    Inputs:
    - cSoilLayerDepth: List[float] Depth to bottom of each soil layer [m]
    - cFirstDayMeanTemp: float Mean air temperature on first day [°C]
    - cDampingDepth: float Initial damping depth of soil [m]
    - cAVT: float Long-term average annual air temperature [°C]

    Returns:
    - SoilTempArray: List[float] Initial soil temperature for each (possibly extended) layer [°C]
    - rSoilTempArrayRate: List[float] Initial daily soil temperature change rates (zeros) [°C/day]
    - pSoilLayerDepth: List[float] Depth to bottom of each (possibly extended) layer [m]
    """
    Z = cSoilLayerDepth
    tProfileDepth = Z[-1]
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(abs(ceil(additionalDepth))) + len(Z)
    tStmp = [0.0] * layers
    tStmpRate = [0.0] * layers
    tz = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = depth
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def process(
    cSoilLayerDepth: List[float],
    cABD: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    cAVT: float,
    SoilTempArray: List[float],
    pSoilLayerDepth: List[float],
    iDoInitialize: bool = False,
    cFirstDayMeanTemp: Optional[float] = None,
    cDampingDepth: Optional[float] = None,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Compute daily soil temperature updates for each layer.

    Inputs:
    - cSoilLayerDepth: List[float] Depth to bottom of each soil layer [m]
    - cABD: float Average bulk density [t/m3]
    - iSoilWaterContent: float Water content, sum over soil profile [mm]
    - iSoilSurfaceTemperature: float Soil surface temperature [°C]
    - cAVT: float Long-term average annual air temperature [°C]
    - SoilTempArray: List[float] Current soil temperature in layers [°C]
    - pSoilLayerDepth: List[float] Depth to bottom of each (possibly extended) layer [m]
    - iDoInitialize: bool If True, re-initialize state before processing
    - cFirstDayMeanTemp: Optional[float] Required if iDoInitialize is True [°C]
    - cDampingDepth: Optional[float] Required if iDoInitialize is True [m]

    Returns:
    - SoilTempArray: List[float] Updated soil temperatures [°C]
    - rSoilTempArrayRate: List[float] Daily temperature change for each layer [°C/day]
    - pSoilLayerDepth: List[float] (Possibly re-initialized) layer depths [m]
    """
    if iDoInitialize:
        if cFirstDayMeanTemp is None or cDampingDepth is None:
            raise ValueError("cFirstDayMeanTemp and cDampingDepth must be provided when iDoInitialize is True.")
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )
    else:
        rSoilTempArrayRate = [0.0] * len(SoilTempArray)

    tZp = pSoilLayerDepth
    tZc = cSoilLayerDepth

    XLAG = 0.8
    XLG1 = 1.0 - XLAG
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    newSoilTempArray = list(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        base = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        old_temp = SoilTempArray[i]
        RATE = XLG1 * (base + iSoilSurfaceTemperature - old_temp)
        Z1 = tZp[i]
        rSoilTempArrayRate[i] = RATE
        newSoilTempArray[i] = old_temp + RATE

    return newSoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def test_STMPsimCalculator_case0() -> dict:
    """
    Test derived from STMPsimCalculator.fillTestVariables (aTestIndex == 0).
    Returns a dictionary with computed and expected outputs for comparison.
    """
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        cSoilLayerDepth=cSoilLayerDepth,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        cAVT=cAVT,
        SoilTempArray=SoilTempArray,
        pSoilLayerDepth=pSoilLayerDepth,
    )

    expected_SoilTempArray = [
        13.624360856350041,
        13.399968504634286,
        12.599999999999845,
        12.2,
        11.4,
        10.6,
        9.799999999999999,
        9.0,
    ]

    return {
        "SoilTempArray": SoilTempArray,
        "expected_SoilTempArray": expected_SoilTempArray,
        "rSoilTempArrayRate": rSoilTempArrayRate,
        "pSoilLayerDepth": pSoilLayerDepth,
    }