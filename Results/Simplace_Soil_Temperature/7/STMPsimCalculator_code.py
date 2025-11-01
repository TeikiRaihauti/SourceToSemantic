from typing import List, Tuple
from math import floor, ceil, exp, log


def init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    cAVT: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize soil temperature profile and extended layer depths.

    Inputs:
    - cSoilLayerDepth: List[float] depths to bottom of each soil layer (m), ascending
    - cFirstDayMeanTemp: float mean air temperature on first day (°C)
    - cDampingDepth: float initial damping depth of soil (m)
    - cAVT: float long-term average annual air temperature (°C)

    Returns:
    - SoilTempArray: List[float] initial soil temperatures for each (extended) layer (°C)
    - rSoilTempArrayRate: List[float] initial rate of daily temperature change (°C d-1), zeros
    - pSoilLayerDepth: List[float] depths to bottom of each (extended) layer (m)
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(abs(ceil(additionalDepth))) + len(Z)

    tStmp: List[float] = [0.0] * layers
    tStmpRate: List[float] = [0.0] * layers
    tz: List[float] = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            depth = tProfileDepth + firstAdditionalLayerHight + (i - len(Z))
        tz[i] = depth
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def process(
    cSoilLayerDepth: List[float],
    cAVT: float,
    cABD: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    iDoInitialize: bool,
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
) -> Tuple[List[float], List[float], List[float]]:
    """
    Update soil temperatures for one time step.

    Inputs:
    - cSoilLayerDepth: List[float] depths to bottom of each soil layer (m), ascending
    - cAVT: float long-term average annual air temperature (°C)
    - cABD: float average soil bulk density (t m-3)
    - iSoilWaterContent: float water stored in soil profile (mm)
    - iSoilSurfaceTemperature: float soil surface temperature DST0 (°C)
    - iDoInitialize: bool if True, re-initialize state using init inputs below
    - cFirstDayMeanTemp: float mean air temperature on first day (°C), used if iDoInitialize
    - cDampingDepth: float initial damping depth of soil (m), used if iDoInitialize
    - SoilTempArray: List[float] current soil temperatures (°C)
    - rSoilTempArrayRate: List[float] current rate of daily temperature change (°C d-1)
    - pSoilLayerDepth: List[float] depths to bottom of each (extended) layer (m)

    Returns:
    - SoilTempArray: List[float] updated soil temperatures (°C)
    - rSoilTempArrayRate: List[float] updated rate of daily temperature change (°C d-1)
    - pSoilLayerDepth: List[float] (unchanged) depths to bottom of each (extended) layer (m)
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )

    tZp = pSoilLayerDepth
    tZc = cSoilLayerDepth
    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    new_r: List[float] = [0.0] * len(SoilTempArray)
    new_T: List[float] = list(SoilTempArray)

    for i in range(len(new_T)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - new_T[i])
        new_r[i] = RATE
        new_T[i] = new_T[i] + RATE
        Z1 = tZp[i]

    return new_T, new_r, pSoilLayerDepth


def test_STMPsimCalculator_case0() -> None:
    """
    Derived test from original fillTestVariables (aTestIndex=0).
    Verifies SoilTempArray after one process step.
    """
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0
    iDoInitialize = False

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        iDoInitialize=iDoInitialize,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
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
    assert len(SoilTempArray) == len(expected)
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) < 1e-6, f"{a} != {b}"