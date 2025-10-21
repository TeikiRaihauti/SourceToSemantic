from typing import List, Tuple
import math


def STMPsimCalculator_init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    cAVT: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize soil temperature profile and depths.

    Inputs:
    - cSoilLayerDepth: List[float] - Depth to the bottom of each soil layer (m).
    - cFirstDayMeanTemp: float - Mean air temperature on first day (°C).
    - cDampingDepth: float - Initial damping depth for soil (m).
    - cAVT: float - Long-term average annual air temperature (°C).

    Returns:
    - SoilTempArray: List[float] - Initial soil temperature per layer (°C).
    - rSoilTempArrayRate: List[float] - Initial daily temperature change per layer (°C/day), initialized to 0.0.
    - pSoilLayerDepth: List[float] - Depth to the bottom of each (including added) soil layer (m).
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - math.floor(additionalDepth)
    layers = int(abs(math.ceil(additionalDepth))) + len(Z)

    tStmp: List[float] = [0.0] * layers
    tStmpRate: List[float] = [0.0] * layers
    tz: List[float] = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = depth
        # Linear approximation for initial temperature profile
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def _compute_damping_depth(
    cABD: float,
    iSoilWaterContent: float,
    profile_depth: float,
) -> float:
    """
    Compute damping depth (DD) based on bulk density and water content.

    Inputs:
    - cABD: float - Average soil bulk density of profile (t m^-3).
    - iSoilWaterContent: float - Water stored in the profile (mm).
    - profile_depth: float - Depth to the bottom of the last original soil layer (m).

    Returns:
    - DD: float - Damping depth (m).
    """
    DP = 1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * profile_depth)
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP
    return DD


def STMPsimCalculator_process(
    iDoInitialize: bool,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    cSoilLayerDepth: List[float],
    cAVT: float,
    cABD: float,
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Update daily soil temperatures across layers.

    Inputs:
    - iDoInitialize: bool - If True, re-initialize state before processing.
    - iSoilWaterContent: float - Water content of the whole soil profile (mm).
    - iSoilSurfaceTemperature: float - Soil surface temperature (°C).
    - cSoilLayerDepth: List[float] - Original depth to bottom of each soil layer (m).
    - cAVT: float - Long-term average annual air temperature (°C).
    - cABD: float - Average soil bulk density of the profile (t m^-3).
    - SoilTempArray: List[float] - Current soil temperature per layer (°C).
    - rSoilTempArrayRate: List[float] - Current daily soil temperature change per layer (°C/day).
    - pSoilLayerDepth: List[float] - Depth to bottom of each (including added) soil layer (m).
    - cFirstDayMeanTemp: float - Mean air temperature on first day (°C) for re-initialization.
    - cDampingDepth: float - Damping depth for initialization (m).

    Returns:
    - SoilTempArray: List[float] - Updated soil temperatures (°C).
    - rSoilTempArrayRate: List[float] - Updated daily temperature change per layer (°C/day).
    - pSoilLayerDepth: List[float] - (Potentially re-initialized) depths of layers (m).
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )

    # Ensure rate array matches current number of layers
    if len(rSoilTempArrayRate) != len(SoilTempArray):
        rSoilTempArrayRate = [0.0] * len(SoilTempArray)

    # Parameters
    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    # Damping depth based on current conditions
    DD = _compute_damping_depth(
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        profile_depth=cSoilLayerDepth[-1],
    )

    Z1 = 0.0
    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + pSoilLayerDepth[i]) / DD
        FZ = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = pSoilLayerDepth[i]
        rSoilTempArrayRate[i] = RATE
        SoilTempArray[i] = SoilTempArray[i] + RATE

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def test_STMPsimCalculator_case0() -> None:
    """
    Derived test from STMPsimCalculator.fillTestVariables(aTestIndex=0).
    Initializes then processes one day and checks SoilTempArray against expected values.
    """
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        iDoInitialize=False,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    assert len(SoilTempArray) == len(expected)
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) < 1e-9, f"Got {a}, expected {b}"