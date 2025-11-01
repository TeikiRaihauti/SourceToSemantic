from typing import List, Tuple
import math


def init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cAVT: float,
    cDampingDepth: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize STMPsimCalculator state.

    Inputs:
    - cSoilLayerDepth: List[float] - Depth to the bottom of each soil layer Z(l) [m]
    - cFirstDayMeanTemp: float - Mean air temperature on first day [°C]
    - cAVT: float - Long-term average annual air temperature (deepest soil layer temp) [°C]
    - cDampingDepth: float - Initial damping depth of soil [m]

    Returns:
    - SoilTempArray: List[float] - Initial soil temperature per layer [°C]
    - rSoilTempArrayRate: List[float] - Initial daily temperature change rate per layer [°C d-1] (zeros)
    - pSoilLayerDepth: List[float] - Depths to the bottom of the layers extended to damping depth [m]
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp
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
        # Linear approximation to the soil temperature as initial value
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz
    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def reInitialize(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cAVT: float,
    cDampingDepth: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Re-initialize state. See init for parameter descriptions and returns.
    """
    return init(cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth)


def process(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cAVT: float,
    cABD: float,
    cDampingDepth: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
    iDoInitialize: bool,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Main daily soil temperature process.

    Inputs:
    - cSoilLayerDepth: List[float] - Depth to the bottom of each soil layer Z(l) [m]
    - cFirstDayMeanTemp: float - Mean air temperature on first day [°C]
    - cAVT: float - Long-term average annual air temperature (deepest soil layer temp) [°C]
    - cABD: float - Average bulk density over all layers [t m-3]
    - cDampingDepth: float - Damping depth of soil [m]
    - iSoilWaterContent: float - Soil water stored in the soil profile [mm]
    - iSoilSurfaceTemperature: float - Soil surface temperature [°C]
    - SoilTempArray: List[float] - Current soil temperature per layer [°C]
    - rSoilTempArrayRate: List[float] - Current daily temperature change rate per layer [°C d-1]
    - pSoilLayerDepth: List[float] - Depths to the bottom of the layers extended to damping depth [m]
    - iDoInitialize: bool - If True, re-initialize state before processing

    Returns:
    - SoilTempArray: List[float] - Updated soil temperature per layer [°C]
    - rSoilTempArrayRate: List[float] - Updated daily temperature change rate per layer [°C d-1]
    - pSoilLayerDepth: List[float] - Layer depths (unchanged unless reinitialized) [m]
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = reInitialize(
            cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth
        )

    # Use local copies to preserve function purity
    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)
    prev_SoilTempArray = list(SoilTempArray)

    # Coefficients and derived parameters
    XLAG = 0.8
    XLG1 = 1.0 - XLAG
    DP = 1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    new_rSoilTempArrayRate: List[float] = [0.0] * len(prev_SoilTempArray)
    new_SoilTempArray: List[float] = [0.0] * len(prev_SoilTempArray)

    Z1 = 0.0
    for i in range(len(prev_SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - prev_SoilTempArray[i])
        Z1 = tZp[i]

        new_rSoilTempArrayRate[i] = RATE
        new_SoilTempArray[i] = prev_SoilTempArray[i] + RATE

    return new_SoilTempArray, new_rSoilTempArrayRate, tZp


def test_STMPsimCalculator_case0() -> bool:
    """
    Test derived from STMPsimCalculator.fillTestVariables for a single day update.
    Returns True if the test passes, otherwise raises AssertionError.
    """
    # Inputs (DEFINE)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize state
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth
    )

    # Process
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cABD=cABD,
        cDampingDepth=cDampingDepth,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        iDoInitialize=False,
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
        assert abs(a - b) <= 1e-9, f"Expected {b}, got {a}"

    return True