from typing import List, Tuple
import math


def init(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cAVT: float,
    cDampingDepth: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialize soil temperature profile and related arrays.

    Inputs:
    - cSoilLayerDepth: List[float], depths to the bottom of each soil layer [m]
    - cFirstDayMeanTemp: float, mean air temperature on first day [°C]
    - cAVT: float, constant temperature of deepest soil layer (long-term mean air temp) [°C]
    - cDampingDepth: float, initial damping depth of soil [m]

    Returns:
    - SoilTempArray: List[float], initialized soil temperatures per layer [°C]
    - rSoilTempArrayRate: List[float], initialized daily temperature change rates per layer [°C d-1] (zeros)
    - pSoilLayerDepth: List[float], depths to the bottom of layers including additional layers up to damping depth [m]
    """
    Z = cSoilLayerDepth
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - math.floor(additionalDepth)
    layers = int(abs(math.ceil(additionalDepth))) + len(Z)

    tStmp = [0.0] * layers
    tStmpRate = [0.0] * layers
    tz = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = float(depth)
        # Linear approximation from first day mean temp at surface to AVT at damping depth
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz
    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def process(
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
    cSoilLayerDepth: List[float],
    cAVT: float,
    cABD: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    iDoInitialize: bool,
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Update soil temperature profile for one time step.

    Inputs:
    - SoilTempArray: List[float], current soil temperatures per layer [°C]
    - rSoilTempArrayRate: List[float], current daily temperature change rates per layer [°C d-1] (ignored, recomputed)
    - pSoilLayerDepth: List[float], depths to the bottom of layers incl. additional layers [m]
    - cSoilLayerDepth: List[float], depths to the bottom of measured soil layers [m]
    - cAVT: float, constant temperature of deepest soil layer (long-term mean air temp) [°C]
    - cABD: float, average soil bulk density [t m-3]
    - iSoilWaterContent: float, sum of water content of the whole soil profile [mm]
    - iSoilSurfaceTemperature: float, soil surface temperature [°C]
    - iDoInitialize: bool, if True re-initialize state before processing
    - cFirstDayMeanTemp: float, mean air temperature on first day [°C] (needed when re-initializing)
    - cDampingDepth: float, damping depth of soil [m] (needed when re-initializing)

    Returns:
    - SoilTempArray: List[float], updated soil temperatures per layer [°C]
    - rSoilTempArrayRate: List[float], updated daily temperature change rates per layer [°C d-1]
    - pSoilLayerDepth: List[float], (potentially re-initialized) depths to the bottom of layers [m]
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cAVT=cAVT,
            cDampingDepth=cDampingDepth,
        )

    tZp = pSoilLayerDepth
    tZc = cSoilLayerDepth

    XLAG = 0.8  # Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG

    # Maximum damping depth (m)
    DP = 1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD))

    # Water content related factor (dimensionless)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])

    # Damping depth (m), multiplied by 2 in exponent as in original implementation
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0  # Depth of bottom of previous soil layer (m)
    new_rates: List[float] = [0.0] * len(SoilTempArray)
    new_temps: List[float] = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Depth factor at layer center
        RATE = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])  # Rate of change (°C d-1)
        Z1 = tZp[i]
        new_rates[i] = RATE
        new_temps[i] = SoilTempArray[i] + RATE

    return new_temps, new_rates, pSoilLayerDepth


def test_STMPsimCalculator_case0() -> None:
    """
    Test derived from STMPsimCalculator.fillTestVariables (aTestIndex = 0).
    """
    # Constants and inputs
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0
    iDoInitialize = False

    # Initialize state
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cDampingDepth=cDampingDepth,
    )

    # Process one time step
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        iDoInitialize=iDoInitialize,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
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
    tol = 1e-6
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) <= tol, f"Expected {b}, got {a}"