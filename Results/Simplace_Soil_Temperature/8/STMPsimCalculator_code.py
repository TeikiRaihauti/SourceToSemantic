import math

def STMPsimCalculator_initialize(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and layer depths up to damping depth.

    Inputs:
    - cSoilLayerDepth: list of floats, depths to the bottom of each initial soil layer (m)
    - cFirstDayMeanTemp: float, mean air temperature on first day (°C)
    - cDampingDepth: float, damping depth of soil (m)
    - cAVT: float, long-term average air temperature / deep soil temperature (°C)

    Returns:
    - SoilTempArray: list of floats, initial soil temperature profile (°C)
    - rSoilTempArrayRate: list of floats, initial rate of change (°C/d), zeros
    - pSoilLayerDepth: list of floats, depths to the bottom of each layer including additional layers (m)
    """
    Z = list(cSoilLayerDepth)
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
        tz[i] = depth
        # Set linear approximation to the soil temperature as initial value
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz
    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def STMPsimCalculator_process(SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth,
                              cSoilLayerDepth, cABD, iSoilWaterContent, iSoilSurfaceTemperature,
                              cAVT, iDoInitialize, cFirstDayMeanTemp, cDampingDepth):
    """
    Main daily soil temperature process.

    Inputs:
    - SoilTempArray: list of floats, current soil temperatures (°C)
    - rSoilTempArrayRate: list of floats, rate of change from previous step (°C/d)
    - pSoilLayerDepth: list of floats, depths to bottom of each layer including additional layers (m)
    - cSoilLayerDepth: list of floats, original configured depths to bottom of each soil layer (m)
    - cABD: float, mean bulk density (t m-3)
    - iSoilWaterContent: float, water content summed over the whole soil profile (mm)
    - iSoilSurfaceTemperature: float, soil surface temperature (°C)
    - cAVT: float, long-term average/deep soil temperature (°C)
    - iDoInitialize: bool, if True re-initialize arrays using initial conditions before processing
    - cFirstDayMeanTemp: float, mean air temperature on first day (°C) (used only if iDoInitialize is True)
    - cDampingDepth: float, damping depth (m) (used only if iDoInitialize is True)

    Returns:
    - SoilTempArray: list of floats, updated soil temperatures (°C)
    - rSoilTempArrayRate: list of floats, rate of change applied today (°C/d)
    - pSoilLayerDepth: list of floats, depths to bottom of each layer (m)
    """
    # Re-initialize if requested
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_initialize(
            cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
        )

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)
    XLAG = 0.8  # Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG

    # Damping depth calculation
    DP = (1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD)))  # Maximum damping depth (m)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP  # multiplied by 2 to increase damping

    Z1 = 0.0
    new_rates = [0.0] * len(SoilTempArray)
    new_temps = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Middle of layer depth divided by damping depth
        FZ = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = tZp[i]
        new_rates[i] = RATE
        new_temps[i] = SoilTempArray[i] + RATE

    return new_temps, new_rates, pSoilLayerDepth


def test_STMPsimCalculator_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables test case aTestIndex == 0.
    """
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_initialize(
        cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
    )

    # Process one day (no re-initialization during process)
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth,
        cSoilLayerDepth, cABD, iSoilWaterContent, iSoilSurfaceTemperature,
        cAVT, False, cFirstDayMeanTemp, cDampingDepth
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    assert len(SoilTempArray) == len(expected), "Length mismatch in SoilTempArray"
    for a, e in zip(SoilTempArray, expected):
        assert abs(a - e) < 1e-6, f"SoilTempArray mismatch: got {a}, expected {e}"