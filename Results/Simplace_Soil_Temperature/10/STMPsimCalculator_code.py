from math import floor, ceil, fabs, exp, log

def STMPsimCalculator_init(cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth):
    """
    Initialize soil temperature profile and layer depths up to the damping depth.

    Inputs:
    - cSoilLayerDepth: list of soil layer bottom depths (m)
    - cFirstDayMeanTemp: mean air temperature on first day (°C)
    - cAVT: long-term average annual air temperature (°C)
    - cDampingDepth: initial value for damping depth of soil (m)

    Returns:
    - SoilTempArray: list of initial soil temperatures in layers (°C)
    - rSoilTempArrayRate: list of initial daily temperature change rates (°C/day), initialized to 0
    - pSoilLayerDepth: list of soil layer depths including additional layers to reach damping depth (m)
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
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = depth
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz
    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def STMPsimCalculator_process(SoilTempArray, pSoilLayerDepth, cSoilLayerDepth, cABD, iSoilWaterContent, iSoilSurfaceTemperature, cAVT):
    """
    Update soil temperature profile for one time step.

    Inputs:
    - SoilTempArray: list of current soil temperatures in layers (°C)
    - pSoilLayerDepth: list of depths to the bottom of each layer including additional layers (m)
    - cSoilLayerDepth: list of original soil layer depths (m)
    - cABD: average soil bulk density of the profile (t m-3)
    - iSoilWaterContent: water stored in the profile (mm)
    - iSoilSurfaceTemperature: soil surface temperature (°C)
    - cAVT: long-term average annual air temperature (°C)

    Returns:
    - SoilTempArray: updated soil temperatures (°C)
    - rSoilTempArrayRate: list of daily temperature change rates (°C/day)
    """
    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)

    XLAG = 0.8
    XLG1 = 1.0 - XLAG
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    new_rates = [0.0] * len(SoilTempArray)
    new_temps = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = tZp[i]
        new_rates[i] = RATE
        new_temps[i] = SoilTempArray[i] + RATE

    return new_temps, new_rates


def test_STMPsimCalculator_process_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables test case aTestIndex == 0.
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
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cDampingDepth=cDampingDepth
    )

    # Process one step
    SoilTempArray, rSoilTempArrayRate = STMPsimCalculator_process(
        SoilTempArray=SoilTempArray,
        pSoilLayerDepth=pSoilLayerDepth,
        cSoilLayerDepth=cSoilLayerDepth,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        cAVT=cAVT
    )

    # Expected after process
    expected_SoilTempArray = [
        13.624360856350041,
        13.399968504634286,
        12.599999999999845,
        12.2,
        11.4,
        10.6,
        9.799999999999999,
        9.0
    ]

    # Assert with tolerance
    assert len(SoilTempArray) == len(expected_SoilTempArray)
    for a, b in zip(SoilTempArray, expected_SoilTempArray):
        assert abs(a - b) < 1e-6, f"Expected {b}, got {a}"