def STMPsimCalculator_initialize(cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth):
    """
    Initialize STMPsimCalculator state.

    Inputs:
    - cSoilLayerDepth: list[float] Depth to the bottom of each soil layer (m). Example: [0.1, 0.5, 1.5]
    - cFirstDayMeanTemp: float Mean air temperature on first day (deg C)
    - cAVT: float Constant temperature at damping depth (deg C), typically long-term average annual air temperature
    - cDampingDepth: float Initial damping depth of soil (m)

    Returns:
    - SoilTempArray: list[float] Initial soil temperature for each (possibly extended) soil layer (deg C)
    - rSoilTempArrayRate: list[float] Initial daily temperature change rates for each layer (deg C d-1), zeros
    - pSoilLayerDepth: list[float] Depths of layers including additional layers down to cDampingDepth (m)
    """
    from math import floor, ceil, fabs

    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
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
            # First additional layer might be smaller than 1 m
            depth = tProfileDepth + firstAdditionalLayerHight + (i - len(Z))
        tz[i] = depth
        # Linear approximation of soil temperature from surface to damping depth
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def STMPsimCalculator_process(
    cSoilLayerDepth,
    cAVT,
    cABD,
    iSoilWaterContent,
    iSoilSurfaceTemperature,
    SoilTempArray,
    pSoilLayerDepth,
    iDoInitialize=False,
    cFirstDayMeanTemp=0.0,
    cDampingDepth=0.0,
):
    """
    Main daily soil temperature process.

    Inputs:
    - cSoilLayerDepth: list[float] Depth to the bottom of each original soil layer (m)
    - cAVT: float Constant temperature at damping depth (deg C), typically long-term average annual air temperature
    - cABD: float Mean bulk density (t m-3)
    - iSoilWaterContent: float Water stored in the soil profile (mm)
    - iSoilSurfaceTemperature: float Soil surface temperature (deg C)
    - SoilTempArray: list[float] Current soil temperature profile (deg C) [state]
    - pSoilLayerDepth: list[float] Depth to bottom of each layer including additional layers (m) [state]
    - iDoInitialize: bool If True, re-initialize state using provided initial parameters
    - cFirstDayMeanTemp: float Required only if iDoInitialize is True; Mean air temperature on first day (deg C)
    - cDampingDepth: float Required only if iDoInitialize is True; Initial damping depth of soil (m)

    Returns:
    - SoilTempArray: list[float] Updated soil temperature profile (deg C)
    - rSoilTempArrayRate: list[float] Daily temperature change rates for each layer (deg C d-1)
    - pSoilLayerDepth: list[float] Depth to bottom of each layer including additional layers (m)
    """
    from math import exp, log

    # Re-initialize if requested
    if iDoInitialize:
        SoilTempArray, _, pSoilLayerDepth = STMPsimCalculator_initialize(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cAVT=cAVT,
            cDampingDepth=cDampingDepth,
        )

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)

    XLAG = 0.8  # Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))  # Maximum damping depth (m)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP  # Damping depth (m)

    Z1 = 0.0  # Depth to the bottom of the previous soil layer (m)
    new_SoilTempArray = list(SoilTempArray)
    rSoilTempArrayRate = [0.0] * len(new_SoilTempArray)

    for i in range(len(new_SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Depth factor
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - new_SoilTempArray[i])  # Daily change (deg C d-1)
        Z1 = tZp[i]
        rSoilTempArrayRate[i] = RATE
        new_SoilTempArray[i] = new_SoilTempArray[i] + RATE

    return new_SoilTempArray, rSoilTempArrayRate, tZp


def test_STMPsimCalculator_0():
    """
    Derived from STMPsimCalculator.fillTestVariables (aTestIndex=0).

    Steps:
    - Initialize with:
        cSoilLayerDepth = [0.1, 0.5, 1.5]
        cFirstDayMeanTemp = 15.0
        cAVT = 9.0
        cDampingDepth = 6.0
    - Process with:
        cABD = 1.4
        iSoilWaterContent = 0.3
        iSoilSurfaceTemperature = 6.0

    Expects SoilTempArray after one process step:
    [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    """
    # Define inputs
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cDampingDepth = 6.0
    cABD = 1.4
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_initialize(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cDampingDepth=cDampingDepth,
    )

    # Process one time step
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        pSoilLayerDepth=pSoilLayerDepth,
        iDoInitialize=False,
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    assert len(SoilTempArray) == len(expected)
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) < 1e-9

    return {
        "SoilTempArray": SoilTempArray,
        "rSoilTempArrayRate": rSoilTempArrayRate,
        "pSoilLayerDepth": pSoilLayerDepth,
        "passed": True,
    }