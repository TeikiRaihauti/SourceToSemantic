from math import floor, ceil, fabs, exp, log


def STMPsimCalculator_init(cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cDampingDepth):
    """
    Initialize soil temperature profile and layer depths for STMP.

    Inputs:
    - cSoilLayerDepth: list of depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: mean air temperature on first day (degC)
    - cAVT: long-term average annual air temperature (degC)
    - cDampingDepth: initial damping depth of soil (m)

    Returns:
    - SoilTempArray: list of initial soil temperatures per layer (degC)
    - rSoilTempArrayRate: list of initial daily temperature change rates per layer (degC/day), initialized to 0.0
    - pSoilLayerDepth: list of layer bottom depths including any added layers to reach the damping depth (m)
    """
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
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = depth
        # Linear approximation to soil temperature as initial value
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def STMPsimCalculator_process(cSoilLayerDepth, cFirstDayMeanTemp, cAVT, cABD, cDampingDepth,
                              iSoilWaterContent, iSoilSurfaceTemperature, iDoInitialize,
                              pSoilLayerDepth, SoilTempArray):
    """
    One time-step update of soil temperature in layers.

    Inputs:
    - cSoilLayerDepth: list of depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: mean air temperature on first day (degC), used only when iDoInitialize is True
    - cAVT: long-term average annual air temperature (degC)
    - cABD: average soil bulk density of the profile (t m-3)
    - cDampingDepth: damping depth of soil (m), used only when iDoInitialize is True
    - iSoilWaterContent: water stored in the profile (mm)
    - iSoilSurfaceTemperature: soil surface temperature (degC)
    - iDoInitialize: if True, re-initialize the model before processing this time step
    - pSoilLayerDepth: list of layer bottom depths including any added layers to reach the damping depth (m)
    - SoilTempArray: list of yesterday's soil temperatures per layer (degC)

    Returns:
    - SoilTempArray: updated list of soil temperatures per layer (degC)
    - rSoilTempArrayRate: list of daily temperature change rates per layer (degC/day)
    - pSoilLayerDepth: potentially re-initialized list of layer depths (m)
    """
    # Optional re-initialization
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cAVT=cAVT,
            cDampingDepth=cDampingDepth
        )
    else:
        rSoilTempArrayRate = [0.0] * len(SoilTempArray)

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)

    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    updated_SoilTempArray = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = tZp[i]
        rSoilTempArrayRate[i] = RATE
        updated_SoilTempArray[i] = SoilTempArray[i] + RATE

    return updated_SoilTempArray, rSoilTempArrayRate, tZp


def test_STMPsimCalculator_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables case 0.
    """
    # Inputs before process
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0
    iDoInitialize = False

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cDampingDepth=cDampingDepth
    )

    # Process one step
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cABD=cABD,
        cDampingDepth=cDampingDepth,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        iDoInitialize=iDoInitialize,
        pSoilLayerDepth=pSoilLayerDepth,
        SoilTempArray=SoilTempArray
    )

    expected_SoilTempArray = [
        13.624360856350041, 13.399968504634286, 12.599999999999845,
        12.2, 11.4, 10.6, 9.799999999999999, 9.0
    ]

    # Assert approximate equality
    for a, b in zip(SoilTempArray, expected_SoilTempArray):
        assert abs(a - b) < 1e-6, f"Mismatch: {a} vs {b}"