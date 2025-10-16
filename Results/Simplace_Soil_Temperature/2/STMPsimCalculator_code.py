def stmp_initialize(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and extended layer depths up to damping depth.

    Inputs:
    - cSoilLayerDepth: list of layer bottom depths (m), increasing, profile bottom at last element
    - cFirstDayMeanTemp: mean air temperature on first day (degC)
    - cDampingDepth: damping depth (m)
    - cAVT: long-term average annual air temperature (degC)

    Returns:
    - SoilTempArray: list of initial soil temperatures for each layer (degC)
    - rSoilTempArrayRate: list of initial daily temperature change rates for each layer (degC/day), zeros
    - pSoilLayerDepth: list of extended layer bottom depths including additional layers to reach damping depth (m)
    """
    import math

    Z = list(cSoilLayerDepth)
    if len(Z) == 0:
        return [], [], []

    tProfileDepth = Z[-1]
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
        # Linear approximation between surface (cFirstDayMeanTemp at depth 0) and cAVT at damping depth
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def stmp_process(
    cSoilLayerDepth,
    cAVT,
    cABD,
    iSoilWaterContent,
    iSoilSurfaceTemperature,
    SoilTempArray,
    pSoilLayerDepth,
    iDoInitialize,
    cFirstDayMeanTemp,
    cDampingDepth
):
    """
    Advance soil temperatures one day using the APEX soil temperature approach (Williams & Izaurralde, 2005).

    Inputs:
    - cSoilLayerDepth: list of layer bottom depths (m) (original configuration layers)
    - cAVT: long-term average annual air temperature (degC)
    - cABD: average soil bulk density (t m-3)
    - iSoilWaterContent: total water stored in profile (mm)
    - iSoilSurfaceTemperature: current soil surface temperature (degC)
    - SoilTempArray: current soil temperatures (degC) for each layer (state before update)
    - pSoilLayerDepth: extended layer depths (m) corresponding to SoilTempArray
    - iDoInitialize: boolean flag; if True, re-initialize states before processing this day
    - cFirstDayMeanTemp: mean air temperature on first day (degC) for re-initialization
    - cDampingDepth: damping depth (m) for re-initialization

    Returns:
    - SoilTempArray_new: updated soil temperatures (degC)
    - rSoilTempArrayRate_new: daily temperature change for each layer (degC/day)
    - pSoilLayerDepth_new: extended depths (m) used for this step (may change if re-initialized)
    """
    import math

    # Optionally re-initialize
    if iDoInitialize:
        SoilTempArray, _, pSoilLayerDepth = stmp_initialize(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)

    if len(SoilTempArray) == 0 or len(tZp) == 0:
        return list(SoilTempArray), [0.0] * len(SoilTempArray), list(pSoilLayerDepth)

    # Parameters
    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    # Damping depth terms
    DP = 1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    # Layer-by-layer update
    Z1 = 0.0
    rSoilTempArrayRate_new = [0.0] * len(SoilTempArray)
    SoilTempArray_new = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        # FZ = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        denom = ZD + math.exp(-0.8669 - 2.0775 * ZD)
        FZ = ZD / denom if denom != 0.0 else 0.0
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])
        Z1 = tZp[i]
        rSoilTempArrayRate_new[i] = RATE
        SoilTempArray_new[i] = SoilTempArray[i] + RATE

    return SoilTempArray_new, rSoilTempArrayRate_new, tZp


def test_stmp_process_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables test case (aTestIndex == 0).
    """
    # Inputs before process (assume init occurs prior to process step)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = stmp_initialize(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    # Process one day
    SoilTempArray_new, rSoilTempArrayRate_new, pSoilLayerDepth_new = stmp_process(
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        pSoilLayerDepth=pSoilLayerDepth,
        iDoInitialize=False,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    assert len(SoilTempArray_new) == len(expected)
    for a, b in zip(SoilTempArray_new, expected):
        assert abs(a - b) < 1e-9, f"{a} != {b}"