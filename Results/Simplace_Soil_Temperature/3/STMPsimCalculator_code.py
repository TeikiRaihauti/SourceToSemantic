def init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and derived layer depths for STMPsimCalculator.

    Inputs:
    - cSoilLayerDepth: list of depths (m) to the bottom of each existing soil layer (Z)
    - cFirstDayMeanTemp: mean air temperature on the first day (degC)
    - cDampingDepth: initial damping depth of soil (m)
    - cAVT: long-term average annual air temperature (degC)

    Returns:
    - SoilTempArray: initialized soil temperature at the center of layers (degC)
    - rSoilTempArrayRate: initialized daily temperature change (degC/day), zeros
    - pSoilLayerDepth: derived depths (m) including additional layers down to cDampingDepth
    """
    from math import floor, ceil, fabs

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
        # Linear approximation between surface (firstDayMeanTemp) and damping depth (cAVT)
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def _compute_damping_depth(cABD, iSoilWaterContent, cSoilLayerDepth):
    """
    Supporting function: compute DP, WC, and DD according to APEX equations.
    """
    from math import exp, log

    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    # WC related factor from total soil water content and profile depth
    profile_depth = cSoilLayerDepth[-1]
    denom = (0.356 - 0.144 * cABD) * profile_depth
    # Avoid division by zero (degenerate soils)
    WC = 0.0 if denom == 0.0 else 0.001 * iSoilWaterContent / denom
    # Multiply exponent by 2 (per original code comment)
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP
    return DP, WC, DD


def process(iDoInitialize,
            cSoilLayerDepth,
            cFirstDayMeanTemp,
            cAVT,
            cABD,
            cDampingDepth,
            iSoilWaterContent,
            iSoilSurfaceTemperature,
            SoilTempArray,
            pSoilLayerDepth):
    """
    Main daily process for STMPsimCalculator.

    Inputs:
    - iDoInitialize: boolean flag; if True, re-initialize states using init
    - cSoilLayerDepth: list of depths (m) to the bottom of each soil layer (Z)
    - cFirstDayMeanTemp: mean air temperature on the first day (degC)
    - cAVT: long-term average annual air temperature (degC)
    - cABD: average soil bulk density (t m-3)
    - cDampingDepth: initial damping depth for initialization (m)
    - iSoilWaterContent: water stored in the profile (mm)
    - iSoilSurfaceTemperature: soil surface temperature (degC)
    - SoilTempArray: current state soil temperature array (degC)
    - pSoilLayerDepth: derived depths (m) including additional layers

    Returns:
    - SoilTempArray_new: updated soil temperature array (degC)
    - rSoilTempArrayRate_new: daily temperature change array (degC/day)
    - pSoilLayerDepth_new: potentially updated layer depths (if re-initialized)
    """
    from math import exp

    # Re-initialize if requested
    if iDoInitialize:
        SoilTempArray, _, pSoilLayerDepth = init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT)

    tZp = pSoilLayerDepth
    tZc = cSoilLayerDepth

    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    _, _, DD = _compute_damping_depth(cABD, iSoilWaterContent, tZc)

    Z1 = 0.0  # Depth of the bottom of the previous soil layer (m)
    rSoilTempArrayRate_new = [0.0] * len(SoilTempArray)
    SoilTempArray_new = [0.0] * len(SoilTempArray)

    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Depth factor
        FZ = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])  # daily change
        Z1 = tZp[i]
        rSoilTempArrayRate_new[i] = RATE
        SoilTempArray_new[i] = SoilTempArray[i] + RATE

    return SoilTempArray_new, rSoilTempArrayRate_new, pSoilLayerDepth


def test_STMPsimCalculator_case0():
    """
    Derived test from STMPsimCalculator.fillTestVariables (aTestIndex == 0).
    Sets up inputs, runs init and one process step, and checks SoilTempArray.
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
        cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        iDoInitialize,
        cSoilLayerDepth,
        cFirstDayMeanTemp,
        cAVT,
        cABD,
        cDampingDepth,
        iSoilWaterContent,
        iSoilSurfaceTemperature,
        SoilTempArray,
        pSoilLayerDepth
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]

    # Simple assertion to verify results within a small tolerance
    assert len(SoilTempArray) == len(expected)
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) < 1e-9, f"Mismatch: {a} vs {b}"

    return SoilTempArray, expected