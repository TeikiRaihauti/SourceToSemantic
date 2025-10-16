from math import floor, ceil, fabs, exp, log


def init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and layer depths.

    Inputs:
    - cSoilLayerDepth: list of depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: mean air temperature on first day (degC)
    - cDampingDepth: damping depth (m)
    - cAVT: long-term average annual air temperature (degC)

    Returns:
    - SoilTempArray: list of initialized soil temperatures per layer (degC)
    - rSoilTempArrayRate: list of daily soil temperature change rates per layer (degC/day), initialized to 0
    - pSoilLayerDepth: list of depths to the bottom of each (possibly extended) soil layer (m)
    """
    Z = cSoilLayerDepth
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
        # Linear approximation to soil temperature as initial value
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def process(cSoilLayerDepth,
            cAVT,
            cABD,
            iSoilWaterContent,
            iSoilSurfaceTemperature,
            SoilTempArray,
            rSoilTempArrayRate,
            pSoilLayerDepth,
            iDoInitialize,
            cFirstDayMeanTemp,
            cDampingDepth):
    """
    Update soil temperature profile for one time step.

    Inputs:
    - cSoilLayerDepth: list of depths to the bottom of each native soil layer (m)
    - cAVT: long-term average annual air temperature (degC)
    - cABD: average soil bulk density (t m-3)
    - iSoilWaterContent: water stored in the soil profile (mm)
    - iSoilSurfaceTemperature: soil surface temperature (degC)
    - SoilTempArray: current soil temperature profile (degC)
    - rSoilTempArrayRate: current rates of change (degC/day)
    - pSoilLayerDepth: list of depths to the bottom of each (possibly extended) soil layer (m)
    - iDoInitialize: boolean flag to re-initialize before processing
    - cFirstDayMeanTemp: mean air temperature on first day (degC), used for re-initialization
    - cDampingDepth: damping depth (m), used for re-initialization

    Returns:
    - SoilTempArray: updated soil temperature profile (degC)
    - rSoilTempArrayRate: updated rates of change (degC/day)
    - pSoilLayerDepth: possibly updated layer depths (m) if re-initialized
    """
    # Optionally re-initialize state
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT
        )

    tZp = pSoilLayerDepth
    tZc = cSoilLayerDepth

    # Parameters and derived variables
    XLAG = 0.8  # lag coefficient
    XLG1 = 1.0 - XLAG
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))  # maximum damping depth (m)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP  # damping depth (m)

    Z1 = 0.0  # depth of bottom of previous soil layer
    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # depth factor
        RATE = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])  # rate of change of temperature at layer i

        Z1 = tZp[i]
        rSoilTempArrayRate[i] = RATE
        SoilTempArray[i] = SoilTempArray[i] + RATE

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


# Tests derived from STMPsimCalculator.fillTestVariables
def test_STMPsimCalculator_case0():
    # Define inputs (from STMPsimCalculator.fillTestVariables, aTestIndex = 0)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0
    iDoInitialize = False  # no re-initialization before process

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT
    )

    # Process one step
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        cSoilLayerDepth=cSoilLayerDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        iDoInitialize=iDoInitialize,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth
    )

    # Expected output after one process call
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

    assert len(SoilTempArray) == len(expected_SoilTempArray)
    for a, b in zip(SoilTempArray, expected_SoilTempArray):
        assert abs(a - b) < 1e-9, f"Mismatch: got {a}, expected {b}"