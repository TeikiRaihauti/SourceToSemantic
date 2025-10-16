from math import exp, log as ln, floor, ceil

def STMPsimCalculator_init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialization function for STMPsimCalculator.

    Inputs:
    - cSoilLayerDepth: list of doubles, depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: double, mean air temperature on first day (°C)
    - cDampingDepth: double, initial value for damping depth of soil (m)
    - cAVT: double, constant temperature of deepest soil layer, use long term mean air temperature (°C)

    Returns:
    - SoilTempArray: list of doubles, initial soil temperatures in layers (°C)
    - rSoilTempArrayRate: list of doubles, initial daily temperature change (°C/d), initialized to zeros
    - pSoilLayerDepth: list of doubles, depth to bottom of soil layers including added layers to reach damping depth (m)
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(abs(ceil(additionalDepth))) + len(Z)
    tStmp = [0.0] * layers
    tStmpRate = [0.0] * layers
    tz = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        tz[i] = depth
        # Linear approximation for initial temperature profile
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def STMPsimCalculator_process(iDoInitialize,
                              cSoilLayerDepth,
                              cFirstDayMeanTemp,
                              cDampingDepth,
                              cAVT,
                              cABD,
                              iSoilWaterContent,
                              iSoilSurfaceTemperature,
                              SoilTempArray,
                              rSoilTempArrayRate,
                              pSoilLayerDepth):
    """
    Main biophysical process function for STMPsimCalculator.

    Inputs:
    - iDoInitialize: boolean, if True re-initialize state before processing
    - cSoilLayerDepth: list of doubles, depths to the bottom of each original soil layer (m)
    - cFirstDayMeanTemp: double, mean air temperature on first day (°C) (used only if iDoInitialize)
    - cDampingDepth: double, initial value for damping depth of soil (m) (used only if iDoInitialize)
    - cAVT: double, constant temperature of deepest soil layer, use long term mean air temperature (°C)
    - cABD: double, mean bulk density (t m-3)
    - iSoilWaterContent: double, water content, sum of whole soil profile (mm)
    - iSoilSurfaceTemperature: double, soil surface temperature (°C)
    - SoilTempArray: list of doubles, previous soil temperatures in layers (°C)
    - rSoilTempArrayRate: list of doubles, previous daily temperature change (°C/d)
    - pSoilLayerDepth: list of doubles, depths to bottom of layers including extensions (m)

    Returns:
    - SoilTempArray: list of doubles, updated soil temperatures in layers (°C)
    - rSoilTempArrayRate: list of doubles, updated daily temperature change (°C/d)
    - pSoilLayerDepth: list of doubles, unchanged depths array
    """
    # Optional re-initialization
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT
        )

    # Use local copies to avoid mutating input lists directly
    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)
    curr_SoilTempArray = list(SoilTempArray)
    curr_rSoilTempArrayRate = [0.0] * len(tZp)

    # Constants
    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    # Damping depth calculations
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(ln(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    updated_SoilTempArray = [0.0] * len(tZp)
    for i in range(len(tZp)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        RATE = (ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - curr_SoilTempArray[i])
        Z1 = tZp[i]
        curr_rSoilTempArrayRate[i] = RATE
        updated_SoilTempArray[i] = curr_SoilTempArray[i] + RATE

    return updated_SoilTempArray, curr_rSoilTempArrayRate, tZp


def test_STMPsimCalculator_case0():
    """
    Derived test based on STMPsimCalculator.fillTestVariables for aTestIndex == 0.

    Procedure:
    - Define inputs as in the Java test 'DEFINE' part
    - Call init()
    - Call process()
    - Check SoilTempArray matches expected values (within tolerance)
    """
    # Define inputs
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
        cDampingDepth=cDampingDepth,
        cAVT=cAVT
    )

    # Process
    iDoInitialize = False
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        iDoInitialize=iDoInitialize,
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth
    )

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

    tol = 1e-9
    assert len(SoilTempArray) == len(expected_SoilTempArray), "Length mismatch in SoilTempArray"
    for i, (got, expv) in enumerate(zip(SoilTempArray, expected_SoilTempArray)):
        assert abs(got - expv) <= tol, f"SoilTempArray[{i}] = {got}, expected {expv}"