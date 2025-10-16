def init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialization function for STMPsimCalculator.

    Inputs:
    - cSoilLayerDepth: list of depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: mean air temperature on first day (°C)
    - cDampingDepth: initial damping depth of soil (m)
    - cAVT: constant temperature of deepest soil layer (long-term mean air temperature, °C)

    Returns:
    - SoilTempArray: list of initial soil temperatures per layer (°C)
    - rSoilTempArrayRate: list of initial daily temperature change rates per layer (°C/day), initialized to 0.0
    - pSoilLayerDepth: list of depths to the bottom of each (possibly extended) soil layer (m)
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
            depth = tProfileDepth + firstAdditionalLayerHight + (i - len(Z))
        tz[i] = depth
        # Linear approximation of soil temperature as initial value
        # STMP(depth) = (T0*(DD-depth) + AVT*depth) / DD
        # where T0 is cFirstDayMeanTemp and DD is cDampingDepth at initialization
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    SoilTempArray = tStmp
    rSoilTempArrayRate = tStmpRate
    pSoilLayerDepth = tz

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def process(iDoInitialize,
            cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cABD, cAVT,
            iSoilWaterContent, iSoilSurfaceTemperature,
            SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth):
    """
    Main biophysical process function for daily soil temperature update (STMPsimCalculator).

    Inputs:
    - iDoInitialize: boolean flag; if True, re-initializes state using init() before processing
    - cSoilLayerDepth: list of depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: mean air temperature on first day (°C) [used if iDoInitialize is True]
    - cDampingDepth: initial damping depth of soil (m) [used if iDoInitialize is True]
    - cABD: mean bulk density (t m-3)
    - cAVT: long-term average annual air temperature (°C)
    - iSoilWaterContent: soil water stored in the profile (mm)
    - iSoilSurfaceTemperature: soil surface temperature, DST0 (°C)
    - SoilTempArray: current soil temperature per layer (°C)
    - rSoilTempArrayRate: current daily temperature change per layer (°C/day)
    - pSoilLayerDepth: depth to bottom of (possibly extended) layers (m)

    Returns:
    - SoilTempArray_new: updated soil temperatures per layer (°C)
    - rSoilTempArrayRate_new: updated daily temperature change per layer (°C/day)
    - pSoilLayerDepth_new: (possibly re-initialized) depth to bottom of layers (m)
    """
    from math import exp, log

    # Optionally re-initialize state
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT
        )

    # Copy inputs to avoid mutating provided references (purity)
    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)
    SoilTempArray_new = list(SoilTempArray)
    rSoilTempArrayRate_new = [0.0] * len(SoilTempArray_new)
    pSoilLayerDepth_new = list(tZp)

    XLAG = 0.8  # Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG

    # Maximum damping depth (m)
    DP = (1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD)))

    # Water content related factor (dimensionless)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])

    # Damping depth (m)
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0  # Depth to bottom of previous layer
    for i in range(len(SoilTempArray_new)):
        # Depth factor: middle of layer depth divided by damping depth
        ZD = 0.5 * (Z1 + tZp[i]) / DD

        # FZ factor from APEX
        FZ = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))

        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray_new[i])

        Z1 = tZp[i]

        rSoilTempArrayRate_new[i] = RATE
        SoilTempArray_new[i] = SoilTempArray_new[i] + RATE

    return SoilTempArray_new, rSoilTempArrayRate_new, pSoilLayerDepth_new


def test_STMPsimCalculator_case0():
    """
    Derived test from STMPsimCalculator.fillTestVariables test case 0.
    Initializes, runs one process step, and checks SoilTempArray against expected reference values.
    """
    # Constants and inputs (DEFINE step)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT
    )

    # Process one timestep
    SoilTempArray_out, rSoilTempArrayRate_out, pSoilLayerDepth_out = process(
        iDoInitialize=False,
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cABD=cABD,
        cAVT=cAVT,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth
    )

    # Expected after process (CHECK step)
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

    assert len(SoilTempArray_out) == len(expected_SoilTempArray), "Layer count mismatch"
    for a, b in zip(SoilTempArray_out, expected_SoilTempArray):
        assert abs(a - b) < 1e-9, f"SoilTempArray mismatch: got {a}, expected {b}"

    return True