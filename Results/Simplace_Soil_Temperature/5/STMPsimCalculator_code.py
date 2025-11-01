def init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and layer depths.

    Parameters
    - cSoilLayerDepth: list[float] - Depth to the bottom of each soil layer (m), ascending, length >= 1.
    - cFirstDayMeanTemp: float - Mean air temperature on first day (°C).
    - cDampingDepth: float - Initial damping depth of soil (m).
    - cAVT: float - Constant temperature at damping depth (°C), long-term mean air temperature.

    Returns
    - SoilTempArray: list[float] - Initialized soil temperatures at each (possibly extended) layer center (°C).
    - rSoilTempArrayRate: list[float] - Initialized daily temperature change rates (°C/day), zeros.
    - pSoilLayerDepth: list[float] - Depth to the bottom of each (possibly extended) layer (m).
    """
    from math import ceil, floor

    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    additionalDepth = cDampingDepth - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(abs(ceil(additionalDepth))) + len(Z)

    SoilTempArray = []
    rSoilTempArrayRate = []
    pSoilLayerDepth = []

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            # first additional layer might be smaller than 1 m
            depth = tProfileDepth + firstAdditionalLayerHight + i - len(Z)
        pSoilLayerDepth.append(depth)
        # Linear approximation from surface (firstDayMeanTemp) to damping depth (AVT)
        SoilTempArray.append(
            (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth
        )
        rSoilTempArrayRate.append(0.0)

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def process(
    iDoInitialize,
    cSoilLayerDepth,
    cFirstDayMeanTemp,
    cAVT,
    cABD,
    cDampingDepth,
    iSoilWaterContent,
    iSoilSurfaceTemperature,
    SoilTempArray,
    rSoilTempArrayRate,
    pSoilLayerDepth,
):
    """
    Perform one daily update of the soil temperature profile.

    Parameters
    - iDoInitialize: bool - If True, re-initialize state before processing.
    - cSoilLayerDepth: list[float] - Depth to bottom of each soil layer (m).
    - cFirstDayMeanTemp: float - Mean air temperature on first day (°C), used only if initializing.
    - cAVT: float - Long-term mean air temperature (°C).
    - cABD: float - Average soil bulk density of the profile (t m^-3).
    - cDampingDepth: float - Initial (maximum) damping depth (m), used for initialization only.
    - iSoilWaterContent: float - Water stored in the soil profile at day t (mm).
    - iSoilSurfaceTemperature: float - Current soil surface temperature (°C).
    - SoilTempArray: list[float] - Current soil temperatures by layer (°C).
    - rSoilTempArrayRate: list[float] - Current rates (°C/day); values will be replaced by newly computed rates.
    - pSoilLayerDepth: list[float] - Depth to the bottom of each (possibly extended) soil layer (m).

    Returns
    - SoilTempArray: list[float] - Updated soil temperatures by layer (°C).
    - rSoilTempArrayRate: list[float] - Updated daily temperature change rates (°C/day).
    - pSoilLayerDepth: list[float] - Possibly re-initialized layer depths (m).
    """
    from math import exp, log

    # Re-initialize if requested
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )
    else:
        SoilTempArray = list(SoilTempArray)
        # Ensure rate vector matches the current profile size
        if rSoilTempArrayRate is None or len(rSoilTempArrayRate) != len(SoilTempArray):
            rSoilTempArrayRate = [0.0] * len(SoilTempArray)
        else:
            rSoilTempArrayRate = list(rSoilTempArrayRate)
        pSoilLayerDepth = list(pSoilLayerDepth)

    tZp = pSoilLayerDepth
    tZc = list(cSoilLayerDepth)

    # Coefficients
    XLAG = 0.8  # Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG

    # Damping depth (DD) dynamics based on bulk density and water content
    DP = 1.0 + 2.5 * cABD / (cABD + exp(6.53 - 5.63 * cABD))  # Maximum damping depth (m)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    # Layer-wise update
    Z1 = 0.0  # Depth to the bottom of the previous layer (m)
    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Depth factor relative to damping depth
        FZ = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray[i])  # °C/day
        rSoilTempArrayRate[i] = RATE
        SoilTempArray[i] = SoilTempArray[i] + RATE
        Z1 = tZp[i]

    return SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth


def test_STMPsimCalculator_case0():
    """
    Derived test from STMPsimCalculator.fillTestVariables (aTestIndex == 0).

    Setup:
    - cSoilLayerDepth = [0.1, 0.5, 1.5]
    - cFirstDayMeanTemp = 15.0
    - cAVT = 9.0
    - cABD = 1.4
    - cDampingDepth = 6.0
    - iSoilWaterContent = 0.3
    - iSoilSurfaceTemperature = 6.0

    Expected after one process step (after init):
    SoilTempArray ≈ [13.624360856350041, 13.399968504634286, 12.599999999999845,
                     12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    """
    import math

    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = process(
        iDoInitialize=False,
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cAVT=cAVT,
        cABD=cABD,
        cDampingDepth=cDampingDepth,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
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
    for a, b in zip(SoilTempArray, expected):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)