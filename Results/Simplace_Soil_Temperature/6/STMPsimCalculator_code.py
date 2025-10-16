import math

def init(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature profile and derived layer depths for STMP.

    Inputs:
    - cSoilLayerDepth: list of layer bottom depths (m)
    - cFirstDayMeanTemp: mean air temperature on first day (°C)
    - cDampingDepth: initial damping depth (m)
    - cAVT: long-term average annual air temperature (°C)

    Returns:
    - SoilTempArray: list of initial soil temperatures per layer (°C)
    - rSoilTempArrayRate: list of initial daily temperature change per layer (°C/day), initialized to 0
    - pSoilLayerDepth: list of profile depths including additional layers up to damping depth (m)
    """
    Z = list(cSoilLayerDepth)
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
        # Linear approximation across damping depth
        tStmp[i] = (cFirstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def process(cAVT, cABD, pSoilLayerDepth, cSoilLayerDepth, SoilTempArray, iSoilSurfaceTemperature,
            iSoilWaterContent, iDoInitialize, cFirstDayMeanTemp, cDampingDepth):
    """
    Daily soil temperature process for STMP.

    Inputs:
    - cAVT: long-term average annual air temperature (°C)
    - cABD: average soil bulk density (t m-3)
    - pSoilLayerDepth: list of profile depths including additional layers up to damping depth (m)
    - cSoilLayerDepth: list of soil layer bottom depths (m)
    - SoilTempArray: list of current soil temperatures per layer (°C)
    - iSoilSurfaceTemperature: soil surface temperature DST0 (°C)
    - iSoilWaterContent: total soil water stored in profile (mm)
    - iDoInitialize: boolean flag to re-initialize before processing
    - cFirstDayMeanTemp: mean air temperature on first day (°C) for re-initialization
    - cDampingDepth: initial damping depth (m) for re-initialization

    Returns:
    - SoilTempArray_new: updated soil temperatures per layer (°C)
    - rSoilTempArrayRate_new: daily temperature change per layer (°C/day)
    - pSoilLayerDepth_new: profile depths (m), possibly re-initialized
    """
    # Re-initialize if requested
    if iDoInitialize:
        SoilTempArray_init, rSoilTempArrayRate_init, pSoilLayerDepth_init = init(
            cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
        )
        SoilTempArray_proc = SoilTempArray_init
        pSoilLayerDepth_proc = pSoilLayerDepth_init
    else:
        SoilTempArray_proc = list(SoilTempArray)
        pSoilLayerDepth_proc = list(pSoilLayerDepth)

    tZp = pSoilLayerDepth_proc
    tZc = cSoilLayerDepth

    XLAG = 0.8
    XLG1 = 1.0 - XLAG

    DP = 1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD))
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP

    Z1 = 0.0
    rSoilTempArrayRate_new = [0.0] * len(SoilTempArray_proc)
    SoilTempArray_new = [0.0] * len(SoilTempArray_proc)

    for i in range(len(SoilTempArray_proc)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD
        FZ = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - SoilTempArray_proc[i])
        Z1 = tZp[i]
        rSoilTempArrayRate_new[i] = RATE
        SoilTempArray_new[i] = SoilTempArray_proc[i] + RATE

    return SoilTempArray_new, rSoilTempArrayRate_new, pSoilLayerDepth_proc


def test_STMPsimCalculator_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables, test index 0.
    Initializes the component and runs one process step; checks SoilTempArray.
    """
    # Define inputs (before process)
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0
    iDoInitialize = False

    # Initialize
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = init(
        cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
    )

    # Process one day
    SoilTempArray_new, rSoilTempArrayRate_new, pSoilLayerDepth_new = process(
        cAVT=cAVT,
        cABD=cABD,
        pSoilLayerDepth=pSoilLayerDepth,
        cSoilLayerDepth=cSoilLayerDepth,
        SoilTempArray=SoilTempArray,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        iSoilWaterContent=iSoilWaterContent,
        iDoInitialize=iDoInitialize,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth
    )

    # Expected output (after process)
    expected_SoilTempArray = [
        13.624360856350041,
        13.399968504634286,
        12.599999999999845,
        12.2,
        11.4,
        10.6,
        9.799999999999999,
        9.0,
    ]

    assert len(SoilTempArray_new) == len(expected_SoilTempArray)
    for a, b in zip(SoilTempArray_new, expected_SoilTempArray):
        assert abs(a - b) < 1e-9, f"Mismatch: {a} vs {b}"