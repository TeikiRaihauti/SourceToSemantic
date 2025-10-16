import math

def initialize_STMPsimCalculator(cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT):
    """
    Initialize soil temperature state and depth arrays for the STMP model.

    Parameters
    - cSoilLayerDepth: list of floats, depths to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: float, mean air temperature on the first day (°C)
    - cDampingDepth: float, initial damping depth of soil (m)
    - cAVT: float, constant temperature of deepest soil layer - long term mean air temperature (°C)

    Returns
    - SoilTempArray: list of floats, initial soil temperature per (possibly extended) layer (°C)
    - rSoilTempArrayRate: list of floats, initial daily temperature change per layer (°C/day), initialized to 0.0
    - pSoilLayerDepth: list of floats, depths to the bottom of each (possibly extended) soil layer (m)
    """
    Z = list(cSoilLayerDepth)
    tProfileDepth = Z[-1]
    firstDayMeanTemp = cFirstDayMeanTemp
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
        # set linear approximation to the soil temperature as initial value
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth
    return tStmp, tStmpRate, tz


def model_STMPsimCalculator(
    cSoilLayerDepth,
    cABD,
    cAVT,
    iSoilWaterContent,
    iSoilSurfaceTemperature,
    iDoInitialize,
    cFirstDayMeanTemp,
    cDampingDepth,
    SoilTempArray,
    rSoilTempArrayRate,
    pSoilLayerDepth,
):
    """
    Main daily soil temperature process for the STMP model.

    Parameters
    - cSoilLayerDepth: list of floats, depths to the bottom of each soil layer (m)
    - cABD: float, mean bulk density (t/m3)
    - cAVT: float, constant temperature of deepest soil layer - long term mean air temperature (°C)
    - iSoilWaterContent: float, water content, sum over whole soil profile (mm)
    - iSoilSurfaceTemperature: float, soil surface temperature (°C)
    - iDoInitialize: bool, if True, reinitialize the model state before processing the day
    - cFirstDayMeanTemp: float, mean air temperature on first day (°C) - used only when iDoInitialize is True
    - cDampingDepth: float, initial damping depth of soil (m) - used only when iDoInitialize is True
    - SoilTempArray: list of floats, current soil temperature per layer (°C)
    - rSoilTempArrayRate: list of floats, current daily temperature change per layer (°C/day)
    - pSoilLayerDepth: list of floats, depths to the bottom of each soil layer, extended if needed (m)

    Returns
    - SoilTempArray: list of floats, updated soil temperature per layer (°C)
    - rSoilTempArrayRate: list of floats, updated daily temperature change per layer (°C/day)
    - pSoilLayerDepth: list of floats, (re)initialized or unchanged depths per layer (m)
    """
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = initialize_STMPsimCalculator(
            cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
        )

    tZp = list(pSoilLayerDepth)
    tZc = list(cSoilLayerDepth)
    XLAG = 0.8  # LAG = Coefficient for weighting yesterday's soil temperature
    XLG1 = 1.0 - XLAG
    DP = (1.0 + 2.5 * cABD / (cABD + math.exp(6.53 - 5.63 * cABD)))  # Maximum damping depth (m)
    WC = 0.001 * iSoilWaterContent / ((0.356 - 0.144 * cABD) * tZc[-1])
    DD = math.exp(math.log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP  # Damping depth (m), factor 2 in exponent

    Z1 = 0.0  # Depth of the bottom of the previous soil layer (m), initialized with 0
    new_rates = [0.0] * len(SoilTempArray)
    new_temps = list(SoilTempArray)
    for i in range(len(SoilTempArray)):
        ZD = 0.5 * (Z1 + tZp[i]) / DD  # Mid-layer depth divided by damping depth
        FZ_component = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD)) * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (FZ_component + iSoilSurfaceTemperature - new_temps[i])  # daily change (°C)
        Z1 = tZp[i]
        new_rates[i] = RATE
        new_temps[i] = new_temps[i] + RATE

    return new_temps, new_rates, tZp


def test_STMPsimCalculator_case0():
    """
    Derived from STMPsimCalculator.fillTestVariables, test index 0.
    Initializes the model and performs one process step, then checks expected SoilTempArray.
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
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = initialize_STMPsimCalculator(
        cSoilLayerDepth, cFirstDayMeanTemp, cDampingDepth, cAVT
    )

    # Process one timestep
    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = model_STMPsimCalculator(
        cSoilLayerDepth=cSoilLayerDepth,
        cABD=cABD,
        cAVT=cAVT,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        iDoInitialize=False,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
    )

    expected = [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    assert len(SoilTempArray) == len(expected), f"Length mismatch: {len(SoilTempArray)} vs {len(expected)}"
    for a, b in zip(SoilTempArray, expected):
        assert abs(a - b) < 1e-6, f"SoilTempArray mismatch: got {a}, expected {b}"