from math import floor, ceil, fabs, exp, log
from typing import List, Tuple


def STMPsimCalculator_initialize(
    cSoilLayerDepth: List[float],
    cFirstDayMeanTemp: float,
    cDampingDepth: float,
    cAVT: float,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Initialization for STMP soil temperature component.

    Inputs:
    - cSoilLayerDepth: List[float] Depth to the bottom of each soil layer (m)
    - cFirstDayMeanTemp: float Mean air temperature on first day (°C)
    - cDampingDepth: float Damping depth of soil (m)
    - cAVT: float Long-term average annual air temperature (°C)

    Returns:
    - SoilTempArray: List[float] Initial soil temperatures per layer (°C)
    - rSoilTempArrayRate: List[float] Initial daily soil temperature change rates (°C/day), zeros
    - pSoilLayerDepth: List[float] Depth to the bottom of each (possibly extended) layer (m)
    """
    Z = [float(z) for z in cSoilLayerDepth]
    tProfileDepth = Z[-1]
    firstDayMeanTemp = float(cFirstDayMeanTemp)
    additionalDepth = float(cDampingDepth) - tProfileDepth
    firstAdditionalLayerHight = additionalDepth - floor(additionalDepth)
    layers = int(fabs(ceil(additionalDepth))) + len(Z)

    tStmp: List[float] = [0.0] * layers
    tStmpRate: List[float] = [0.0] * layers
    tz: List[float] = [0.0] * layers

    for i in range(layers):
        if i < len(Z):
            depth = Z[i]
        else:
            # first additional layer may be smaller than 1 m; subsequent in 1 m steps
            depth = tProfileDepth + firstAdditionalLayerHight + (i - len(Z))
        tz[i] = depth
        # Linear approximation from surface (firstDayMeanTemp) to damping depth (cAVT)
        tStmp[i] = (firstDayMeanTemp * (cDampingDepth - depth) + cAVT * depth) / cDampingDepth

    return tStmp, tStmpRate, tz


def _compute_damping_depth(
    cABD: float,
    iSoilWaterContent: float,
    cSoilLayerDepth: List[float],
) -> float:
    """
    Supporting function: compute damping depth DD (m).

    Inputs:
    - cABD: float Mean bulk density (t m-3)
    - iSoilWaterContent: float Water content, sum of whole soil profile (mm)
    - cSoilLayerDepth: List[float] Depth to the bottom of each soil layer (m)

    Returns:
    - DD: float Damping depth (m)
    """
    ABD = float(cABD)
    ZN = float(cSoilLayerDepth[-1])
    DP = 1.0 + 2.5 * ABD / (ABD + exp(6.53 - 5.63 * ABD))
    WC = 0.001 * float(iSoilWaterContent) / ((0.356 - 0.144 * ABD) * ZN)
    DD = exp(log(0.5 / DP) * ((1.0 - WC) / (1.0 + WC)) * 2.0) * DP
    return DD


def STMPsimCalculator_process(
    SoilTempArray: List[float],
    rSoilTempArrayRate: List[float],
    pSoilLayerDepth: List[float],
    cSoilLayerDepth: List[float],
    cABD: float,
    iSoilWaterContent: float,
    iSoilSurfaceTemperature: float,
    cAVT: float,
    iDoInitialize: bool = False,
    cFirstDayMeanTemp: float = 0.0,
    cDampingDepth: float = 0.0,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Main daily process for STMP soil temperature component.

    Inputs:
    - SoilTempArray: List[float] Current soil temperatures per layer (°C)
    - rSoilTempArrayRate: List[float] Current daily soil temperature change rates (°C/day)
    - pSoilLayerDepth: List[float] Depth to the bottom of each (possibly extended) layer (m)
    - cSoilLayerDepth: List[float] Depth to the bottom of original soil layers (m)
    - cABD: float Mean bulk density (t m-3)
    - iSoilWaterContent: float Water content of whole soil profile (mm)
    - iSoilSurfaceTemperature: float Soil surface temperature (°C)
    - cAVT: float Long-term average annual air temperature (°C)
    - iDoInitialize: bool If True, re-initialize state using cFirstDayMeanTemp and cDampingDepth before processing
    - cFirstDayMeanTemp: float Required if iDoInitialize is True; mean air temperature on first day (°C)
    - cDampingDepth: float Required if iDoInitialize is True; damping depth (m)

    Returns:
    - SoilTempArray: List[float] Updated soil temperatures per layer (°C)
    - rSoilTempArrayRate: List[float] Updated daily soil temperature change rates (°C/day)
    - pSoilLayerDepth: List[float] Potentially updated layer depths if re-initialized (m)
    """
    # Handle re-initialization request
    if iDoInitialize:
        SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_initialize(
            cSoilLayerDepth=cSoilLayerDepth,
            cFirstDayMeanTemp=cFirstDayMeanTemp,
            cDampingDepth=cDampingDepth,
            cAVT=cAVT,
        )

    # Ensure local copies to avoid mutating inputs in place
    temps = list(SoilTempArray)
    rates = [0.0] * len(temps)
    depths_p = list(pSoilLayerDepth)
    depths_c = list(cSoilLayerDepth)

    # Parameters and damping depth
    XLAG = 0.8
    XLG1 = 1.0 - XLAG
    DD = _compute_damping_depth(cABD=cABD, iSoilWaterContent=iSoilWaterContent, cSoilLayerDepth=depths_c)

    Z1 = 0.0
    for i in range(len(temps)):
        ZD = 0.5 * (Z1 + depths_p[i]) / DD
        FZ = ZD / (ZD + exp(-0.8669 - 2.0775 * ZD))
        RATE = FZ * (cAVT - iSoilSurfaceTemperature)
        RATE = XLG1 * (RATE + iSoilSurfaceTemperature - temps[i])
        Z1 = depths_p[i]
        rates[i] = RATE
        temps[i] = temps[i] + RATE

    return temps, rates, depths_p


def test_STMPsimCalculator_case0() -> None:
    """
    Derived test from STMPsimCalculator.fillTestVariables for a single daily step.
    Verifies SoilTempArray after one process step following initialization.

    Inputs (from test):
    - cSoilLayerDepth = [0.1, 0.5, 1.5] (m)
    - cFirstDayMeanTemp = 15 (°C)
    - cAVT = 9 (°C)
    - cABD = 1.4 (t m-3)
    - cDampingDepth = 6 (m)
    - iSoilWaterContent = 0.3 (mm)
    - iSoilSurfaceTemperature = 6 (°C)

    Expected SoilTempArray after process:
    [13.624360856350041, 13.399968504634286, 12.599999999999845, 12.2, 11.4, 10.6, 9.799999999999999, 9.0]
    """
    cSoilLayerDepth = [0.1, 0.5, 1.5]
    cFirstDayMeanTemp = 15.0
    cAVT = 9.0
    cABD = 1.4
    cDampingDepth = 6.0
    iSoilWaterContent = 0.3
    iSoilSurfaceTemperature = 6.0

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_initialize(
        cSoilLayerDepth=cSoilLayerDepth,
        cFirstDayMeanTemp=cFirstDayMeanTemp,
        cDampingDepth=cDampingDepth,
        cAVT=cAVT,
    )

    SoilTempArray, rSoilTempArrayRate, pSoilLayerDepth = STMPsimCalculator_process(
        SoilTempArray=SoilTempArray,
        rSoilTempArrayRate=rSoilTempArrayRate,
        pSoilLayerDepth=pSoilLayerDepth,
        cSoilLayerDepth=cSoilLayerDepth,
        cABD=cABD,
        iSoilWaterContent=iSoilWaterContent,
        iSoilSurfaceTemperature=iSoilSurfaceTemperature,
        cAVT=cAVT,
        iDoInitialize=False,
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
    for a, e in zip(SoilTempArray, expected):
        assert abs(a - e) < 1e-9, f"{a} != {e}"