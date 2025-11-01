from typing import List
import math

def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers array.

    Inputs:
    - states_LayerThickness: List[float], units: m

    Returns:
    - SoilTemperatureByLayers: List[float], initialized to 15 degC for each layer
    """
    SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
    return SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: List[float],
    states_BulkDensity: List[float],
    states_VolumetricWaterContent: List[float],
    states_LayerThickness: List[float],
) -> List[float]:
    """
    SoilTemperatureSWAT main process function.

    Inputs:
    - LagCoefficient: float (default 0.8), dimensionless
    - states_SoilProfileDepth: float, units: m
    - states_SurfaceSoilTemperature: float, units: degC
    - exogenous_AirTemperatureAnnualAverage: float, units: degC
    - states_SoilTemperatureByLayers: List[float], units: degC (previous day's)
    - states_BulkDensity: List[float], units: t m-3
    - states_VolumetricWaterContent: List[float], units: m3 m-3
    - states_LayerThickness: List[float], units: m

    Returns:
    - SoilTemperatureByLayers: List[float], updated soil temperature by layers (degC)
    """
    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    # internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / (
        (0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm
    )
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (
                exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature
            )
            + states_SurfaceSoilTemperature
        )
    )

    # other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (
            states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth))
            * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (
                    exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature
                )
                + states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers