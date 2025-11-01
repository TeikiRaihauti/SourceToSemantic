from typing import List
import math


def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize soil temperatures by layer.

    Inputs:
    - states_LayerThickness: List[float], soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: List[float], initialized soil temperature of each layer (degC)
    """
    states_SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


def Estimate(
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: List[float],
    states_BulkDensity: List[float],
    states_VolumetricWaterContent: List[float],
    states_LayerThickness: List[float],
    LagCoefficient: float,
) -> List[float]:
    """
    Update soil temperature by layers using the SWAT method.

    Inputs:
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: List[float], previous day's soil temperature per layer (degC)
    - states_BulkDensity: List[float], bulk density per layer (t m-3)
    - states_VolumetricWaterContent: List[float], volumetric soil water content per layer (m3 m-3)
    - states_LayerThickness: List[float], layer thicknesses (m)
    - LagCoefficient: float, lag coefficient (dimensionless)

    Returns:
    - states_SoilTemperatureByLayers: List[float], updated soil temperature per layer (degC)
    """
    new_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    if len(states_LayerThickness) == 0:
        return new_SoilTemperatureByLayers

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
        * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    new_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
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
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        new_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return new_SoilTemperatureByLayers