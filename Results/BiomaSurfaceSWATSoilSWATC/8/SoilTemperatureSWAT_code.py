from typing import List
import math

def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: List[float]
        Soil layer thickness for each layer [m]

    Returns:
    - SoilTemperatureByLayers: List[float]
        Initialized soil temperature of each layer [degC], set to 15 for all layers
    """
    n_layers = len(states_LayerThickness)
    SoilTemperatureByLayers = [15.0 for _ in range(n_layers)]
    return SoilTemperatureByLayers


def Estimate(
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: List[float],
    states_BulkDensity: List[float],
    states_VolumetricWaterContent: List[float],
    states_LayerThickness: List[float],
    LagCoefficient: float = 0.8,
) -> List[float]:
    """
    Update soil temperature by layers using SWAT method.

    Inputs:
    - states_SoilProfileDepth: float
        Soil profile depth [m]
    - states_SurfaceSoilTemperature: float
        Average surface soil temperature [degC]
    - exogenous_AirTemperatureAnnualAverage: float
        Annual average air temperature [degC]
    - states_SoilTemperatureByLayers: List[float]
        Previous day's soil temperature of each layer [degC]
    - states_BulkDensity: List[float]
        Bulk density for each layer [t m-3]
    - states_VolumetricWaterContent: List[float]
        Volumetric soil water content for each layer [m3 m-3]
    - states_LayerThickness: List[float]
        Soil layer thickness for each layer [m]
    - LagCoefficient: float (default=0.8)
        Lag coefficient that controls the influence of the previous day's temperature
        on the current day's temperature [-]

    Returns:
    - SoilTemperatureByLayers: List[float]
        Updated soil temperature of each layer [degC]
    """
    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content [mm]
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables initialization
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    SoilTemperatureByLayers_out = [0.0] * len(states_LayerThickness)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers_out[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
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
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers_out[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers_out