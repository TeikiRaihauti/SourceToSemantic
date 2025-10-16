import math
from typing import List


def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialization function for SoilTemperatureSWAT.
    Initializes states_SoilTemperatureByLayers to 15 degC for each soil layer.

    Inputs:
    - states_LayerThickness: list of soil layer thicknesses [m]

    Returns:
    - states_SoilTemperatureByLayers: list of initialized soil temperatures by layer [degC]
    """
    states_SoilTemperatureByLayers = [15.0 for _ in states_LayerThickness]
    return states_SoilTemperatureByLayers


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
    Main biophysical process function for SoilTemperatureSWAT.
    Computes soil temperature by layers using the SWAT method (Neitsch et al., 2000).

    Inputs:
    - LagCoefficient: dimensionless lag coefficient controlling the influence of the previous day's temperature [0-1]
    - states_SoilProfileDepth: soil profile depth [m]
    - states_SurfaceSoilTemperature: average surface soil temperature [degC]
    - exogenous_AirTemperatureAnnualAverage: annual average air temperature [degC]
    - states_SoilTemperatureByLayers: previous day's soil temperature by layer [degC]
    - states_BulkDensity: bulk density by layer [t m-3]
    - states_VolumetricWaterContent: volumetric soil water content by layer [m3 m-3]
    - states_LayerThickness: soil layer thickness by layer [m]

    Returns:
    - updated_states_SoilTemperatureByLayers: updated soil temperature by layer [degC]
    """
    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content over profile in mm
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    updated_states_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    updated_states_SoilTemperatureByLayers[0] = (
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
            (math.log(500.0 / _MaximumDumpingDepth))
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        updated_states_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return updated_states_SoilTemperatureByLayers