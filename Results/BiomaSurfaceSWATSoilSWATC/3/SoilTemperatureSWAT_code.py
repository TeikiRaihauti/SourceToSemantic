from typing import List, Sequence
import math


def Init(states_LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialize soil temperature profile by layers.

    Inputs:
    - states_LayerThickness: Sequence[float]
        Soil layer thickness for each layer [m].

    Returns:
    - states_SoilTemperatureByLayers: List[float]
        Initialized soil temperature for each layer [degC], set to 15 for all layers.
    """
    n_layers = len(states_LayerThickness)
    states_SoilTemperatureByLayers = [15.0] * n_layers
    return states_SoilTemperatureByLayers


def Estimate(
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: Sequence[float],
    states_BulkDensity: Sequence[float],
    states_VolumetricWaterContent: Sequence[float],
    states_LayerThickness: Sequence[float],
    LagCoefficient: float = 0.8,
) -> List[float]:
    """
    Soil temperature by layers using SWAT method (Neitsch et al., 2000).

    Inputs:
    - states_SoilProfileDepth: float
        Soil profile depth [m].
    - states_SurfaceSoilTemperature: float
        Average surface soil temperature [degC].
    - exogenous_AirTemperatureAnnualAverage: float
        Annual average air temperature [degC].
    - states_SoilTemperatureByLayers: Sequence[float]
        Previous day's soil temperature of each layer [degC].
    - states_BulkDensity: Sequence[float]
        Bulk density for each layer [t m-3].
    - states_VolumetricWaterContent: Sequence[float]
        Volumetric soil water content for each layer [m3 m-3].
    - states_LayerThickness: Sequence[float]
        Soil layer thickness for each layer [m].
    - LagCoefficient: float (default 0.8)
        Lag coefficient controlling influence of previous day's temperature [dimensionless].

    Returns:
    - updated_states_SoilTemperatureByLayers: List[float]
        Updated soil temperature for each layer [degC].
    """
    # Convert profile depth to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content [mm]
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    updated_states_SoilTemperatureByLayers: List[float] = list(states_SoilTemperatureByLayers)

    _DepthBottom = 0.0  # mm

    # First layer center depth [mm]
    if len(states_LayerThickness) > 0:
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

        updated_states_SoilTemperatureByLayers[0] = (
            LagCoefficient * states_SoilTemperatureByLayers[0]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (
                    exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature
                )
                + states_SurfaceSoilTemperature
            )
        )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom += states_LayerThickness[i - 1] * 1000.0
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

        updated_states_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (
                    exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature
                )
                + states_SurfaceSoilTemperature
            )
        )

    return updated_states_SoilTemperatureByLayers