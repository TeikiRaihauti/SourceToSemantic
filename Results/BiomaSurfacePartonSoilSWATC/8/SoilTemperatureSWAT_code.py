from typing import List
import math

def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: List[float]
      Unit: m
      Description: Soil layer thickness for each layer.

    Returns:
    - states_SoilTemperatureByLayers: List[float]
      Unit: degC
      Description: Initialized soil temperature of each layer (set to 15 degC).
    """
    n_layers = 0 if states_LayerThickness is None else len(states_LayerThickness)
    return [15.0 for _ in range(n_layers)]


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
    SoilTemperatureSWAT main process: update SoilTemperatureByLayers using SWAT method.

    Inputs:
    - LagCoefficient: float
      Unit: dimensionless
      Description: Lag coefficient that controls influence of previous day's temperature.
    - states_SoilProfileDepth: float
      Unit: m
      Description: Soil profile depth.
    - states_SurfaceSoilTemperature: float
      Unit: degC
      Description: Average surface soil temperature.
    - exogenous_AirTemperatureAnnualAverage: float
      Unit: degC
      Description: Annual average air temperature.
    - states_SoilTemperatureByLayers: List[float]
      Unit: degC
      Description: Previous day's soil temperature by layer.
    - states_BulkDensity: List[float]
      Unit: t m-3
      Description: Bulk density by layer.
    - states_VolumetricWaterContent: List[float]
      Unit: m3 m-3
      Description: Volumetric soil water content by layer.
    - states_LayerThickness: List[float]
      Unit: m
      Description: Soil layer thickness by layer.

    Returns:
    - states_SoilTemperatureByLayers: List[float]
      Unit: degC
      Description: Updated soil temperature by layer.
    """
    n_layers = len(states_LayerThickness)
    if n_layers == 0:
        return []

    # Create a copy to avoid mutating the input
    updated_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content in mm (sum over layers of theta * thickness[m] * 1000)
    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer center depth [mm]
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

    updated_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n_layers):
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

        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return updated_SoilTemperatureByLayers