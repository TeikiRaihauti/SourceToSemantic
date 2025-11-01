from typing import List, Sequence
import math

def Init(states_LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialize state variables for the SoilTemperatureSWAT model.

    Inputs:
    - states_LayerThickness: Sequence[float]
      Soil layer thickness for each layer [m]

    Returns:
    - states_SoilTemperatureByLayers: List[float]
      Initialized soil temperature for each layer [degC]
    """
    states_SoilTemperatureByLayers: List[float] = [15.0 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: Sequence[float],
    states_BulkDensity: Sequence[float],
    states_VolumetricWaterContent: Sequence[float],
    states_LayerThickness: Sequence[float],
) -> List[float]:
    """
    Compute soil temperature by layers using the SWAT approach (Neitsch et al., 2000).

    Inputs:
    - LagCoefficient: float
      Lag coefficient controlling influence of previous day's temperature [-] (default in original model: 0.8)
    - states_SoilProfileDepth: float
      Soil profile depth [m]
    - states_SurfaceSoilTemperature: float
      Average surface soil temperature [degC]
    - exogenous_AirTemperatureAnnualAverage: float
      Annual average air temperature [degC]
    - states_SoilTemperatureByLayers: Sequence[float]
      Previous day's soil temperature of each layer [degC]
    - states_BulkDensity: Sequence[float]
      Bulk density for each layer [t m-3]
    - states_VolumetricWaterContent: Sequence[float]
      Volumetric soil water content for each layer [m3 m-3]
    - states_LayerThickness: Sequence[float]
      Soil layer thickness for each layer [m]

    Returns:
    - states_SoilTemperatureByLayers: List[float]
      Updated soil temperature of each layer [degC]
    """
    # Prepare output as a copy to keep function pure
    states_SoilTemperatureByLayers_out: List[float] = list(states_SoilTemperatureByLayers)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content [mm]
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

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
        (math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    states_SoilTemperatureByLayers_out[0] = (
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

        states_SoilTemperatureByLayers_out[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return states_SoilTemperatureByLayers_out