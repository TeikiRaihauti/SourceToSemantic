import math
from typing import List, Sequence


def Init(states_LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: sequence of soil layer thicknesses (m)

    Returns:
    - SoilTemperatureByLayers: list of initial soil temperatures per layer (degC), all set to 15
    """
    return [15.0] * len(states_LayerThickness)


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
    Compute soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: lag coefficient (dimensionless)
    - states_SoilProfileDepth: soil profile depth (m)
    - states_SurfaceSoilTemperature: average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: annual average air temperature (degC)
    - states_SoilTemperatureByLayers: previous day's soil temperature of each layer (degC)
    - states_BulkDensity: bulk density per layer (t m-3)
    - states_VolumetricWaterContent: volumetric soil water content per layer (m3 m-3)
    - states_LayerThickness: soil layer thickness per layer (m)

    Returns:
    - SoilTemperatureByLayers: updated soil temperature of each layer (degC)
    """
    # Convert soil profile depth to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _DepthBottom = 0.0

    n = len(states_LayerThickness)
    SoilTemperatureByLayers = [0.0] * n

    if n == 0:
        return SoilTemperatureByLayers

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
        * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n):
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
            * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers