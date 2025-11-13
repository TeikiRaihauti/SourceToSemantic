from typing import List
import math

def Init(LayerThickness: List[float]) -> List[float]:
    """
    Initialize soil temperature by layers.

    Inputs:
    - LayerThickness: List[float] (m), thickness of each soil layer

    Returns:
    - SoilTemperatureByLayers: List[float] (degC), initialized temperature for each layer
    """
    SoilTemperatureByLayers: List[float] = [15.0 for _ in range(len(LayerThickness))]
    return SoilTemperatureByLayers


def CalculateModel(
    VolumetricWaterContent: List[float],
    SurfaceSoilTemperature: float,
    SoilTemperatureByLayers: List[float],
    LayerThickness: List[float],
    LagCoefficient: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: List[float],
    SoilProfileDepth: float
) -> List[float]:
    """
    Update soil temperature by layers using the SWAT method.

    Inputs:
    - VolumetricWaterContent: List[float] (m3 m-3), volumetric water content per soil layer
    - SurfaceSoilTemperature: float (degC), current surface soil temperature
    - SoilTemperatureByLayers: List[float] (degC), previous day's soil temperature per layer
    - LayerThickness: List[float] (m), thickness of each soil layer
    - LagCoefficient: float (dimensionless), daily temperature lag coefficient [0,1]
    - AirTemperatureAnnualAverage: float (degC), annual average air temperature
    - BulkDensity: List[float] (t m-3), bulk density per soil layer
    - SoilProfileDepth: float (m), total soil profile depth

    Returns:
    - SoilTemperatureByLayers: List[float] (degC), updated soil temperature per layer
    """
    # Copy previous day's temperatures to ensure purity
    prev_SoilTemperatureByLayers: List[float] = list(SoilTemperatureByLayers)
    new_SoilTemperatureByLayers: List[float] = list(SoilTemperatureByLayers)

    i: int
    _SoilProfileDepthmm: float
    _TotalWaterContentmm: float
    _MaximumDumpingDepth: float
    _DumpingDepth: float
    _ScalingFactor: float
    _DepthBottom: float
    _RatioCenter: float
    _DepthFactor: float
    _DepthCenterLayer: float

    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    _TotalWaterContentmm = 0.0
    for i in range(len(LayerThickness)):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0] / (BulkDensity[0] + (686.0 * math.exp(-5.630 * BulkDensity[0]))))
    _ScalingFactor = _TotalWaterContentmm / ((0.3560 - (0.1440 * BulkDensity[0])) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(math.log(500.0 / _MaximumDumpingDepth) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.8670 - (2.0780 * _RatioCenter)))
    new_SoilTemperatureByLayers[0] = (
        LagCoefficient * prev_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient) * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    # Remaining layers
    for i in range(1, len(LayerThickness)):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + (LayerThickness[i] * 1000.0 / 2.0)
        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i] / (BulkDensity[i] + (686.0 * math.exp(-5.630 * BulkDensity[i]))))
        _ScalingFactor = _TotalWaterContentmm / ((0.3560 - (0.1440 * BulkDensity[i])) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(math.log(500.0 / _MaximumDumpingDepth) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.8670 - (2.0780 * _RatioCenter)))
        new_SoilTemperatureByLayers[i] = (
            LagCoefficient * prev_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient) * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return new_SoilTemperatureByLayers