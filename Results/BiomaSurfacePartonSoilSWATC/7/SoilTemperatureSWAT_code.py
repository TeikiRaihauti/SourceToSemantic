from typing import List, Sequence
import math

def Init(LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialization function for SoilTemperatureSWAT.

    Inputs:
    - LayerThickness: List[float], m

    Returns:
    - SoilTemperatureByLayers: List[float], degC
    """
    n = len(LayerThickness)
    SoilTemperatureByLayers = [15.0 for _ in range(n)]
    return SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    SoilProfileDepth: float,
    SurfaceSoilTemperature: float,
    AirTemperatureAnnualAverage: float,
    SoilTemperatureByLayers: Sequence[float],
    BulkDensity: Sequence[float],
    VolumetricWaterContent: Sequence[float],
    LayerThickness: Sequence[float],
) -> List[float]:
    """
    Main biophysical process function for SoilTemperatureSWAT.

    Inputs:
    - LagCoefficient: float, dimensionless
    - SoilProfileDepth: float, m
    - SurfaceSoilTemperature: float, degC
    - AirTemperatureAnnualAverage: float, degC
    - SoilTemperatureByLayers: List[float], degC (previous day's values)
    - BulkDensity: List[float], t m-3
    - VolumetricWaterContent: List[float], m3 m-3
    - LayerThickness: List[float], m

    Returns:
    - SoilTemperatureByLayers: List[float], degC (updated values)
    """
    n = len(LayerThickness)
    if not (len(BulkDensity) == n and len(VolumetricWaterContent) == n and len(SoilTemperatureByLayers) == n):
        raise ValueError("All layer-based inputs must have the same length")

    # Conversion to mm
    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    new_SoilTemperatureByLayers = list(SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    bd0 = BulkDensity[0]
    _MaximumDumpingDepth = 1000.0 + (2500.0 * bd0) / (bd0 + 686.0 * math.exp(-5.63 * bd0))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * bd0) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    new_SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    # Other layers
    for i in range(1, n):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0
        bdi = BulkDensity[i]
        _MaximumDumpingDepth = 1000.0 + (2500.0 * bdi) / (bdi + 686.0 * math.exp(-5.63 * bdi))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * bdi) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        new_SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return new_SoilTemperatureByLayers