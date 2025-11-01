from typing import List
import math


def Init(LayerThickness: List[float]) -> List[float]:
    """
    Initialize state variables for the SoilTemperatureSWAT strategy.

    Inputs:
    - LayerThickness: List[float]
        Soil layer thickness for each layer [m].

    Returns:
    - SoilTemperatureByLayers: List[float]
        Initial soil temperature for each layer [degC], set to 15.0.
    """
    n = len(LayerThickness)
    SoilTemperatureByLayers = [15.0 for _ in range(n)]
    return SoilTemperatureByLayers


def CalculateModel(
    ex_VolumetricWaterContent: List[float],
    a_SurfaceSoilTemperature: float,
    s_SoilTemperatureByLayers: List[float],
    LayerThickness: List[float],
    LagCoefficient: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: List[float],
    SoilProfileDepth: float,
) -> List[float]:
    """
    Compute soil temperature by layers using the SWAT method.

    Inputs:
    - ex_VolumetricWaterContent: List[float]
        Volumetric soil water content by layer [m3 m-3].
    - a_SurfaceSoilTemperature: float
        Average surface soil temperature [degC].
    - s_SoilTemperatureByLayers: List[float]
        Previous day's soil temperature by layer [degC].
    - LayerThickness: List[float]
        Soil layer thickness by layer [m].
    - LagCoefficient: float
        Lag coefficient controlling influence of previous day's temperature [0-1].
    - AirTemperatureAnnualAverage: float
        Annual average air temperature [degC].
    - BulkDensity: List[float]
        Bulk density by layer [t m-3].
    - SoilProfileDepth: float
        Total soil profile depth [m].

    Returns:
    - SoilTemperatureByLayers: List[float]
        Updated soil temperature by layer [degC].
    """
    VolumetricWaterContent = ex_VolumetricWaterContent
    SurfaceSoilTemperature = a_SurfaceSoilTemperature
    prev_SoilTemperatureByLayers = s_SoilTemperatureByLayers

    n = len(LayerThickness)
    SoilTemperatureByLayers = list(prev_SoilTemperatureByLayers)

    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    _TotalWaterContentmm = 0.0
    for i in range(n):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # Layer 0
    if n > 0:
        _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (
            2500.0
            * BulkDensity[0]
            / (BulkDensity[0] + (686.0 * math.exp(-5.630 * BulkDensity[0])))
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.3560 - (0.1440 * BulkDensity[0])) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            math.log(500.0 / _MaximumDumpingDepth)
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.8670 - (2.0780 * _RatioCenter)))
        SoilTemperatureByLayers[0] = (
            LagCoefficient * prev_SoilTemperatureByLayers[0]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    # Layers 1..n-1
    for i in range(1, n):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + (LayerThickness[i] * 1000.0 / 2.0)
        _MaximumDumpingDepth = 1000.0 + (
            2500.0
            * BulkDensity[i]
            / (BulkDensity[i] + (686.0 * math.exp(-5.630 * BulkDensity[i])))
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.3560 - (0.1440 * BulkDensity[i])) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            math.log(500.0 / _MaximumDumpingDepth)
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.8670 - (2.0780 * _RatioCenter)))
        SoilTemperatureByLayers[i] = (
            LagCoefficient * prev_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return SoilTemperatureByLayers