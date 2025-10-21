from typing import List
import math

def Init(LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - LayerThickness: List[float] (m), thickness of each soil layer.

    Returns:
    - SoilTemperatureByLayers: List[float] (degC), initialized to 15 degC for each layer.
    """
    return [15.0 for _ in LayerThickness]


def Estimate(
    LagCoefficient: float,
    SoilProfileDepth: float,
    SurfaceSoilTemperature: float,
    AirTemperatureAnnualAverage: float,
    SoilTemperatureByLayers: List[float],
    BulkDensity: List[float],
    VolumetricWaterContent: List[float],
    LayerThickness: List[float],
) -> List[float]:
    """
    Update soil temperature by layers using SWAT method.

    Inputs:
    - LagCoefficient: float (dimensionless), default behavior in original model was 0.8 if unset.
    - SoilProfileDepth: float (m)
    - SurfaceSoilTemperature: float (degC)
    - AirTemperatureAnnualAverage: float (degC)
    - SoilTemperatureByLayers: List[float] (degC), previous day's soil temperature by layer
    - BulkDensity: List[float] (t m-3)
    - VolumetricWaterContent: List[float] (m3 m-3)
    - LayerThickness: List[float] (m)

    Returns:
    - SoilTemperatureByLayers: List[float] (degC), updated temperatures by layer.
    """
    n_layers = len(LayerThickness)
    if n_layers == 0:
        return []

    # Conversion to mm
    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    SoilTemperatureByLayers_out = [0.0] * n_layers

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0]) / (
        BulkDensity[0] + 686.0 * math.exp(-5.63 * BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / (
        (0.356 - 0.144 * BulkDensity[0]) * _SoilProfileDepthmm
    )
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers_out[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature)
            + SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n_layers):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i]) / (
            BulkDensity[i] + 686.0 * math.exp(-5.63 * BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.356 - 0.144 * BulkDensity[i]) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth))
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers_out[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature)
                + SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers_out