from typing import List
import math

def Init(LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers to default value.

    Inputs:
    - LayerThickness: List[float] (m)

    Returns:
    - SoilTemperatureByLayers: List[float] (degC)
    """
    SoilTemperatureByLayers = [15.0 for _ in range(len(LayerThickness))]
    return SoilTemperatureByLayers


def Estimate(
    SoilProfileDepth: float,
    SurfaceSoilTemperature: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: List[float],
    VolumetricWaterContent: List[float],
    LayerThickness: List[float],
    SoilTemperatureByLayers: List[float],
    LagCoefficient: float,
) -> List[float]:
    """
    Compute soil temperature by layers using the SWAT method (daily time step).

    Inputs:
    - SoilProfileDepth: float (m)
    - SurfaceSoilTemperature: float (degC)
    - AirTemperatureAnnualAverage: float (degC)
    - BulkDensity: List[float] (t m-3), per layer
    - VolumetricWaterContent: List[float] (m3 m-3), per layer
    - LayerThickness: List[float] (m), per layer
    - SoilTemperatureByLayers: List[float] (degC), previous day's layer temperatures
    - LagCoefficient: float (dimensionless, 0-1)

    Returns:
    - SoilTemperatureByLayers: List[float] (degC), updated layer temperatures
    """
    # Create a copy to avoid mutating the input list
    SoilTemperatureByLayers = SoilTemperatureByLayers.copy()

    # Conversion to mm
    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(LayerThickness)):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0]) / (
        BulkDensity[0] + 686.0 * math.exp(-5.63 * BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature)
            + SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, len(LayerThickness)):
        _DepthBottom = _DepthBottom + LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i]) / (
            BulkDensity[i] + 686.0 * math.exp(-5.63 * BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth))
            * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature)
                + SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers