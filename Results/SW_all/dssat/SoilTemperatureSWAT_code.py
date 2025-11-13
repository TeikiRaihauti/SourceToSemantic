from typing import List


def init_soiltemperatureswat(
    VolumetricWaterContent: List[float],
    LayerThickness: List[float],
    LagCoefficient: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: List[float],
    SoilProfileDepth: float,
    SoilTemperatureByLayers: List[float],
) -> List[float]:
    """
    Initialization function for SoilTemperatureSWAT.

    Inputs:
    - VolumetricWaterContent: list[float] (m3 m-3)
    - LayerThickness: list[float] (m)
    - LagCoefficient: float (dimensionless)
    - AirTemperatureAnnualAverage: float (degC)
    - BulkDensity: list[float] (t m-3)
    - SoilProfileDepth: float (m)
    - SoilTemperatureByLayers: list[float] (degC) [ignored; output initialized]

    Returns:
    - SoilTemperatureByLayers: list[float] (degC)
    """
    # Initialize all soil layer temperatures to 15 degC
    SoilTemperatureByLayers_out = [15.0 for _ in LayerThickness]
    return SoilTemperatureByLayers_out


def model_soiltemperatureswat(
    VolumetricWaterContent: List[float],
    SurfaceSoilTemperature: float,
    LayerThickness: List[float],
    LagCoefficient: float,
    SoilTemperatureByLayers: List[float],
    AirTemperatureAnnualAverage: float,
    BulkDensity: List[float],
    SoilProfileDepth: float,
) -> List[float]:
    """
    Main process function for SoilTemperatureSWAT.

    Inputs:
    - VolumetricWaterContent: list[float] (m3 m-3)
    - SurfaceSoilTemperature: float (degC)
    - LayerThickness: list[float] (m)
    - LagCoefficient: float (dimensionless)
    - SoilTemperatureByLayers: list[float] (degC), previous state
    - AirTemperatureAnnualAverage: float (degC)
    - BulkDensity: list[float] (t m-3)
    - SoilProfileDepth: float (m)

    Returns:
    - SoilTemperatureByLayers: list[float] (degC), updated state
    """
    import math

    n_layers = len(LayerThickness)
    SoilTemperatureByLayers_out = [0.0] * n_layers

    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0
    if n_layers == 0:
        return SoilTemperatureByLayers_out

    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    bd0 = BulkDensity[0]
    _MaximumDumpingDepth = 1000.0 + (2500.0 * bd0 / (bd0 + (686.0 * math.exp((-5.63) * bd0))))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * bd0)) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        math.log(500.0 / _MaximumDumpingDepth) * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
    SoilTemperatureByLayers_out[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    for i in range(1, n_layers):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + (LayerThickness[i] * 1000.0 / 2.0)
        bdi = BulkDensity[i]
        _MaximumDumpingDepth = 1000.0 + (2500.0 * bdi / (bdi + (686.0 * math.exp((-5.63) * bdi))))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * bdi)) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            math.log(500.0 / _MaximumDumpingDepth) * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
        SoilTemperatureByLayers_out[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return SoilTemperatureByLayers_out