def Init(LayerThickness: list[float]) -> list[float]:
    """
    initialization of the SoilTemperatureSWAT component

    Inputs:
    - LayerThickness: list[float]
        Soil layer thickness (m)

    Returns:
    - SoilTemperatureByLayers: list[float]
        Initialized soil temperature of each layer (degC)
    """
    SoilTemperatureByLayers = [15.0 for _ in range(len(LayerThickness))]
    return SoilTemperatureByLayers


def CalculateModel(
    ex_VolumetricWaterContent: list[float],
    a_SurfaceSoilTemperature: float,
    s_SoilTemperatureByLayers: list[float],
    LayerThickness: list[float],
    LagCoefficient: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: list[float],
    SoilProfileDepth: float,
) -> list[float]:
    """
    Algorithm of the SoilTemperatureSWAT component

    Inputs:
    - ex_VolumetricWaterContent: list[float]
        Volumetric soil water content (m3 m-3)
    - a_SurfaceSoilTemperature: float
        Average surface soil temperature (degC)
    - s_SoilTemperatureByLayers: list[float]
        Previous day's soil temperature of each layer (degC)
    - LayerThickness: list[float]
        Soil layer thickness (m)
    - LagCoefficient: float
        Lag coefficient (dimensionless)
    - AirTemperatureAnnualAverage: float
        Annual average air temperature (degC)
    - BulkDensity: list[float]
        Bulk density (t m-3)
    - SoilProfileDepth: float
        Soil profile depth (m)

    Returns:
    - SoilTemperatureByLayers: list[float]
        Updated soil temperature of each layer (degC)
    """
    import math

    VolumetricWaterContent = ex_VolumetricWaterContent
    SurfaceSoilTemperature = a_SurfaceSoilTemperature
    SoilTemperatureByLayers = list(s_SoilTemperatureByLayers)

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

    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (
        2500.0 * BulkDensity[0] / (BulkDensity[0] + (686.0 * math.exp(-5.63 * BulkDensity[0])))
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * BulkDensity[0])) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        math.log(500.0 / _MaximumDumpingDepth)
        * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
    SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    for i in range(1, len(LayerThickness)):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + (LayerThickness[i] * 1000.0 / 2.0)
        _MaximumDumpingDepth = 1000.0 + (
            2500.0 * BulkDensity[i] / (BulkDensity[i] + (686.0 * math.exp(-5.63 * BulkDensity[i])))
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * BulkDensity[i])) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            math.log(500.0 / _MaximumDumpingDepth)
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
        SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return SoilTemperatureByLayers