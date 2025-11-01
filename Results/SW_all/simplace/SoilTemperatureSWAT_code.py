def init(LayerThickness: list[float]) -> list[float]:
    """
    Initialize soil temperature by layers for the SWAT soil temperature module.

    Parameters:
    - LayerThickness: list[float]
        Soil layer thicknesses [m] for each layer.

    Returns:
    - SoilTemperatureByLayers: list[float]
        Initialized soil temperature [degC] for each layer (set to 15 degC).
    """
    t_LayerThickness = LayerThickness
    t_SoilTemperatureByLayers = [0.0] * len(t_LayerThickness)
    for i in range(len(t_LayerThickness)):
        t_SoilTemperatureByLayers[i] = 15.0
    return t_SoilTemperatureByLayers


def process(
    VolumetricWaterContent: list[float],
    SurfaceSoilTemperature: float,
    LayerThickness: list[float],
    LagCoefficient: float,
    SoilTemperatureByLayers: list[float],
    AirTemperatureAnnualAverage: float,
    BulkDensity: list[float],
    SoilProfileDepth: float,
) -> list[float]:
    """
    Update soil temperature by layers using the SWAT approach for one timestep.

    Parameters:
    - VolumetricWaterContent: list[float]
        Volumetric soil water content [m3 m-3] per layer.
    - SurfaceSoilTemperature: float
        Surface soil temperature [degC].
    - LayerThickness: list[float]
        Soil layer thickness [m] per layer.
    - LagCoefficient: float
        Lag coefficient controlling influence of previous day (0-1).
    - SoilTemperatureByLayers: list[float]
        Current soil temperature [degC] per layer (state to be updated).
    - AirTemperatureAnnualAverage: float
        Annual average air temperature [degC].
    - BulkDensity: list[float]
        Bulk density [t m-3] per layer.
    - SoilProfileDepth: float
        Total soil profile depth [m].

    Returns:
    - SoilTemperatureByLayers: list[float]
        Updated soil temperature [degC] per layer.
    """
    import math

    t_VolumetricWaterContent = VolumetricWaterContent
    t_SurfaceSoilTemperature = SurfaceSoilTemperature
    t_LayerThickness = LayerThickness
    t_LagCoefficient = LagCoefficient
    t_SoilTemperatureByLayers = list(SoilTemperatureByLayers)
    t_AirTemperatureAnnualAverage = AirTemperatureAnnualAverage
    t_BulkDensity = BulkDensity
    t_SoilProfileDepth = SoilProfileDepth

    _SoilProfileDepthmm = t_SoilProfileDepth * 1000.0
    _TotalWaterContentmm = 0.0
    for i in range(len(t_LayerThickness)):
        _TotalWaterContentmm = _TotalWaterContentmm + (t_VolumetricWaterContent[i] * t_LayerThickness[i])
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0
    if len(t_LayerThickness) == 0:
        return t_SoilTemperatureByLayers

    _DepthCenterLayer = t_LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * t_BulkDensity[0] / (t_BulkDensity[0] + (686.0 * math.exp(-5.63 * t_BulkDensity[0]))))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * t_BulkDensity[0])) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(math.log(500.0 / _MaximumDumpingDepth) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
    t_SoilTemperatureByLayers[0] = (
        t_LagCoefficient * t_SoilTemperatureByLayers[0]
        + ((1.0 - t_LagCoefficient) * (_DepthFactor * (t_AirTemperatureAnnualAverage - t_SurfaceSoilTemperature) + t_SurfaceSoilTemperature))
    )

    for i in range(1, len(t_LayerThickness)):
        _DepthBottom = _DepthBottom + (t_LayerThickness[i - 1] * 1000.0)
        _DepthCenterLayer = _DepthBottom + (t_LayerThickness[i] * 1000.0 / 2.0)
        _MaximumDumpingDepth = 1000.0 + (2500.0 * t_BulkDensity[i] / (t_BulkDensity[i] + (686.0 * math.exp(-5.63 * t_BulkDensity[i]))))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - (0.144 * t_BulkDensity[i])) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(math.log(500.0 / _MaximumDumpingDepth) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - (2.078 * _RatioCenter)))
        t_SoilTemperatureByLayers[i] = (
            t_LagCoefficient * t_SoilTemperatureByLayers[i]
            + ((1.0 - t_LagCoefficient) * (_DepthFactor * (t_AirTemperatureAnnualAverage - t_SurfaceSoilTemperature) + t_SurfaceSoilTemperature))
        )

    return t_SoilTemperatureByLayers