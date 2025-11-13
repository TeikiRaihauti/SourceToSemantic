def init_soiltemperatureswat(
    VolumetricWaterContent: list[float],
    LayerThickness: list[float],
    LagCoefficient: float,
    AirTemperatureAnnualAverage: float,
    BulkDensity: list[float],
    SoilProfileDepth: float,
) -> list[float]:
    """
    Initialize soil temperature by layers.

    Inputs:
    - VolumetricWaterContent: list[float], m3 m-3
    - LayerThickness: list[float], m
    - LagCoefficient: float, dimensionless
    - AirTemperatureAnnualAverage: float, degC
    - BulkDensity: list[float], t m-3
    - SoilProfileDepth: float, m

    Returns:
    - SoilTemperatureByLayers: list[float], degC
    """
    # Set each layer temperature to 15 degC, same length as LayerThickness
    SoilTemperatureByLayers = [15.0 for _ in range(len(LayerThickness))]
    return SoilTemperatureByLayers


def model_soiltemperatureswat(
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
    Compute soil temperature by layers using SWAT method.

    Inputs:
    - VolumetricWaterContent: list[float], m3 m-3
    - SurfaceSoilTemperature: float, degC
    - LayerThickness: list[float], m
    - LagCoefficient: float, dimensionless
    - SoilTemperatureByLayers: list[float], degC (previous state)
    - AirTemperatureAnnualAverage: float, degC
    - BulkDensity: list[float], t m-3
    - SoilProfileDepth: float, m

    Returns:
    - SoilTemperatureByLayers: list[float], degC (updated state)
    """
    from math import exp, log

    n_layers = len(LayerThickness)
    updated_SoilTemperatureByLayers = SoilTemperatureByLayers[:]

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

    # First layer (index 0)
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (
        2500.0
        * BulkDensity[0]
        / (BulkDensity[0] + (686.0 * exp((-5.63) * BulkDensity[0])))
    )
    _ScalingFactor = _TotalWaterContentmm / (
        (0.356 - (0.144 * BulkDensity[0])) * _SoilProfileDepthmm
    )
    _DumpingDepth = _MaximumDumpingDepth * exp(
        log(500.0 / _MaximumDumpingDepth)
        * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + exp(-0.867 - (2.078 * _RatioCenter)))
    updated_SoilTemperatureByLayers[0] = LagCoefficient * SoilTemperatureByLayers[0] + (
        (1.0 - LagCoefficient)
        * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    # Remaining layers
    _DepthBottom = 0.0
    for i in range(1, n_layers):
        _DepthBottom += LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + (LayerThickness[i] * 1000.0 / 2.0)
        _MaximumDumpingDepth = 1000.0 + (
            2500.0
            * BulkDensity[i]
            / (BulkDensity[i] + (686.0 * exp((-5.63) * BulkDensity[i])))
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.356 - (0.144 * BulkDensity[i])) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * exp(
            log(500.0 / _MaximumDumpingDepth)
            * (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (
            _RatioCenter + exp(-0.867 - (2.078 * _RatioCenter))
        )
        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (
                (1.0 - LagCoefficient)
                * (
                    _DepthFactor
                    * (AirTemperatureAnnualAverage - SurfaceSoilTemperature)
                    + SurfaceSoilTemperature
                )
            )
        )

    return updated_SoilTemperatureByLayers