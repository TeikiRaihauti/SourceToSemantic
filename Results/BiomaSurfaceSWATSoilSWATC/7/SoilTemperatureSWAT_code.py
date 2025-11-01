def Init(states_LayerThickness: list) -> list:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: list[float], soil layer thickness (m) for each layer

    Returns:
    - states_SoilTemperatureByLayers: list[float], initialized soil temperature of each layer (degC)
    """
    n = len(states_LayerThickness) if states_LayerThickness is not None else 0
    states_SoilTemperatureByLayers = [15.0] * n
    return states_SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: list,
    states_BulkDensity: list,
    states_VolumetricWaterContent: list,
    states_LayerThickness: list,
) -> list:
    """
    Update SoilTemperatureByLayers using SWAT method.

    Inputs:
    - LagCoefficient: float, lag coefficient (dimensionless)
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list[float], previous soil temperature by layer (degC)
    - states_BulkDensity: list[float], bulk density by layer (t m-3)
    - states_VolumetricWaterContent: list[float], volumetric water content by layer (m3 m-3)
    - states_LayerThickness: list[float], soil layer thickness by layer (m)

    Returns:
    - states_SoilTemperatureByLayers: list[float], updated soil temperature by layer (degC)
    """
    import math

    n_layers = len(states_LayerThickness) if states_LayerThickness is not None else 0
    if n_layers == 0:
        return list(states_SoilTemperatureByLayers) if states_SoilTemperatureByLayers is not None else []

    # Copy to avoid mutating input
    out_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    bd0 = states_BulkDensity[0]
    _MaximumDumpingDepth = 1000.0 + (2500.0 * bd0) / (bd0 + 686.0 * math.exp(-5.63 * bd0))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * bd0) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    out_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n_layers):
        _DepthBottom += states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        bdi = states_BulkDensity[i]
        _MaximumDumpingDepth = 1000.0 + (2500.0 * bdi) / (bdi + 686.0 * math.exp(-5.63 * bdi))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * bdi) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        out_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return out_SoilTemperatureByLayers