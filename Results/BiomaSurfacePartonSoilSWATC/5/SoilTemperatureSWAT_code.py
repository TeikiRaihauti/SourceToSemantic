def Init(states_LayerThickness: list[float]) -> list[float]:
    """
    Initialize soil temperature by layers.

    Inputs:
    - states_LayerThickness: list[float], soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list[float], initialized soil temperatures by layer (degC)
    """
    states_SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: list[float],
    states_BulkDensity: list[float],
    states_VolumetricWaterContent: list[float],
    states_LayerThickness: list[float],
) -> list[float]:
    """
    Update soil temperature by layers using SWAT method.

    Inputs:
    - LagCoefficient: float, lag coefficient (dimensionless)
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list[float], previous day's soil temperature by layer (degC)
    - states_BulkDensity: list[float], bulk density by layer (t m-3)
    - states_VolumetricWaterContent: list[float], volumetric water content by layer (m3 m-3)
    - states_LayerThickness: list[float], layer thickness by layer (m)

    Returns:
    - states_SoilTemperatureByLayers: list[float], updated soil temperatures by layer (degC)
    """
    import math

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # Result copy
    result_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / (
        (0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm
    )
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    result_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0

        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (
            states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth))
            * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        result_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return result_SoilTemperatureByLayers