def Init(states_LayerThickness):
    """
    Initializes the SoilTemperatureByLayers array.

    Inputs:
    - states_LayerThickness: list[float]
        Soil layer thickness (m) for each layer.

    Returns:
    - states_SoilTemperatureByLayers: list[float]
        Initialized soil temperature for each layer (degC), set to 15.0.
    """
    n = len(states_LayerThickness)
    states_SoilTemperatureByLayers = [15.0 for _ in range(n)]
    return states_SoilTemperatureByLayers


def Estimate(
    states_SoilProfileDepth,
    states_SurfaceSoilTemperature,
    exogenous_AirTemperatureAnnualAverage,
    states_SoilTemperatureByLayers,
    states_BulkDensity,
    states_VolumetricWaterContent,
    states_LayerThickness,
    LagCoefficient=0.8,
):
    """
    Updates soil temperature by layers using the SWAT approach.

    Inputs:
    - states_SoilProfileDepth: float
        Soil profile depth (m).
    - states_SurfaceSoilTemperature: float
        Average surface soil temperature (degC).
    - exogenous_AirTemperatureAnnualAverage: float
        Annual average air temperature (degC).
    - states_SoilTemperatureByLayers: list[float]
        Previous day's soil temperature for each layer (degC).
    - states_BulkDensity: list[float]
        Bulk density for each layer (t m-3).
    - states_VolumetricWaterContent: list[float]
        Volumetric soil water content for each layer (m3 m-3).
    - states_LayerThickness: list[float]
        Soil layer thickness for each layer (m).
    - LagCoefficient: float (default 0.8)
        Lag coefficient controlling previous day's temperature influence (dimensionless).

    Returns:
    - states_SoilTemperatureByLayers_updated: list[float]
        Updated soil temperature for each layer (degC).
    """
    import math

    n_layers = len(states_LayerThickness)
    states_SoilTemperatureByLayers_updated = list(states_SoilTemperatureByLayers)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    states_SoilTemperatureByLayers_updated[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n_layers):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0

        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (
            states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        states_SoilTemperatureByLayers_updated[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return states_SoilTemperatureByLayers_updated