def Init(states_LayerThickness):
    """
    Initialize the soil temperature profile by layers.

    Inputs:
    - states_LayerThickness: list of soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of initialized soil temperatures (degC) for each layer
    """
    n_layers = 0 if states_LayerThickness is None else len(states_LayerThickness)
    states_SoilTemperatureByLayers = [15.0] * n_layers
    return states_SoilTemperatureByLayers


def Estimate(
    LagCoefficient,
    states_SoilProfileDepth,
    states_SurfaceSoilTemperature,
    exogenous_AirTemperatureAnnualAverage,
    states_SoilTemperatureByLayers,
    states_BulkDensity,
    states_VolumetricWaterContent,
    states_LayerThickness,
):
    """
    Estimate soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: lag coefficient controlling influence of previous day's temperature (dimensionless)
    - states_SoilProfileDepth: soil profile depth (m)
    - states_SurfaceSoilTemperature: average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of previous soil temperatures by layer (degC)
    - states_BulkDensity: list of bulk densities by layer (t m-3)
    - states_VolumetricWaterContent: list of volumetric soil water content by layer (m3 m-3)
    - states_LayerThickness: list of layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: updated list of soil temperatures by layer (degC)
    """
    import math

    # Ensure inputs are valid lists
    if states_LayerThickness is None or len(states_LayerThickness) == 0:
        return []

    n_layers = len(states_LayerThickness)
    # Defensive checks to avoid index errors; assume inputs are consistent as per original behavior
    assert len(states_SoilTemperatureByLayers) >= n_layers
    assert len(states_BulkDensity) >= n_layers
    assert len(states_VolumetricWaterContent) >= n_layers

    # Conversions to mm
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

    updated_SoilTemperatureByLayers = [0.0] * n_layers

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    updated_SoilTemperatureByLayers[0] = (
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
            (math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return updated_SoilTemperatureByLayers