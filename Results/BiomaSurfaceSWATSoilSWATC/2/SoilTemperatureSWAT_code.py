def Init(states_LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Creates and returns the initial SoilTemperatureByLayers array, with each layer set to 15 degC.

    Inputs:
    - states_LayerThickness: list or sequence of layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of soil temperatures by layer (degC)
    """
    if states_LayerThickness is None:
        return []
    return [15.0 for _ in range(len(states_LayerThickness))]


def Estimate(LagCoefficient,
             states_SoilProfileDepth,
             states_SurfaceSoilTemperature,
             exogenous_AirTemperatureAnnualAverage,
             states_SoilTemperatureByLayers,
             states_BulkDensity,
             states_VolumetricWaterContent,
             states_LayerThickness):
    """
    Main biophysical process function: SoilTemperatureSWAT.
    Updates the soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: float, lag coefficient (dimensionless)
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of float, previous day's soil temperature by layer (degC)
    - states_BulkDensity: list of float, bulk density by layer (t m-3)
    - states_VolumetricWaterContent: list of float, volumetric water content by layer (m3 m-3)
    - states_LayerThickness: list of float, layer thickness by layer (m)

    Returns:
    - states_SoilTemperatureByLayers: list of float, updated soil temperature by layer (degC)
    """
    import math

    n_layers = 0 if states_LayerThickness is None else len(states_LayerThickness)
    if n_layers == 0:
        return []

    # Use copies to avoid mutating inputs
    SoilTemperatureByLayers_prev = list(states_SoilTemperatureByLayers) if states_SoilTemperatureByLayers is not None else [15.0] * n_layers
    BulkDensity = list(states_BulkDensity)
    VolumetricWaterContent = list(states_VolumetricWaterContent)
    LayerThickness = list(states_LayerThickness)

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content in mm
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

    # Output temperatures
    SoilTemperatureByLayers_new = [0.0] * n_layers

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0]) / (BulkDensity[0] + 686.0 * math.exp(-5.63 * BulkDensity[0]))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers_new[0] = (
        LagCoefficient * SoilTemperatureByLayers_prev[0] +
        (1.0 - LagCoefficient) * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) +
            states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, n_layers):
        _DepthBottom = _DepthBottom + LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i]) / (BulkDensity[i] + 686.0 * math.exp(-5.63 * BulkDensity[i]))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers_new[i] = (
            LagCoefficient * SoilTemperatureByLayers_prev[i] +
            (1.0 - LagCoefficient) * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) +
                states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers_new