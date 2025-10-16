def Init(states_LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Sets SoilTemperatureByLayers to 15 degC for each layer, with length equal to states_LayerThickness.

    Inputs:
    - states_LayerThickness: list or sequence of layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of soil temperatures by layer (degC)
    """
    states_SoilTemperatureByLayers = [15 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


def Estimate(LagCoefficient,
             states_SoilProfileDepth,
             states_SurfaceSoilTemperature,
             exogenous_AirTemperatureAnnualAverage,
             states_SoilTemperatureByLayers,
             states_BulkDensity,
             states_VolumetricWaterContent,
             states_LayerThickness):
    """
    Main biophysical process function for SoilTemperatureSWAT.
    Updates soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: float, lag coefficient controlling influence of previous day's temperature (dimensionless)
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of float, previous soil temperature of each layer (degC)
    - states_BulkDensity: list of float, bulk density for each layer (t m-3)
    - states_VolumetricWaterContent: list of float, volumetric soil water content for each layer (m3 m-3)
    - states_LayerThickness: list of float, soil layer thickness for each layer (m)

    Returns:
    - states1_SoilTemperatureByLayers: list of float, updated soil temperature of each layer (degC)
    """
    import math

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    n_layers = len(states_LayerThickness)
    states1_SoilTemperatureByLayers = [0.0] * n_layers

    if n_layers == 0:
        return states1_SoilTemperatureByLayers

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

    states1_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
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

        states1_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return states1_SoilTemperatureByLayers