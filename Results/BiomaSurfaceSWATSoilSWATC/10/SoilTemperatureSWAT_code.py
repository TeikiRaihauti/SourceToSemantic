def Init(states_LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Initializes soil temperature by layers to 15 degC for each layer.

    Inputs:
    - states_LayerThickness: list or sequence of floats, soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of floats, initialized soil temperatures by layer (degC)
    """
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
    Main biophysical process function for SoilTemperatureSWAT.
    Updates soil temperature for each soil layer based on SWAT methodology.

    Inputs:
    - LagCoefficient: float, lag coefficient (0-1) controlling previous day's temperature influence
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of floats, current soil temperature by layer (degC)
    - states_BulkDensity: list of floats, bulk density by layer (t m-3)
    - states_VolumetricWaterContent: list of floats, volumetric water content by layer (m3 m-3)
    - states_LayerThickness: list of floats, layer thickness by layer (m)

    Returns:
    - states_SoilTemperatureByLayers: list of floats, updated soil temperature by layer (degC)
    """
    import math

    # Work on a copy to keep function pure
    SoilTemperatureByLayers_out = list(states_SoilTemperatureByLayers)

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

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0]))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    SoilTemperatureByLayers_out[0] = (
        LagCoefficient * SoilTemperatureByLayers_out[0]
        + (1.0 - LagCoefficient) * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom += states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i]))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        SoilTemperatureByLayers_out[i] = (
            LagCoefficient * SoilTemperatureByLayers_out[i]
            + (1.0 - LagCoefficient) * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers_out