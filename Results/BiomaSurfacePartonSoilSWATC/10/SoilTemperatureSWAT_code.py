def Init(states_LayerThickness):
    """
    Initializes the SoilTemperatureByLayers state variable.

    Inputs:
    - states_LayerThickness: list or sequence of floats (m), soil layer thickness for each layer.

    Returns:
    - states_SoilTemperatureByLayers: list of floats (degC), initialized to 15 degC for each layer.
    """
    n_layers = len(states_LayerThickness)
    states_SoilTemperatureByLayers = [15.0 for _ in range(n_layers)]
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
    Soil temperature by layers computed using the SWAT method.
    Reference: Neitsch et al., SWAT Theoretical Documentation, Version 2000.

    Inputs:
    - LagCoefficient: float, dimensionless lag coefficient [0-1]
    - states_SoilProfileDepth: float (m)
    - states_SurfaceSoilTemperature: float (degC)
    - exogenous_AirTemperatureAnnualAverage: float (degC)
    - states_SoilTemperatureByLayers: list of floats (degC), previous day's soil temperature by layer
    - states_BulkDensity: list of floats (t m-3), for each soil layer
    - states_VolumetricWaterContent: list of floats (m3 m-3), for each soil layer
    - states_LayerThickness: list of floats (m), for each soil layer

    Returns:
    - states_SoilTemperatureByLayers: list of floats (degC), updated soil temperature by layer
    """
    import math

    # Copy to avoid mutating the input list in place
    states_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

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

    states_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i]))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        states_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
        )

    return states_SoilTemperatureByLayers