def Init(states_LayerThickness):
    """
    Initialize soil temperature by layers.

    Inputs:
    - states_LayerThickness: list or array of layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of initialized soil temperatures (degC) for each layer
    """
    return [15.0 for _ in range(len(states_LayerThickness))]


def Estimate(
    states_SoilProfileDepth,
    states_SurfaceSoilTemperature,
    exogenous_AirTemperatureAnnualAverage,
    states_SoilTemperatureByLayers,
    states_BulkDensity,
    states_VolumetricWaterContent,
    states_LayerThickness,
    LagCoefficient,
):
    """
    Compute soil temperature by layers using the SWAT method.

    Inputs:
    - states_SoilProfileDepth: soil profile depth (m)
    - states_SurfaceSoilTemperature: average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: annual average air temperature (degC)
    - states_SoilTemperatureByLayers: previous day's soil temperature of each layer (degC), list/array length = number of layers
    - states_BulkDensity: bulk density per layer (t m-3), list/array length = number of layers
    - states_VolumetricWaterContent: volumetric water content per layer (m3 m-3), list/array length = number of layers
    - states_LayerThickness: layer thickness (m), list/array length = number of layers
    - LagCoefficient: lag coefficient (dimensionless, 0-1)

    Returns:
    - states1_SoilTemperatureByLayers: updated soil temperature of each layer (degC), list length = number of layers
    """
    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(states_LayerThickness)):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Prepare output without mutating input temperatures
    states1_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # First layer
    if len(states_LayerThickness) > 0:
        _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
            states_BulkDensity[0] + 686.0 * (2.718281828459045 ** (-5.63 * states_BulkDensity[0]))
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * (
            (500.0 / _MaximumDumpingDepth) ** (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + (2.718281828459045 ** (-0.867 - 2.078 * _RatioCenter)))

        states1_SoilTemperatureByLayers[0] = (
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
            states_BulkDensity[i] + 686.0 * (2.718281828459045 ** (-5.63 * states_BulkDensity[i]))
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * (
            (500.0 / _MaximumDumpingDepth) ** (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + (2.718281828459045 ** (-0.867 - 2.078 * _RatioCenter)))

        states1_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return states1_SoilTemperatureByLayers