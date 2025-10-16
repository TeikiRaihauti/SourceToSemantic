def Init(states_LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Creates the SoilTemperatureByLayers array with default value 15 degC for each soil layer.

    Inputs:
    - states_LayerThickness: list of float, soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of float, initialized soil temperature for each layer (degC)
    """
    states_SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
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
    Computes the soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: float (dimensionless), lag coefficient controlling previous day's influence [0..1]
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of float, previous soil temperature by layers (degC)
    - states_BulkDensity: list of float, bulk density by layers (t m-3)
    - states_VolumetricWaterContent: list of float, volumetric soil water content by layers (m3 m-3)
    - states_LayerThickness: list of float, soil layer thickness by layers (m)

    Returns:
    - states_SoilTemperatureByLayers: list of float, updated soil temperature by layers (degC)
    """
    import math

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm): sum over layers of (theta * thickness) in meters, then convert to mm
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

    n = len(states_LayerThickness)
    # Work on a copy to keep function pure (no in-place modification of inputs)
    new_soil_temp = list(states_SoilTemperatureByLayers)

    if n > 0:
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
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        new_soil_temp[0] = (
            LagCoefficient * states_SoilTemperatureByLayers[0]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    # Other layers
    for i in range(1, n):
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
            * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        new_soil_temp[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    states_SoilTemperatureByLayers = new_soil_temp
    return states_SoilTemperatureByLayers