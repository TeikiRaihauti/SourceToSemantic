import math

def Init(states_LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.

    Inputs:
    - states_LayerThickness: list of soil layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list initialized with 15.0 degC for each layer
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

    Inputs:
    - LagCoefficient: float, lag coefficient controlling previous day's influence [0-1]
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_SoilTemperatureByLayers: list of floats, previous day's soil temperature per layer (degC)
    - states_BulkDensity: list of floats, bulk density per layer (t m-3)
    - states_VolumetricWaterContent: list of floats, volumetric soil water content per layer (m3 m-3)
    - states_LayerThickness: list of floats, soil layer thickness per layer (m)

    Returns:
    - states_SoilTemperatureByLayers: list of floats, updated soil temperature per layer (degC)
    """

    # Copy to avoid in-place modification of input lists
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

    # First layer center depth (mm)
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    # Maximum dumping depth (mm)
    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0]))
    # Scaling factor based on water content and porosity
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    # Actual dumping depth (mm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    # Ratio of depth center to dumping depth
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    # Depth factor
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    states_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient) * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
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

        states_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient) * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
        )

    return states_SoilTemperatureByLayers