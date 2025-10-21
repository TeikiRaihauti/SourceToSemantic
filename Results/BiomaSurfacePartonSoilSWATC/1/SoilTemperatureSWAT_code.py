from typing import List


def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize soil temperature by layers.

    Inputs:
    - states_LayerThickness: List[float]
        Soil layer thickness for each layer (m)

    Returns:
    - states_SoilTemperatureByLayers: List[float]
        Initialized soil temperature for each layer (degC), set to 15 for all layers
    """
    states_SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float = 0.8,
    states_SoilProfileDepth: float = 0.0,
    states_SurfaceSoilTemperature: float = 0.0,
    exogenous_AirTemperatureAnnualAverage: float = 0.0,
    states_SoilTemperatureByLayers: List[float] = None,
    states_BulkDensity: List[float] = None,
    states_VolumetricWaterContent: List[float] = None,
    states_LayerThickness: List[float] = None,
) -> List[float]:
    """
    Compute soil temperature by layers using the SWAT method.

    Inputs:
    - LagCoefficient: float
        Lag coefficient that controls the influence of the previous day's temperature on the current day's temperature (dimensionless) [default: 0.8]
    - states_SoilProfileDepth: float
        Soil profile depth (m)
    - states_SurfaceSoilTemperature: float
        Average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float
        Annual average air temperature (degC)
    - states_SoilTemperatureByLayers: List[float]
        Soil temperature of each layer from previous day (degC)
    - states_BulkDensity: List[float]
        Bulk density of each layer (t m-3)
    - states_VolumetricWaterContent: List[float]
        Volumetric soil water content of each layer (m3 m-3)
    - states_LayerThickness: List[float]
        Soil layer thickness for each layer (m)

    Returns:
    - states_SoilTemperatureByLayers_new: List[float]
        Updated soil temperature of each layer (degC)
    """
    if states_SoilTemperatureByLayers is None:
        states_SoilTemperatureByLayers = []
    if states_BulkDensity is None:
        states_BulkDensity = []
    if states_VolumetricWaterContent is None:
        states_VolumetricWaterContent = []
    if states_LayerThickness is None:
        states_LayerThickness = []

    n_layers = len(states_LayerThickness)
    states_SoilTemperatureByLayers_new = [0.0] * n_layers

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
        states_BulkDensity[0] + 686.0 * pow(2.718281828, -5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * pow(
        500.0 / _MaximumDumpingDepth,
        pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0),
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + pow(2.718281828, -0.867 - 2.078 * _RatioCenter))

    states_SoilTemperatureByLayers_new[0] = (
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
            states_BulkDensity[i] + 686.0 * pow(2.718281828, -5.63 * states_BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * pow(
            500.0 / _MaximumDumpingDepth,
            pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0),
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + pow(2.718281828, -0.867 - 2.078 * _RatioCenter))

        states_SoilTemperatureByLayers_new[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return states_SoilTemperatureByLayers_new