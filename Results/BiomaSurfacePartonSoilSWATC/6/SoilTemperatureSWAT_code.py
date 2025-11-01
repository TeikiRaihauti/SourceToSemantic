from typing import List, Sequence


def Init(states_LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers state.

    Inputs:
    - states_LayerThickness: sequence of float, m; soil layer thicknesses

    Returns:
    - SoilTemperatureByLayers: list of float, degC; initialized temperature for each layer
    """
    n = len(states_LayerThickness)
    SoilTemperatureByLayers = [15.0 for _ in range(n)]
    return SoilTemperatureByLayers


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: Sequence[float],
    states_BulkDensity: Sequence[float],
    states_VolumetricWaterContent: Sequence[float],
    states_LayerThickness: Sequence[float],
) -> List[float]:
    """
    SoilTemperatureSWAT main process. Update soil temperature by layers using SWAT method.

    Inputs:
    - LagCoefficient: float, dimensionless; lag coefficient [0-1]
    - states_SoilProfileDepth: float, m
    - states_SurfaceSoilTemperature: float, degC
    - exogenous_AirTemperatureAnnualAverage: float, degC
    - states_SoilTemperatureByLayers: sequence of float, degC; previous day's layer temperatures
    - states_BulkDensity: sequence of float, t m-3
    - states_VolumetricWaterContent: sequence of float, m3 m-3
    - states_LayerThickness: sequence of float, m

    Returns:
    - SoilTemperatureByLayers: list of float, degC; updated layer temperatures
    """
    n = len(states_LayerThickness)
    if n == 0:
        return []

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n):
        _TotalWaterContentmm += states_VolumetricWaterContent[i] * states_LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # Prepare output list
    SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * (2.718281828459045 ** (-5.63 * states_BulkDensity[0]))
    )
    _ScalingFactor = _TotalWaterContentmm / (
        (0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm
    )
    _DumpingDepth = _MaximumDumpingDepth * (
        (500.0 / _MaximumDumpingDepth) ** (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + (2.718281828459045 ** (-0.867 - 2.078 * _RatioCenter)))

    SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor
            * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    _DepthBottom = 0.0
    for i in range(1, n):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (
            states_BulkDensity[i] + 686.0 * (2.718281828459045 ** (-5.63 * states_BulkDensity[i]))
        )
        _ScalingFactor = _TotalWaterContentmm / (
            (0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm
        )
        _DumpingDepth = _MaximumDumpingDepth * (
            (500.0 / _MaximumDumpingDepth) ** (((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + (2.718281828459045 ** (-0.867 - 2.078 * _RatioCenter)))

        SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor
                * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return SoilTemperatureByLayers