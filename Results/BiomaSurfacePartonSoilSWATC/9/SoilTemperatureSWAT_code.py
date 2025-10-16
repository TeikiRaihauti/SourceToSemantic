from typing import List, Sequence
import math


def Init(states_LayerThickness: Sequence[float]) -> List[float]:
    """
    Initialization function for SoilTemperatureSWAT.
    Creates and initializes the SoilTemperatureByLayers array to 15 degC for each soil layer.

    Inputs:
    - states_LayerThickness: list/sequence of layer thicknesses (m)

    Returns:
    - states_SoilTemperatureByLayers: list of initialized soil temperatures by layer (degC)
    """
    states_SoilTemperatureByLayers = [15.0 for _ in range(len(states_LayerThickness))]
    return states_SoilTemperatureByLayers


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
    Main biophysical process function for SoilTemperatureSWAT.
    Computes soil temperature by layer using the SWAT method.

    Inputs:
    - LagCoefficient: lag coefficient (dimensionless, 0-1)
    - states_SoilProfileDepth: soil profile depth (m)
    - states_SurfaceSoilTemperature: average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: annual average air temperature (degC)
    - states_SoilTemperatureByLayers: previous day's soil temperature by layer (degC)
    - states_BulkDensity: bulk density by layer (t m-3)
    - states_VolumetricWaterContent: volumetric water content by layer (m3 m-3)
    - states_LayerThickness: layer thickness by layer (m)

    Returns:
    - states_SoilTemperatureByLayers: updated soil temperature by layer (degC)
    """

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content in profile (mm)
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

    # Use a copy for pure functional behavior while reading previous values
    prev = list(states_SoilTemperatureByLayers)
    new_temps = list(states_SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[0]) / (
        states_BulkDensity[0] + 686.0 * math.exp(-5.63 * states_BulkDensity[0])
    )
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp(
        (math.log(500.0 / _MaximumDumpingDepth))
        * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
    )
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    new_temps[0] = LagCoefficient * prev[0] + (1.0 - LagCoefficient) * (
        _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
        + states_SurfaceSoilTemperature
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom = _DepthBottom + states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * states_BulkDensity[i]) / (
            states_BulkDensity[i] + 686.0 * math.exp(-5.63 * states_BulkDensity[i])
        )
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * states_BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp(
            (math.log(500.0 / _MaximumDumpingDepth))
            * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0)
        )
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        new_temps[i] = LagCoefficient * prev[i] + (1.0 - LagCoefficient) * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )

    return new_temps