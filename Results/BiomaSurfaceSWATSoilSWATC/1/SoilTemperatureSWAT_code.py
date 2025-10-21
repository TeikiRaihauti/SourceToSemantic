from typing import List, Sequence
import math


def Init(
    states_LayerThickness: Sequence[float],
    initial_temperature: float = 15.0,
) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: Sequence[float], m, soil layer thicknesses
    - initial_temperature: float, degC, initial temperature for all layers (default 15.0)

    Returns:
    - states_SoilTemperatureByLayers: List[float], degC, initialized soil temperature per layer
    """
    return [float(initial_temperature) for _ in states_LayerThickness]


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
    SWAT soil temperature by layers update.

    Inputs:
    - LagCoefficient: float, dimensionless
    - states_SoilProfileDepth: float, m
    - states_SurfaceSoilTemperature: float, degC
    - exogenous_AirTemperatureAnnualAverage: float, degC
    - states_SoilTemperatureByLayers: Sequence[float], degC, previous day's soil temperature by layer
    - states_BulkDensity: Sequence[float], t m-3
    - states_VolumetricWaterContent: Sequence[float], m3 m-3
    - states_LayerThickness: Sequence[float], m

    Returns:
    - states_SoilTemperatureByLayers: List[float], degC, updated soil temperature by layers
    """

    # Conversion to mm
    _SoilProfileDepthmm = states_SoilProfileDepth * 1000.0

    # Total water content in mm
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

    # Prepare output as a copy to keep function pure
    updated_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    bd0 = states_BulkDensity[0]

    _MaximumDumpingDepth = 1000.0 + (2500.0 * bd0) / (bd0 + 686.0 * math.exp(-5.63 * bd0))
    denom = (0.356 - 0.144 * bd0) * _SoilProfileDepthmm
    _ScalingFactor = _TotalWaterContentmm / denom if denom != 0.0 else 0.0
    log_term = math.log(500.0 / _MaximumDumpingDepth) if _MaximumDumpingDepth != 0.0 else 0.0
    frac = (1.0 - _ScalingFactor) / (1.0 + _ScalingFactor) if (1.0 + _ScalingFactor) != 0.0 else 0.0
    _DumpingDepth = _MaximumDumpingDepth * math.exp(log_term * (frac ** 2))
    _RatioCenter = _DepthCenterLayer / _DumpingDepth if _DumpingDepth != 0.0 else 0.0
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) if (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) != 0.0 else 0.0

    updated_SoilTemperatureByLayers[0] = (
        LagCoefficient * states_SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (
            _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
            + states_SurfaceSoilTemperature
        )
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom += states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0
        bdi = states_BulkDensity[i]

        _MaximumDumpingDepth = 1000.0 + (2500.0 * bdi) / (bdi + 686.0 * math.exp(-5.63 * bdi))
        denom = (0.356 - 0.144 * bdi) * _SoilProfileDepthmm
        _ScalingFactor = _TotalWaterContentmm / denom if denom != 0.0 else 0.0
        log_term = math.log(500.0 / _MaximumDumpingDepth) if _MaximumDumpingDepth != 0.0 else 0.0
        frac = (1.0 - _ScalingFactor) / (1.0 + _ScalingFactor) if (1.0 + _ScalingFactor) != 0.0 else 0.0
        _DumpingDepth = _MaximumDumpingDepth * math.exp(log_term * (frac ** 2))
        _RatioCenter = _DepthCenterLayer / _DumpingDepth if _DumpingDepth != 0.0 else 0.0
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) if (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) != 0.0 else 0.0

        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return updated_SoilTemperatureByLayers