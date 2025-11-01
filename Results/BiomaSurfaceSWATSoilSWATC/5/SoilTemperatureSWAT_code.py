from typing import List
from math import exp, log


def Init(states_LayerThickness: List[float]) -> List[float]:
    """
    Initialize SoilTemperatureByLayers.

    Inputs:
    - states_LayerThickness: List[float], m

    Returns:
    - states_SoilTemperatureByLayers: List[float], degC
    """
    return [15.0 for _ in range(len(states_LayerThickness))]


def Estimate(
    LagCoefficient: float,
    states_SoilProfileDepth: float,
    states_SurfaceSoilTemperature: float,
    exogenous_AirTemperatureAnnualAverage: float,
    states_SoilTemperatureByLayers: List[float],
    states_BulkDensity: List[float],
    states_VolumetricWaterContent: List[float],
    states_LayerThickness: List[float],
) -> List[float]:
    """
    Update SoilTemperatureByLayers using SWAT method.

    Inputs:
    - LagCoefficient: float, dimensionless
    - states_SoilProfileDepth: float, m
    - states_SurfaceSoilTemperature: float, degC
    - exogenous_AirTemperatureAnnualAverage: float, degC
    - states_SoilTemperatureByLayers: List[float], degC (previous day's values)
    - states_BulkDensity: List[float], t m-3
    - states_VolumetricWaterContent: List[float], m3 m-3
    - states_LayerThickness: List[float], m

    Returns:
    - states_SoilTemperatureByLayers: List[float], degC (updated)
    """
    # Copy input to avoid in-place mutation (functional purity)
    new_SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

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

    if len(states_LayerThickness) == 0:
        return new_SoilTemperatureByLayers

    # First layer
    _DepthCenterLayer = states_LayerThickness[0] * 1000.0 / 2.0
    _MaximumDumpingDepth = _compute_maximum_dumping_depth(states_BulkDensity[0])
    _ScalingFactor = _compute_scaling_factor(_TotalWaterContentmm, _SoilProfileDepthmm, states_BulkDensity[0])
    _DumpingDepth = _compute_dumping_depth(_MaximumDumpingDepth, _ScalingFactor)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _compute_depth_factor(_RatioCenter)

    new_SoilTemperatureByLayers[0] = (
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
        _MaximumDumpingDepth = _compute_maximum_dumping_depth(states_BulkDensity[i])
        _ScalingFactor = _compute_scaling_factor(_TotalWaterContentmm, _SoilProfileDepthmm, states_BulkDensity[i])
        _DumpingDepth = _compute_dumping_depth(_MaximumDumpingDepth, _ScalingFactor)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _compute_depth_factor(_RatioCenter)

        new_SoilTemperatureByLayers[i] = (
            LagCoefficient * states_SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (
                _DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature)
                + states_SurfaceSoilTemperature
            )
        )

    return new_SoilTemperatureByLayers


# Supporting domain-logic functions

def _compute_maximum_dumping_depth(bulk_density: float) -> float:
    # 1000 + (2500 * BD) / (BD + 686 * exp(-5.63 * BD))
    denominator = bulk_density + 686.0 * exp(-5.63 * bulk_density)
    return 1000.0 + (2500.0 * bulk_density) / denominator if denominator != 0.0 else 1000.0


def _compute_scaling_factor(total_water_content_mm: float, soil_profile_depth_mm: float, bulk_density: float) -> float:
    # SF = TWCmm / ((0.356 - 0.144 * BD) * SPDmm)
    denom = (0.356 - 0.144 * bulk_density) * soil_profile_depth_mm
    return total_water_content_mm / denom if denom != 0.0 else 0.0


def _compute_dumping_depth(maximum_dumping_depth: float, scaling_factor: float) -> float:
    # DD = MDD * exp(ln(500/MDD) * ((1 - SF) / (1 + SF))^2)
    if maximum_dumping_depth == 0.0:
        return 0.0
    ratio = 500.0 / maximum_dumping_depth
    if ratio <= 0.0:
        return maximum_dumping_depth
    inner = ((1.0 - scaling_factor) / (1.0 + scaling_factor)) ** 2
    return maximum_dumping_depth * exp(log(ratio) * inner)


def _compute_depth_factor(ratio_center: float) -> float:
    # DF = RC / (RC + exp(-0.867 - 2.078 * RC))
    return ratio_center / (ratio_center + exp(-0.867 - 2.078 * ratio_center)) if ratio_center >= 0.0 else 0.0