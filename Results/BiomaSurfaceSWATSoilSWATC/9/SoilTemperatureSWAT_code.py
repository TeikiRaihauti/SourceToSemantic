def Init(LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Inputs:
    - LayerThickness: list of soil layer thicknesses (m)
    Returns:
    - SoilTemperatureByLayers: list initialized to 15.0 degC for each layer
    """
    return [15.0 for _ in range(len(LayerThickness))]


def Estimate(LagCoefficient,
             SoilProfileDepth,
             SurfaceSoilTemperature,
             AirTemperatureAnnualAverage,
             SoilTemperatureByLayers,
             BulkDensity,
             VolumetricWaterContent,
             LayerThickness):
    """
    Main biophysical process function for SoilTemperatureSWAT.

    Inputs:
    - LagCoefficient: float, dimensionless (0-1)
    - SoilProfileDepth: float, m
    - SurfaceSoilTemperature: float, degC
    - AirTemperatureAnnualAverage: float, degC
    - SoilTemperatureByLayers: list of floats, degC (previous day's values)
    - BulkDensity: list of floats, t m-3 (same length as layers)
    - VolumetricWaterContent: list of floats, m3 m-3 (same length as layers)
    - LayerThickness: list of floats, m (same length as layers)

    Returns:
    - SoilTemperatureByLayers: list of floats, degC (updated values)
    """
    import math

    n_layers = len(LayerThickness)
    if n_layers == 0:
        return []

    # Conversion to mm
    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(n_layers):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm *= 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # Prepare output using previous day's values as base
    updated_SoilTemperatureByLayers = list(SoilTemperatureByLayers)

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0]) / (BulkDensity[0] + 686.0 * math.exp(-5.63 * BulkDensity[0]))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    updated_SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient) * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    # Other layers
    for i in range(1, n_layers):
        _DepthBottom = _DepthBottom + LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0

        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i]) / (BulkDensity[i] + 686.0 * math.exp(-5.63 * BulkDensity[i]))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * ((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor)) ** 2)
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient) * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return updated_SoilTemperatureByLayers