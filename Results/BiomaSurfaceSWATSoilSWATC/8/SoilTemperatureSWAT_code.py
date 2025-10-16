def Init(LayerThickness):
    """
    Initialization function for SoilTemperatureSWAT.
    Sets the initial soil temperature for each layer to 15 degC.

    Inputs:
    - LayerThickness: list of layer thicknesses (m)

    Returns:
    - SoilTemperatureByLayers: list of initial soil temperatures by layer (degC)
    """
    SoilTemperatureByLayers = [15.0 for _ in range(len(LayerThickness))]
    return SoilTemperatureByLayers


def Estimate(SoilProfileDepth,
             SurfaceSoilTemperature,
             AirTemperatureAnnualAverage,
             SoilTemperatureByLayers,
             BulkDensity,
             VolumetricWaterContent,
             LayerThickness,
             LagCoefficient):
    """
    Main biophysical process function for SoilTemperatureSWAT.
    Updates soil temperature for each layer based on SWAT methodology.

    Inputs:
    - SoilProfileDepth: soil profile depth (m)
    - SurfaceSoilTemperature: average surface soil temperature (degC)
    - AirTemperatureAnnualAverage: annual average air temperature (degC)
    - SoilTemperatureByLayers: list of current soil temperatures by layer (degC) [state at t-1]
    - BulkDensity: list of bulk density by layer (t m-3)
    - VolumetricWaterContent: list of volumetric soil water content by layer (m3 m-3)
    - LayerThickness: list of layer thickness by layer (m)
    - LagCoefficient: lag coefficient (dimensionless)

    Returns:
    - SoilTemperatureByLayers: list of updated soil temperatures by layer (degC) [state at t]
    """
    import math

    # Conversion to mm
    _SoilProfileDepthmm = SoilProfileDepth * 1000.0

    # Total water content (mm)
    _TotalWaterContentmm = 0.0
    for i in range(len(LayerThickness)):
        _TotalWaterContentmm += VolumetricWaterContent[i] * LayerThickness[i]
    _TotalWaterContentmm = _TotalWaterContentmm * 1000.0

    # Internal variables
    _MaximumDumpingDepth = 0.0
    _DumpingDepth = 0.0
    _ScalingFactor = 0.0
    _DepthBottom = 0.0
    _RatioCenter = 0.0
    _DepthFactor = 0.0

    # Prepare output as a copy to avoid in-place modification of input state
    updated_SoilTemperatureByLayers = list(SoilTemperatureByLayers)

    if len(LayerThickness) == 0:
        return updated_SoilTemperatureByLayers

    # First layer
    _DepthCenterLayer = LayerThickness[0] * 1000.0 / 2.0

    _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[0]) / (BulkDensity[0] + 686.0 * math.exp(-5.63 * BulkDensity[0]))
    _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[0]) * _SoilProfileDepthmm)
    _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0))
    _RatioCenter = _DepthCenterLayer / _DumpingDepth
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

    updated_SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
    )

    # Other layers
    for i in range(1, len(LayerThickness)):
        _DepthBottom = _DepthBottom + LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + LayerThickness[i] * 1000.0 / 2.0
        _MaximumDumpingDepth = 1000.0 + (2500.0 * BulkDensity[i]) / (BulkDensity[i] + 686.0 * math.exp(-5.63 * BulkDensity[i]))
        _ScalingFactor = _TotalWaterContentmm / ((0.356 - 0.144 * BulkDensity[i]) * _SoilProfileDepthmm)
        _DumpingDepth = _MaximumDumpingDepth * math.exp((math.log(500.0 / _MaximumDumpingDepth)) * math.pow((1.0 - _ScalingFactor) / (1.0 + _ScalingFactor), 2.0))
        _RatioCenter = _DepthCenterLayer / _DumpingDepth
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter))

        updated_SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (AirTemperatureAnnualAverage - SurfaceSoilTemperature) + SurfaceSoilTemperature)
        )

    return updated_SoilTemperatureByLayers