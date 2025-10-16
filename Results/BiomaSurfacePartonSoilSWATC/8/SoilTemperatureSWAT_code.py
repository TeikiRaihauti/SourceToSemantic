def Init(states_LayerThickness):
    """
    Initialization function corresponding to SoilTemperatureSWAT.Init.
    Initializes SoilTemperatureByLayers to 15 degC for each layer.
    
    Inputs:
    - states_LayerThickness: list or array of layer thicknesses (m)
    
    Returns:
    - SoilTemperatureByLayers: list of initialized soil temperatures by layer (degC)
    """
    return [15.0 for _ in range(len(states_LayerThickness))]


def Estimate(LagCoefficient,
             states_SoilProfileDepth,
             states_SurfaceSoilTemperature,
             exogenous_AirTemperatureAnnualAverage,
             states_BulkDensity,
             states_VolumetricWaterContent,
             states_LayerThickness,
             states_SoilTemperatureByLayers):
    """
    Main biophysical process function corresponding to SoilTemperatureSWAT.Estimate.
    Updates soil temperature by layers using the SWAT method.
    
    Inputs:
    - LagCoefficient: float, lag coefficient controlling previous day's influence (0-1)
    - states_SoilProfileDepth: float, soil profile depth (m)
    - states_SurfaceSoilTemperature: float, average surface soil temperature (degC)
    - exogenous_AirTemperatureAnnualAverage: float, annual average air temperature (degC)
    - states_BulkDensity: list/array of bulk density by layer (t m-3)
    - states_VolumetricWaterContent: list/array of volumetric water content by layer (m3 m-3)
    - states_LayerThickness: list/array of soil layer thicknesses (m)
    - states_SoilTemperatureByLayers: list/array of previous soil temperatures by layer (degC)
    
    Returns:
    - SoilTemperatureByLayers: list of updated soil temperatures by layer (degC)
    """
    import math

    # Copy input temperatures to avoid mutating inputs
    SoilTemperatureByLayers = list(states_SoilTemperatureByLayers)

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

    # First layer computations
    bd0 = states_BulkDensity[0]
    _MaximumDumpingDepth = 1000.0 + (2500.0 * bd0) / (bd0 + 686.0 * math.exp(-5.63 * bd0))
    denom0 = (0.356 - 0.144 * bd0) * _SoilProfileDepthmm
    _ScalingFactor = _TotalWaterContentmm / denom0 if denom0 != 0.0 else 0.0
    log_term0 = math.log(500.0 / _MaximumDumpingDepth)
    frac0 = (1.0 - _ScalingFactor) / (1.0 + _ScalingFactor) if (1.0 + _ScalingFactor) != 0.0 else 0.0
    _DumpingDepth = _MaximumDumpingDepth * math.exp(log_term0 * (frac0 ** 2))
    _RatioCenter = _DepthCenterLayer / _DumpingDepth if _DumpingDepth != 0.0 else 0.0
    _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) if (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) != 0.0 else 0.0

    SoilTemperatureByLayers[0] = (
        LagCoefficient * SoilTemperatureByLayers[0]
        + (1.0 - LagCoefficient)
        * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
    )

    # Other layers
    for i in range(1, len(states_LayerThickness)):
        _DepthBottom += states_LayerThickness[i - 1] * 1000.0
        _DepthCenterLayer = _DepthBottom + states_LayerThickness[i] * 1000.0 / 2.0

        bdi = states_BulkDensity[i]
        _MaximumDumpingDepth = 1000.0 + (2500.0 * bdi) / (bdi + 686.0 * math.exp(-5.63 * bdi))
        denomi = (0.356 - 0.144 * bdi) * _SoilProfileDepthmm
        _ScalingFactor = _TotalWaterContentmm / denomi if denomi != 0.0 else 0.0
        log_termi = math.log(500.0 / _MaximumDumpingDepth)
        fraci = (1.0 - _ScalingFactor) / (1.0 + _ScalingFactor) if (1.0 + _ScalingFactor) != 0.0 else 0.0
        _DumpingDepth = _MaximumDumpingDepth * math.exp(log_termi * (fraci ** 2))
        _RatioCenter = _DepthCenterLayer / _DumpingDepth if _DumpingDepth != 0.0 else 0.0
        _DepthFactor = _RatioCenter / (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) if (_RatioCenter + math.exp(-0.867 - 2.078 * _RatioCenter)) != 0.0 else 0.0

        SoilTemperatureByLayers[i] = (
            LagCoefficient * SoilTemperatureByLayers[i]
            + (1.0 - LagCoefficient)
            * (_DepthFactor * (exogenous_AirTemperatureAnnualAverage - states_SurfaceSoilTemperature) + states_SurfaceSoilTemperature)
        )

    return SoilTemperatureByLayers