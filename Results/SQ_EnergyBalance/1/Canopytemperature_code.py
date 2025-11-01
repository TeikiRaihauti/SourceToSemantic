def CalculateModel(
    minTair: float,
    maxTair: float,
    cropHeatFlux: float,
    conductance: float,
    lambdaV: float = 2.454,
    rhoDensityAir: float = 1.225,
    specificHeatCapacityAir: float = 0.00101,
) -> (float, float):
    """
    Compute minimal and maximal canopy temperature from crop heat flux and boundary layer conductance.

    Inputs:
    - minTair: float, minimum air temperature (degC)
    - maxTair: float, maximum air temperature (degC)
    - cropHeatFlux: float, crop heat flux
    - conductance: float, boundary layer conductance
    - lambdaV: float, latent heat of vaporization of water (MJ/kg), default 2.454
    - rhoDensityAir: float, density of air (kg/m^3), default 1.225
    - specificHeatCapacityAir: float, specific heat capacity of dry air (MJ/kg/degC), default 0.00101

    Returns:
    - minCanopyTemperature: float (degC)
    - maxCanopyTemperature: float (degC)
    """
    minCanopyTemperature: float = minTair + (
        cropHeatFlux
        / (rhoDensityAir * specificHeatCapacityAir * conductance / lambdaV * 1000.0)
    )
    maxCanopyTemperature: float = maxTair + (
        cropHeatFlux
        / (rhoDensityAir * specificHeatCapacityAir * conductance / lambdaV * 1000.0)
    )
    return minCanopyTemperature, maxCanopyTemperature