def CalculateModel(
    a_netRadiationEquivalentEvaporation: float,
    r_soilHeatFlux: float,
    r_potentialTranspiration: float,
) -> float:
    """
    Calculate crop heat flux.

    Inputs:
    - a_netRadiationEquivalentEvaporation: float
      Net radiation in equivalent evaporation (g m-2 d-1)
    - r_soilHeatFlux: float
      Soil heat flux (g m-2 d-1)
    - r_potentialTranspiration: float
      Potential transpiration (g m-2 d-1)

    Returns:
    - cropHeatFlux: float
      Crop heat flux (g m-2 d-1)
    """
    netRadiationEquivalentEvaporation = a_netRadiationEquivalentEvaporation
    soilHeatFlux = r_soilHeatFlux
    potentialTranspiration = r_potentialTranspiration
    cropHeatFlux = netRadiationEquivalentEvaporation - soilHeatFlux - potentialTranspiration
    return cropHeatFlux