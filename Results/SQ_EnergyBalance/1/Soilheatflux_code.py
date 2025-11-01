def Estimate(tau: float, s_netRadiationEquivalentEvaporation: float, s_soilEvaporation: float) -> float:
    """
    Main biophysical process function for SoilHeatFlux.

    Inputs:
    - tau: float
      Description: plant cover factor (dimensionless)
    - s_netRadiationEquivalentEvaporation: float
      Description: net radiation in equivalent evaporation (g m-2 d-1)
    - s_soilEvaporation: float
      Description: soil evaporation (g m-2 d-1)

    Returns:
    - soilHeatFlux: float
      Description: soil heat flux (g m-2 d-1)
    """
    return CalculateModel(s_netRadiationEquivalentEvaporation, s_soilEvaporation, tau)


def CalculateModel(s_netRadiationEquivalentEvaporation: float, s_soilEvaporation: float, tau: float) -> float:
    """
    Supporting domain logic function for SoilHeatFlux.

    Inputs:
    - s_netRadiationEquivalentEvaporation: float
      Description: net radiation in equivalent evaporation (g m-2 d-1)
    - s_soilEvaporation: float
      Description: soil evaporation (g m-2 d-1)
    - tau: float
      Description: plant cover factor (dimensionless)

    Returns:
    - soilHeatFlux: float
      Description: soil heat flux (g m-2 d-1)
    """
    netRadiationEquivalentEvaporation = s_netRadiationEquivalentEvaporation
    soilEvaporation = s_soilEvaporation
    soilHeatFlux = tau * netRadiationEquivalentEvaporation - soilEvaporation
    return soilHeatFlux