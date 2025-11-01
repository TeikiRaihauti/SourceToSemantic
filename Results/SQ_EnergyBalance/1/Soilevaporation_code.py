def CalculateModel(diffusionLimitedEvaporation: float, energyLimitedEvaporation: float) -> float:
    """
    Compute soil evaporation as the minimum of diffusion-limited and energy-limited evaporation.

    Inputs:
    - diffusionLimitedEvaporation: float
      Evaporation rate limited by diffusion from soil (e.g., g m-2 d-1).
    - energyLimitedEvaporation: float
      Evaporation rate limited by available energy (e.g., g m-2 d-1).

    Returns:
    - soilEvaporation: float
      Soil evaporation (e.g., g m-2 d-1).
    """
    soilEvaporation = min(diffusionLimitedEvaporation, energyLimitedEvaporation)
    return soilEvaporation