def CalculateModel(a_deficitOnTopLayers: float, soilDiffusionConstant: float = 4.2) -> float:
    """
    Compute diffusion limited evaporation from the soil.

    Inputs:
    - a_deficitOnTopLayers: float
        Deficit on top soil layers (g m-2).
    - soilDiffusionConstant: float (default 4.2)
        Soil diffusion constant (dimensionless as per original code context).

    Returns:
    - diffusionLimitedEvaporation: float
        Evaporation limited by diffusion (g m-2 d-1).
    """
    deficitOnTopLayers = a_deficitOnTopLayers
    diffusionLimitedEvaporation: float
    if (deficitOnTopLayers / 1000.0) <= 0.0:
        diffusionLimitedEvaporation = 8.3 * 1000.0
    else:
        if (deficitOnTopLayers / 1000.0) < 25.0:
            diffusionLimitedEvaporation = (
                2.0 * soilDiffusionConstant * soilDiffusionConstant / (deficitOnTopLayers / 1000.0) * 1000.0
            )
        else:
            diffusionLimitedEvaporation = 0.0
    return diffusionLimitedEvaporation