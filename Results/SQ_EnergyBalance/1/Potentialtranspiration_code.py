def CalculateModel(r_evapoTranspiration: float, tau: float = 0.9983) -> float:
    """
    Compute potentialTranspiration from evapoTranspiration and plant cover factor (tau).

    Inputs:
    - r_evapoTranspiration: float, evapotranspiration rate.
    - tau: float, plant cover factor (0-1). Default is 0.9983.

    Returns:
    - potentialTranspiration: float, computed potential transpiration.
    """
    potentialTranspiration = r_evapoTranspiration * (1.0 - tau)
    return potentialTranspiration