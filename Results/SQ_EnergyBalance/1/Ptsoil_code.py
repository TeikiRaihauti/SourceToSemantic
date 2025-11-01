def CalculateModel(
    r_evapoTranspirationPriestlyTaylor: float,
    Alpha: float = 1.5,
    tau: float = 0.9983,
    tauAlpha: float = 0.3,
) -> float:
    """
    Calculate energy-limited soil evaporation (PtSoil).

    Inputs:
    - r_evapoTranspirationPriestlyTaylor: float
      Evapotranspiration from Priestly-Taylor (g m-2 d-1)
    - Alpha: float
      Priestley-Taylor evapotranspiration proportionality constant (default 1.5)
    - tau: float
      Plant cover factor (default 0.9983)
    - tauAlpha: float
      Fraction of the total net radiation exchanged at the soil surface when AlphaE = 1 (default 0.3)

    Returns:
    - energyLimitedEvaporation: float
      Energy-limited soil evaporation (g m-2 d-1)
    """
    evapoTranspirationPriestlyTaylor: float = r_evapoTranspirationPriestlyTaylor
    energyLimitedEvaporation: float
    AlphaE: float
    if tau < tauAlpha:
        AlphaE = 1.0
    else:
        AlphaE = Alpha - ((Alpha - 1.0) * (1.0 - tau) / (1.0 - tauAlpha))
    energyLimitedEvaporation = evapoTranspirationPriestlyTaylor / Alpha * AlphaE * tau
    return energyLimitedEvaporation