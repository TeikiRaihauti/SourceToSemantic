def Estimate(evapoTranspirationPriestlyTaylor: float, evapoTranspirationPenman: float, isWindVpDefined: int = 1) -> float:
    """
    Main biophysical process function for EvapoTranspiration.

    Inputs:
    - evapoTranspirationPriestlyTaylor: float
      Priestly-Taylor evapotranspiration rate.
    - evapoTranspirationPenman: float
      Penman-Monteith evapotranspiration rate.
    - isWindVpDefined: int (default=1)
      Indicator (1 or 0) whether wind and vapour pressure are defined (1 -> use Penman, 0 -> use Priestly-Taylor).

    Returns:
    - evapoTranspiration: float
      Selected evapotranspiration rate based on data availability.
    """
    return CalculateModel(evapoTranspirationPriestlyTaylor, evapoTranspirationPenman, isWindVpDefined)


def CalculateModel(evapoTranspirationPriestlyTaylor: float, evapoTranspirationPenman: float, isWindVpDefined: int) -> float:
    """
    Supporting function implementing the evapotranspiration selection logic.

    Inputs:
    - evapoTranspirationPriestlyTaylor: float
    - evapoTranspirationPenman: float
    - isWindVpDefined: int

    Returns:
    - evapoTranspiration: float
    """
    if isWindVpDefined == 1:
        evapoTranspiration = evapoTranspirationPenman
    else:
        evapoTranspiration = evapoTranspirationPriestlyTaylor
    return evapoTranspiration