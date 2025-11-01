def CalculateModel(s_netRadiation: float, lambdaV: float) -> float:
    """
    Compute netRadiationEquivalentEvaporation from netRadiation and lambdaV.
    Inputs:
    - s_netRadiation: float (MJ m-2 d-1), net radiation at canopy surface
    - lambdaV: float (MJ kg-1), latent heat of vaporization of water
    Returns:
    - netRadiationEquivalentEvaporation: float (g m-2 d-1)
    """
    netRadiation = s_netRadiation
    netRadiationEquivalentEvaporation = netRadiation / lambdaV * 1000.0
    return netRadiationEquivalentEvaporation


def Estimate(s_netRadiation: float, lambdaV: float = 2.454) -> float:
    """
    Main biophysical process function for NetRadiationEquivalentEvaporation.
    Inputs:
    - s_netRadiation: float (MJ m-2 d-1), net radiation at canopy surface
    - lambdaV: float (MJ kg-1), latent heat of vaporization of water (default 2.454)
    Returns:
    - netRadiationEquivalentEvaporation: float (g m-2 d-1)
    """
    return CalculateModel(s_netRadiation, lambdaV)