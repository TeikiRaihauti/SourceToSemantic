def CalculateModel(
    psychrometricConstant: float,
    Alpha: float,
    a_netRadiationEquivalentEvaporation: float,
    a_hslope: float,
) -> float:
    """
    PriestlyTaylor main process function.

    Computes evapoTranspirationPriestlyTaylor as:
        max(Alpha * hslope * netRadiationEquivalentEvaporation / (hslope + psychrometricConstant), 0)

    Inputs:
    - psychrometricConstant: float
    - Alpha: float
    - a_netRadiationEquivalentEvaporation: float
    - a_hslope: float

    Returns:
    - r_evapoTranspirationPriestlyTaylor: float
    """
    netRadiationEquivalentEvaporation = a_netRadiationEquivalentEvaporation
    hslope = a_hslope
    evapoTranspirationPriestlyTaylor = max(
        Alpha * hslope * netRadiationEquivalentEvaporation / (hslope + psychrometricConstant),
        0.0,
    )
    return evapoTranspirationPriestlyTaylor