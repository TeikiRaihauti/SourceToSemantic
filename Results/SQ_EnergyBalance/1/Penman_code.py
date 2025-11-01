from typing import Tuple


def CalculateModel(
    evapoTranspirationPriestlyTaylor: float,
    hslope: float,
    VPDair: float,
    conductance: float,
    psychrometricConstant: float = 0.66,
    Alpha: float = 1.5,
    lambdaV: float = 2.454,
    rhoDensityAir: float = 1.225,
    specificHeatCapacityAir: float = 0.00101,
) -> float:
    """
    Compute evapoTranspirationPenman using the Penman-Monteith formulation.

    Inputs:
    - evapoTranspirationPriestlyTaylor: float
      Priestly-Taylor evapotranspiration (g m-2 d-1)
    - hslope: float
      Slope of saturated vapor pressure curve at air temperature (hPa °C-1)
    - VPDair: float
      Vapor Pressure Deficit of air (hPa)
    - conductance: float
      Boundary layer conductance (m d-1)
    - psychrometricConstant: float, default=0.66
      Psychrometric constant (hPa °C-1)
    - Alpha: float, default=1.5
      Priestley-Taylor proportionality constant (-)
    - lambdaV: float, default=2.454
      Latent heat of vaporization of water (MJ kg-1)
    - rhoDensityAir: float, default=1.225
      Density of air (kg m-3)
    - specificHeatCapacityAir: float, default=0.00101
      Specific heat capacity of dry air (MJ kg-1 °C-1)

    Returns:
    - evapoTranspirationPenman: float
      Penman-Monteith evapotranspiration (g m-2 d-1)

    Notes:
    - May raise ZeroDivisionError if Alpha == 0 or (hslope + psychrometricConstant) == 0.
    """
    evapoTranspirationPenman = evapoTranspirationPriestlyTaylor / Alpha + (
        1000.0
        * (
            rhoDensityAir
            * specificHeatCapacityAir
            * VPDair
            * conductance
            / (lambdaV * (hslope + psychrometricConstant))
        )
    )
    return evapoTranspirationPenman