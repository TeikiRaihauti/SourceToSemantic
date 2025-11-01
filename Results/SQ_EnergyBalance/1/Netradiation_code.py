def CalculateModel(
    a_minTair: float,
    a_maxTair: float,
    a_solarRadiation: float,
    a_vaporPressure: float,
    a_extraSolarRadiation: float,
    albedoCoefficient: float = 0.23,
    stefanBoltzman: float = 4.903e-09,
    elevation: float = 0.0
) -> tuple:
    """
    Compute net radiation at the canopy surface.

    Domain logic:
    - Shortwave net radiation (Nsr) = (1 - albedoCoefficient) * solarRadiation
    - Clear sky solar radiation = (0.75 + 2e-5 * elevation) * extraSolarRadiation
    - Average fourth-power absolute temperature = ((Tmax+273.16)^4 + (Tmin+273.16)^4) / 2
    - Surface emissivity = 0.34 - 0.14 * sqrt(vaporPressure / 10)
    - Cloud cover factor = 1.35 * (solarRadiation / clearSkySolarRadiation) - 0.35
    - Net outgoing longwave radiation (Nolr) = stefanBoltzman * averageT * surfaceEmissivity * cloudCoverFactor
    - Net radiation = Nsr - Nolr

    Inputs (all floats):
    - a_minTair: minimum air temperature (degC)
    - a_maxTair: maximum air temperature (degC)
    - a_solarRadiation: incoming solar radiation (MJ m-2 d-1)
    - a_vaporPressure: vapor pressure (hPa)
    - a_extraSolarRadiation: extraterrestrial solar radiation (MJ m-2 d-1)
    - albedoCoefficient: canopy albedo (fraction), default 0.23
    - stefanBoltzman: Stefan-Boltzmann constant (MJ K-4 m-2 d-1), default 4.903e-09
    - elevation: site elevation (m), default 0.0

    Returns:
    - netRadiation (float): MJ m-2 d-1
    - netOutGoingLongWaveRadiation (float): MJ m-2 d-1
    """
    minTair = a_minTair
    maxTair = a_maxTair
    solarRadiation = a_solarRadiation
    vaporPressure = a_vaporPressure
    extraSolarRadiation = a_extraSolarRadiation

    Nsr = (1.0 - albedoCoefficient) * solarRadiation
    clearSkySolarRadiation = (0.75 + (2 * (10.0 ** -5) * elevation)) * extraSolarRadiation
    averageT = (((maxTair + 273.16) ** 4) + ((minTair + 273.16) ** 4)) / 2.0
    surfaceEmissivity = 0.34 - (0.14 * ((vaporPressure / 10.0) ** 0.5))
    cloudCoverFactor = 1.35 * (solarRadiation / clearSkySolarRadiation) - 0.35
    Nolr = stefanBoltzman * averageT * surfaceEmissivity * cloudCoverFactor
    netRadiation = Nsr - Nolr
    netOutGoingLongWaveRadiation = Nolr

    return netRadiation, netOutGoingLongWaveRadiation