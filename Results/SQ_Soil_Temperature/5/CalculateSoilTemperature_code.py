def Init(exogenous_meanAnnualAirTemp: float) -> float:
    """
    Initialize deep soil layer temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: float, annual mean air temperature (°C)

    Returns:
    - deepLayerT: float, initialized deep soil layer temperature (°C)
    """
    deepLayerT = exogenous_meanAnnualAirTemp
    return deepLayerT


def Estimate(
    deeplayerstates_deepLayerT: float,
    exogenous_maxTAir: float,
    exogenous_meanTAir: float,
    exogenous_minTAir: float,
    ratesexternal_heatFlux: float,
    lambda_: float = 2.454,
) -> tuple:
    """
    Calculate minimum and maximum soil temperature and update deep soil temperature.

    Inputs:
    - deeplayerstates_deepLayerT: float, previous deep soil layer temperature (°C)
    - exogenous_maxTAir: float, maximum air temperature (°C)
    - exogenous_meanTAir: float, mean air temperature (°C)
    - exogenous_minTAir: float, minimum air temperature (°C)
    - ratesexternal_heatFlux: float, soil heat flux (g m-2 d-1)
    - lambda_: float, latent heat of water vaporization at 20°C (MJ kg-1), default 2.454

    Returns (tuple):
    - minTSoil: float, minimum soil temperature (°C)
    - maxTSoil: float, maximum soil temperature (°C)
    - deepLayerT: float, updated deep soil layer temperature (°C)
    """
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        minTSoil = 999.0
        maxTSoil = -999.0
        deepLayerT = 0.0
    else:
        minTSoil = SoilMinimumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT,
        )
        maxTSoil = SoilMaximumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT,
        )
        deepLayerT = UpdateTemperature(minTSoil, maxTSoil, deeplayerstates_deepLayerT)

    return minTSoil, maxTSoil, deepLayerT


def SoilMinimumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float,
) -> float:
    """
    Soil minimum temperature.
    OUTPUT UNITS: °C
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature),
    )


def SoilMaximumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float,
) -> float:
    """
    Soil maximum temperature.
    OUTPUT UNITS: °C
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature),
    )


def SoilTempA(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    soilHeatFlux: float,
    lambda_: float,
) -> float:
    """
    Soil temperature A.
    """
    import math

    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp: float, deepTemperature: float) -> float:
    """
    Soil temperature B.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp: float, maxSoilTemp: float, Temperature: float) -> float:
    """
    Update the deep soil layer temperature.

    Inputs:
    - minSoilTemp: float, minimum soil temperature (°C)
    - maxSoilTemp: float, maximum soil temperature (°C)
    - Temperature: float, previous deep soil temperature (°C)

    Returns:
    - Temperature: float, updated deep soil temperature (°C)
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature