import math

def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function for deep soil layer temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: Annual mean air temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: Initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function.

    Calculates daily minimum and maximum soil temperature and updates the deep soil
    layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: Previous deep soil layer temperature (°C)
    - exogenous_maxTAir: Maximum daily air temperature (°C)
    - exogenous_meanTAir: Mean daily air temperature (°C)
    - exogenous_minTAir: Minimum daily air temperature (°C)
    - ratesexternal_heatFlux: Soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)

    Returns:
    - deeplayerstates_deepLayerT: Updated deep soil layer temperature (°C)
    - states_minTSoil: Daily minimum soil temperature (°C)
    - states_maxTSoil: Daily maximum soil temperature (°C)
    """
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999
        states_maxTSoil = -999
        deeplayerstates_deepLayerT = 0.0
    else:
        states_minTSoil = SoilMinimumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT
        )
        states_maxTSoil = SoilMaximumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT
        )
        deeplayerstates_deepLayerT = UpdateTemperature(
            states_minTSoil,
            states_maxTSoil,
            deeplayerstates_deepLayerT
        )

    return deeplayerstates_deepLayerT, states_minTSoil, states_maxTSoil


def SoilMinimumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Minimum soil temperature (°C) as the minimum of SoilTempA and SoilTempB.
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Maximum soil temperature (°C) as the maximum of SoilTempA and SoilTempB.
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature component A based on maximum air temperature, mean air temperature,
    and soil available energy.
    """
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature component B as the average of minimum air temperature and deep soil temperature.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep soil temperature as a running mean with current mean soil temperature.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature