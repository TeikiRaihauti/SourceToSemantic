import math

def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function.
    Sets the deep soil layer temperature to the annual mean air temperature.
    
    Inputs:
    - exogenous_meanAnnualAirTemp: Annual Mean Air Temperature (°C)
    
    Returns:
    - deeplayerstates_deepLayerT: Initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def CalculateModel(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function.
    Calculates minimum and maximum soil temperature and updates deep soil layer temperature.
    
    Inputs:
    - deeplayerstates_deepLayerT: Previous deep soil layer temperature (°C)
    - exogenous_maxTAir: Maximum Air Temperature (°C)
    - exogenous_meanTAir: Mean Air Temperature (°C)
    - exogenous_minTAir: Minimum Air Temperature (°C)
    - ratesexternal_heatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    
    Returns:
    - states_minTSoil: Minimum Soil Temperature (°C)
    - states_maxTSoil: Maximum Soil Temperature (°C)
    - deeplayerstates_deepLayerT: Updated deep soil layer temperature (°C)
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
        deeplayerstates_deepLayerT = UpdateTemperature(states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT)

    return states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT


def SoilMinimumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil minimum temperature (°C)
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil maximum temperature (°C)
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature A
    """
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature B
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep soil layer temperature (°C) from min and max soil temperatures.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature