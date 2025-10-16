import math

def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function.
    Sets the deep soil layer temperature to the annual mean air temperature.
    Inputs:
    - exogenous_meanAnnualAirTemp: Annual Mean Air Temperature (°C)
    Returns:
    - deeplayerstates_deepLayerT: Initial deep layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature A component based on maximum air temperature and soil available energy.
    Inputs:
    - weatherMaxTemp: Maximum Air Temperature (°C)
    - weatherMeanTemp: Mean Air Temperature (°C)
    - soilHeatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    Returns:
    - SoilTempA value (°C)
    """
    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature B component based on minimum air temperature and deep soil temperature.
    Inputs:
    - weatherMinTemp: Minimum Air Temperature (°C)
    - deepTemperature: Deep soil layer temperature (°C)
    Returns:
    - SoilTempB value (°C)
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def SoilMinimumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil minimum temperature.
    Inputs:
    - weatherMaxTemp: Maximum Air Temperature (°C)
    - weatherMeanTemp: Mean Air Temperature (°C)
    - weatherMinTemp: Minimum Air Temperature (°C)
    - soilHeatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    - deepTemperature: Deep soil layer temperature (°C)
    Returns:
    - Minimum soil temperature (°C)
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil maximum temperature.
    Inputs:
    - weatherMaxTemp: Maximum Air Temperature (°C)
    - weatherMeanTemp: Mean Air Temperature (°C)
    - weatherMinTemp: Minimum Air Temperature (°C)
    - soilHeatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    - deepTemperature: Deep soil layer temperature (°C)
    Returns:
    - Maximum soil temperature (°C)
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep layer temperature based on current min/max soil temperature.
    Inputs:
    - minSoilTemp: Minimum soil temperature (°C)
    - maxSoilTemp: Maximum soil temperature (°C)
    - Temperature: Previous deep soil layer temperature (°C)
    Returns:
    - Updated deep soil layer temperature (°C)
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature


def CalculateModel(exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_, deeplayerstates_deepLayerT):
    """
    Main biophysical process function for soil temperature.
    Computes min and max soil temperature and updates the deep layer temperature.
    Inputs:
    - exogenous_maxTAir: Maximum Air Temperature (°C)
    - exogenous_meanTAir: Mean Air Temperature (°C)
    - exogenous_minTAir: Minimum Air Temperature (°C)
    - ratesexternal_heatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    - deeplayerstates_deepLayerT: Current deep soil layer temperature (°C)
    Returns:
    - states_minTSoil: Minimum soil temperature (°C)
    - states_maxTSoil: Maximum soil temperature (°C)
    - deeplayerstates_deepLayerT: Updated deep soil layer temperature (°C)
    """
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999.0
        states_maxTSoil = -999.0
        deeplayerstates_deepLayerT = 0.0
    else:
        states_minTSoil = SoilMinimumTemperature(
            exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir,
            ratesexternal_heatFlux, lambda_, deeplayerstates_deepLayerT
        )
        states_maxTSoil = SoilMaximumTemperature(
            exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir,
            ratesexternal_heatFlux, lambda_, deeplayerstates_deepLayerT
        )
        deeplayerstates_deepLayerT = UpdateTemperature(
            states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT
        )
    return states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Wrapper for the main process function to preserve original structure.
    Inputs:
    - deeplayerstates_deepLayerT: Current deep soil layer temperature (°C)
    - exogenous_maxTAir: Maximum Air Temperature (°C)
    - exogenous_meanTAir: Mean Air Temperature (°C)
    - exogenous_minTAir: Minimum Air Temperature (°C)
    - ratesexternal_heatFlux: Soil Heat Flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)
    Returns:
    - states_minTSoil: Minimum soil temperature (°C)
    - states_maxTSoil: Maximum soil temperature (°C)
    - deeplayerstates_deepLayerT: Updated deep soil layer temperature (°C)
    """
    return CalculateModel(
        exogenous_maxTAir,
        exogenous_meanTAir,
        exogenous_minTAir,
        ratesexternal_heatFlux,
        lambda_,
        deeplayerstates_deepLayerT
    )