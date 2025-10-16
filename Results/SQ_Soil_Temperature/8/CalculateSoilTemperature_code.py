import math

def Init(deeplayerstates_deepLayerT, exogenous_meanAnnualAirTemp):
    """
    Initialization function setting the deep soil layer temperature to the annual mean air temperature.
    Inputs:
    - deeplayerstates_deepLayerT: current deep soil layer temperature (ignored, replaced)
    - exogenous_meanAnnualAirTemp: annual mean air temperature

    Returns:
    - deeplayerstates_deepLayerT: initialized deep soil layer temperature
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function calculating minimum and maximum soil temperatures and updating deep layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: previous deep soil layer temperature
    - exogenous_maxTAir: daily maximum air temperature
    - exogenous_meanTAir: daily mean air temperature
    - exogenous_minTAir: daily minimum air temperature
    - ratesexternal_heatFlux: soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of vaporization of water (MJ kg-1)

    Returns:
    - states_minTSoil: computed minimum soil temperature
    - states_maxTSoil: computed maximum soil temperature
    - deeplayerstates_deepLayerT: updated deep soil layer temperature
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
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature