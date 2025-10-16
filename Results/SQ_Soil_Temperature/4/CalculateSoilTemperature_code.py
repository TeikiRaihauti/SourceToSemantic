import math

def Init(exogenous_meanAnnualAirTemp):
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    return (weatherMinTemp + deepTemperature) / 2.0


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


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999
        states_maxTSoil = -999
        deeplayerstates_deepLayerT = 0.0
    else:
        states_minTSoil = SoilMinimumTemperature(
            exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_, deeplayerstates_deepLayerT
        )
        states_maxTSoil = SoilMaximumTemperature(
            exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_, deeplayerstates_deepLayerT
        )
        deeplayerstates_deepLayerT = UpdateTemperature(states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT)

    return deeplayerstates_deepLayerT, states_maxTSoil, states_minTSoil