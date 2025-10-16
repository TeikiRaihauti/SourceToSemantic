def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function corresponding to CalculateSoilTemperature.Init.
    Sets the deep soil layer temperature to the mean annual air temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: Annual Mean Air Temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: Initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature A component.
    """
    import math
    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature B component.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def SoilMinimumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil minimum temperature.
    OUTPUT UNITS: °C
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil maximum temperature.
    OUTPUT UNITS: °C
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep layer temperature from soil temperatures.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function corresponding to CalculateSoilTemperature.Estimate.

    Inputs:
    - deeplayerstates_deepLayerT: Current deep soil layer temperature (°C)
    - exogenous_maxTAir: Maximum daily air temperature (°C)
    - exogenous_meanTAir: Mean daily air temperature (°C)
    - exogenous_minTAir: Minimum daily air temperature (°C)
    - ratesexternal_heatFlux: Soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)

    Returns:
    - states_minTSoil: Computed minimum soil temperature (°C)
    - states_maxTSoil: Computed maximum soil temperature (°C)
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
        deeplayerstates_deepLayerT = UpdateTemperature(
            states_minTSoil,
            states_maxTSoil,
            deeplayerstates_deepLayerT
        )

    return states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT