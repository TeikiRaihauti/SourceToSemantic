def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function corresponding to CalculateSoilTemperature.Init.
    Sets the deep soil layer temperature to the annual mean air temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: Annual Mean Air Temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function corresponding to CalculateSoilTemperature.Estimate/CalculateModel.
    Calculates daily minimum and maximum soil temperatures and updates deep soil layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: current deep soil layer temperature (°C)
    - exogenous_maxTAir: maximum air temperature (°C)
    - exogenous_meanTAir: mean air temperature (°C)
    - exogenous_minTAir: minimum air temperature (°C)
    - ratesexternal_heatFlux: soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)

    Returns:
    - states_minTSoil: minimum soil temperature (°C)
    - states_maxTSoil: maximum soil temperature (°C)
    - deeplayerstates_deepLayerT: updated deep soil layer temperature (°C)
    """
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999
        states_maxTSoil = -999
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


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature A component driven by available energy and air temperature.
    """
    import math
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature B component as mean of minimum air and deep soil temperature.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep soil layer temperature as a running mean of soil mean temperature.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature