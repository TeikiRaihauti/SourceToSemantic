def Init(exogenous_meanAnnualAirTemp):
    """
    Initialization function corresponding to CalculateSoilTemperature.Init.
    Sets the deep soil layer temperature to the annual mean air temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: Annual mean air temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: Initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT, exogenous_maxTAir, exogenous_meanTAir, exogenous_minTAir, ratesexternal_heatFlux, lambda_):
    """
    Main biophysical process function corresponding to CalculateSoilTemperature.CalculateModel/Estimate.
    Computes daily minimum and maximum soil temperature and updates the deep soil layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: Previous deep soil layer temperature (°C)
    - exogenous_maxTAir: Maximum air temperature (°C)
    - exogenous_meanTAir: Mean air temperature (°C)
    - exogenous_minTAir: Minimum air temperature (°C)
    - ratesexternal_heatFlux: Soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)

    Returns:
    - states_minTSoil: Minimum soil temperature (°C)
    - states_maxTSoil: Maximum soil temperature (°C)
    - deeplayerstates_deepLayerT: Updated deep soil layer temperature (°C)
    """
    # Missing data sentinel handling
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999
        states_maxTSoil = -999
        deeplayerstates_deepLayerT = 0.0
        return states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT

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


def SoilMinimumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Computes daily minimum soil temperature.
    Returns the minimum of SoilTempA and SoilTempB.
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Computes daily maximum soil temperature.
    Returns the maximum of SoilTempA and SoilTempB.
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature component A.

    Inputs:
    - weatherMaxTemp: Maximum air temperature (°C)
    - weatherMeanTemp: Mean air temperature (°C)
    - soilHeatFlux: Soil heat flux (g m-2 d-1)
    - lambda_: Latent heat of water vaporization at 20°C (MJ kg-1)

    Returns:
    - Soil temperature estimate (°C)
    """
    import math
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature component B.

    Inputs:
    - weatherMinTemp: Minimum air temperature (°C)
    - deepTemperature: Deep soil layer temperature (°C)

    Returns:
    - Soil temperature estimate (°C)
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Updates the deep soil layer temperature from minimum and maximum soil temperatures.

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