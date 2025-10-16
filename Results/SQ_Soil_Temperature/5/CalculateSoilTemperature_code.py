import math

def Init(deeplayerstates_deepLayerT, exogenous_meanAnnualAirTemp):
    """
    Initialize deep soil layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: current deep layer temperature (°C)
    - exogenous_meanAnnualAirTemp: annual mean air temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: initialized deep layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(deeplayerstates_deepLayerT,
             exogenous_maxTAir,
             exogenous_meanTAir,
             exogenous_minTAir,
             ratesexternal_heatFlux,
             exogenous_meanAnnualAirTemp,
             lambda_=2.454):
    """
    Main biophysical process to compute daily minimum and maximum soil temperature
    and update deep layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: previous deep layer temperature (°C)
    - exogenous_maxTAir: maximum air temperature (°C)
    - exogenous_meanTAir: mean air temperature (°C)
    - exogenous_minTAir: minimum air temperature (°C)
    - ratesexternal_heatFlux: soil heat flux (g m-2 d-1)
    - exogenous_meanAnnualAirTemp: annual mean air temperature (°C) [not used in computation]
    - lambda_: latent heat of water vaporization at 20°C (MJ kg-1), default 2.454

    Returns:
    - states_minTSoil: computed minimum soil temperature (°C)
    - states_maxTSoil: computed maximum soil temperature (°C)
    - deeplayerstates_deepLayerT: updated deep layer temperature (°C)
    """
    # Mimic the original code's read of meanAnnualAirTemp (unused in computations)
    _tmp_meanAnnualAirTemp = exogenous_meanAnnualAirTemp  # noqa: F841

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
    Soil minimum temperature (°C).
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(weatherMaxTemp, weatherMeanTemp, weatherMinTemp, soilHeatFlux, lambda_, deepTemperature):
    """
    Soil maximum temperature (°C).
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_):
    """
    Soil temperature A component.
    """
    TempAdjustment = (-0.5 * weatherMeanTemp + 4.0) if (weatherMeanTemp < 8.0) else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp, deepTemperature):
    """
    Soil temperature B component.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp, maxSoilTemp, Temperature):
    """
    Update the deep soil layer temperature (°C) using a weighted average of previous
    deep temperature and current mean soil temperature.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature