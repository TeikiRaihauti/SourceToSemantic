def Init(exogenous_meanAnnualAirTemp: float) -> float:
    """
    Initialize deep soil layer temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: float, °C

    Returns:
    - deeplayerstates_deepLayerT: float, °C
    """
    return exogenous_meanAnnualAirTemp


def CalculateModel(
    deeplayerstates_deepLayerT: float,
    exogenous_maxTAir: float,
    exogenous_meanTAir: float,
    exogenous_minTAir: float,
    ratesexternal_heatFlux: float,
    lambda_: float = 2.454,
) -> tuple:
    """
    Calculate minimum and maximum soil temperature and update deep layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: float, °C (previous deep soil layer temperature)
    - exogenous_maxTAir: float, °C
    - exogenous_meanTAir: float, °C
    - exogenous_minTAir: float, °C
    - ratesexternal_heatFlux: float, g m-2 d-1 (soil heat flux as water mass equivalent)
    - lambda_: float, MJ kg-1 (latent heat of vaporization at 20°C), default 2.454

    Returns:
    - states_minTSoil: float, °C
    - states_maxTSoil: float, °C
    - deeplayerstates_deepLayerT: float, °C (updated)
    """
    if exogenous_maxTAir == -999 and exogenous_minTAir == 999:
        states_minTSoil = 999.0
        states_maxTSoil = -999.0
        deeplayerstates_deepLayerT = 0.0
    else:
        states_minTSoil = SoilMinimumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT,
        )
        states_maxTSoil = SoilMaximumTemperature(
            exogenous_maxTAir,
            exogenous_meanTAir,
            exogenous_minTAir,
            ratesexternal_heatFlux,
            lambda_,
            deeplayerstates_deepLayerT,
        )
        deeplayerstates_deepLayerT = UpdateTemperature(
            states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT
        )

    return (states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT)


def SoilMinimumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float,
) -> float:
    """
    Compute soil minimum temperature (°C).
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
    Compute soil maximum temperature (°C).
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature),
    )


def SoilTempA(
    weatherMaxTemp: float, weatherMeanTemp: float, soilHeatFlux: float, lambda_: float
) -> float:
    """
    Soil temperature contribution A (°C), driven by available energy and mean air temperature.
    """
    import math

    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp: float, deepTemperature: float) -> float:
    """
    Soil temperature contribution B (°C), average of minimum air temp and deep layer temp.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp: float, maxSoilTemp: float, Temperature: float) -> float:
    """
    Update deep soil layer temperature (°C) using a 10-day exponential smoothing equivalent.

    Inputs:
    - minSoilTemp: float, °C
    - maxSoilTemp: float, °C
    - Temperature: float, °C (previous deep layer temperature)

    Returns:
    - Temperature: float, °C (updated deep layer temperature)
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature