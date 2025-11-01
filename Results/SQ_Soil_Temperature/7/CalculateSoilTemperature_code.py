from typing import Tuple
import math

def Init(exogenous_meanAnnualAirTemp: float) -> float:
    """
    Initialization function.

    Inputs:
    - exogenous_meanAnnualAirTemp: float, Annual Mean Air Temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: float, initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(
    exogenous_maxTAir: float,
    exogenous_meanTAir: float,
    exogenous_minTAir: float,
    ratesexternal_heatFlux: float,
    deeplayerstates_deepLayerT: float,
    lambda_: float = 2.454
) -> Tuple[float, float, float]:
    """
    Main biophysical process function.

    Inputs:
    - exogenous_maxTAir: float, Maximum Air Temperature (°C)
    - exogenous_meanTAir: float, Mean Air Temperature (°C)
    - exogenous_minTAir: float, Minimum Air Temperature (°C)
    - ratesexternal_heatFlux: float, Soil Heat Flux (g m-2 d-1)
    - deeplayerstates_deepLayerT: float, previous deep soil layer temperature (°C)
    - lambda_: float, Latent heat of water vaporization at 20°C (MJ kg-1), default 2.454

    Returns:
    - states_minTSoil: float, Minimum Soil Temperature (°C)
    - states_maxTSoil: float, Maximum Soil Temperature (°C)
    - deeplayerstates_deepLayerT: float, updated deep soil layer temperature (°C)
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


def SoilMinimumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float
) -> float:
    """
    Compute soil minimum temperature (°C).
    """
    return min(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilMaximumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float
) -> float:
    """
    Compute soil maximum temperature (°C).
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature)
    )


def SoilTempA(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    soilHeatFlux: float,
    lambda_: float
) -> float:
    """
    Soil temperature component A.
    """
    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp: float, deepTemperature: float) -> float:
    """
    Soil temperature component B.
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp: float, maxSoilTemp: float, Temperature: float) -> float:
    """
    Update the deep soil layer temperature (°C) based on min and max soil temperatures.
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature