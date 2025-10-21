from typing import Tuple
import math


def Init(exogenous_meanAnnualAirTemp: float) -> float:
    """
    Initialize deep soil layer temperature.

    Inputs:
    - exogenous_meanAnnualAirTemp: float, Annual Mean Air Temperature (°C)

    Returns:
    - deeplayerstates_deepLayerT: float, initialized deep soil layer temperature (°C)
    """
    deeplayerstates_deepLayerT = exogenous_meanAnnualAirTemp
    return deeplayerstates_deepLayerT


def Estimate(
    deeplayerstates_deepLayerT: float,
    exogenous_maxTAir: float,
    exogenous_meanTAir: float,
    exogenous_minTAir: float,
    ratesexternal_heatFlux: float,
    lambda_: float = 2.454,
) -> Tuple[float, float, float]:
    """
    Main biophysical process: calculate daily minimum and maximum soil temperature
    and update deep soil layer temperature.

    Inputs:
    - deeplayerstates_deepLayerT: float, previous deep soil layer temperature (°C)
    - exogenous_maxTAir: float, Maximum Air Temperature (°C)
    - exogenous_meanTAir: float, Mean Air Temperature (°C)
    - exogenous_minTAir: float, Minimum Air Temperature (°C)
    - ratesexternal_heatFlux: float, Soil Heat Flux (g m-2 d-1)
    - lambda_: float, Latent heat of water vaporization at 20°C (MJ kg-1) [default 2.454]

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

    return states_minTSoil, states_maxTSoil, deeplayerstates_deepLayerT


def SoilMinimumTemperature(
    weatherMaxTemp: float,
    weatherMeanTemp: float,
    weatherMinTemp: float,
    soilHeatFlux: float,
    lambda_: float,
    deepTemperature: float,
) -> float:
    """
    Compute daily minimum soil temperature.

    Inputs:
    - weatherMaxTemp: float, maximum air temperature (°C)
    - weatherMeanTemp: float, mean air temperature (°C)
    - weatherMinTemp: float, minimum air temperature (°C)
    - soilHeatFlux: float, soil heat flux (g m-2 d-1)
    - lambda_: float, latent heat of vaporization (MJ kg-1)
    - deepTemperature: float, deep layer temperature (°C)

    Returns:
    - min soil temperature (°C)
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
    Compute daily maximum soil temperature.

    Inputs:
    - weatherMaxTemp: float, maximum air temperature (°C)
    - weatherMeanTemp: float, mean air temperature (°C)
    - weatherMinTemp: float, minimum air temperature (°C)
    - soilHeatFlux: float, soil heat flux (g m-2 d-1)
    - lambda_: float, latent heat of vaporization (MJ kg-1)
    - deepTemperature: float, deep layer temperature (°C)

    Returns:
    - max soil temperature (°C)
    """
    return max(
        SoilTempA(weatherMaxTemp, weatherMeanTemp, soilHeatFlux, lambda_),
        SoilTempB(weatherMinTemp, deepTemperature),
    )


def SoilTempA(
    weatherMaxTemp: float, weatherMeanTemp: float, soilHeatFlux: float, lambda_: float
) -> float:
    """
    Soil temperature component A.

    Inputs:
    - weatherMaxTemp: float, maximum air temperature (°C)
    - weatherMeanTemp: float, mean air temperature (°C)
    - soilHeatFlux: float, soil heat flux (g m-2 d-1)
    - lambda_: float, latent heat of vaporization (MJ kg-1)

    Returns:
    - soil temperature contribution (°C)
    """
    TempAdjustment = -0.5 * weatherMeanTemp + 4.0 if weatherMeanTemp < 8.0 else 0.0
    SoilAvailableEnergy = soilHeatFlux * lambda_ / 1000.0
    return weatherMaxTemp + 11.2 * (1.0 - math.exp(-0.07 * (SoilAvailableEnergy - 5.5))) + TempAdjustment


def SoilTempB(weatherMinTemp: float, deepTemperature: float) -> float:
    """
    Soil temperature component B.

    Inputs:
    - weatherMinTemp: float, minimum air temperature (°C)
    - deepTemperature: float, deep layer temperature (°C)

    Returns:
    - soil temperature contribution (°C)
    """
    return (weatherMinTemp + deepTemperature) / 2.0


def UpdateTemperature(minSoilTemp: float, maxSoilTemp: float, Temperature: float) -> float:
    """
    Update deep soil temperature using exponential smoothing.

    Inputs:
    - minSoilTemp: float, minimum soil temperature (°C)
    - maxSoilTemp: float, maximum soil temperature (°C)
    - Temperature: float, previous deep layer temperature (°C)

    Returns:
    - updated deep layer temperature (°C)
    """
    meanTemp = (minSoilTemp + maxSoilTemp) / 2.0
    Temperature = (9.0 * Temperature + meanTemp) / 10.0
    return Temperature