def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute soil temperature profile.

    Inputs:
    - temp_amp: float, temperature amplitude (degC)
    - prev_temp_profile: list[float], previous soil temperature profile (degC) for 1 cm layers
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)

    Returns:
    - temp_profile: list[float], current soil temperature profile (degC) for 1 cm layers
    """
    import math

    n = len(prev_temp_profile)
    temp_profile = [0.0] * n
    vexp = [0.0] * n

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    for z in range(1, n + 1):
        vexp[z - 1] = math.exp(-z * therm_amp)

    for z in range(1, n + 1):
        temp_profile[z - 1] = (
            prev_temp_profile[z - 1]
            - vexp[z - 1] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[z - 1])
            + (temp_amp * vexp[z - 1]) / 2.0
        )

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: list[int],
) -> tuple[list[float], float]:
    """
    Initialize soil temperature profile and previous canopy temperature.

    Inputs:
    - air_temp_day1: float, initial air temperature (degC)
    - layer_thick: list[int], layer thicknesses (cm)

    Returns:
    - prev_temp_profile: list[float], initialized soil temperature profile (degC) for 1 cm layers
    - prev_canopy_temp: float, initialized previous canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp