def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute soil temperature profile for 1-cm layers.

    Inputs:
    - temp_amp: float, temperature amplitude (degC)
    - prev_temp_profile: list[float], previous soil temperature profile (degC) for 1-cm layers
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)

    Returns:
    - temp_profile: list[float], current soil temperature profile (degC) for 1-cm layers
    """
    import math

    n = len(prev_temp_profile)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    temp_profile = [0.0] * n
    for idx in range(n):
        z = idx + 1  # Fortran 1-based index
        vexp = math.exp(-z * therm_amp)
        temp_profile[idx] = (
            prev_temp_profile[idx]
            - vexp * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[idx])
            + (temp_amp * vexp) / 2.0
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
    - layer_thick: list[int], soil layers thickness (cm)

    Returns:
    - prev_temp_profile: list[float], initialized soil temperature profile (degC) for 1-cm layers
    - prev_canopy_temp: float, initialized canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp