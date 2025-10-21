def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute soil temperature profile.

    Parameters:
    - temp_amp: float
        Current temperature amplitude (degC).
    - prev_temp_profile: list[float]
        Previous soil temperature profile (1 cm layers, degC).
    - prev_canopy_temp: float
        Previous canopy temperature (degC).
    - min_air_temp: float
        Current minimum air temperature (degC).

    Returns:
    - temp_profile: list[float]
        Updated soil temperature profile (1 cm layers, degC).
    """
    from math import sqrt, exp

    n = len(prev_temp_profile)

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = sqrt(temp_freq / 2.0 / therm_diff)

    temp_profile = [0.0] * n

    # Fortran loops are 1-based; preserve the exact formulation with z starting at 1
    for z in range(1, n + 1):
        vexp = exp(-z * therm_amp)
        pprev = prev_temp_profile[z - 1]
        temp_profile[z - 1] = (
            pprev
            - vexp * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - pprev)
            + (temp_amp * vexp) / 2.0
        )

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: list[int],
) -> tuple[list[float], float]:
    """
    Initialize soil temperature profile and canopy temperature.

    Parameters:
    - air_temp_day1: float
        Air temperature on day 1 (degC).
    - layer_thick: list[int]
        Soil layer thicknesses (cm).

    Returns:
    - prev_temp_profile: list[float]
        Initialized previous soil temperature profile (1 cm layers, degC).
    - prev_canopy_temp: float
        Initialized previous canopy temperature (degC).
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp