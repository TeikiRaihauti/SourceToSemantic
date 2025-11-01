def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute soil temperature profile for 1-cm layers.

    Inputs:
    - temp_amp: float
    - prev_temp_profile: list[float]
    - prev_canopy_temp: float
    - min_air_temp: float

    Returns:
    - temp_profile: list[float]
    """
    # Local constants and variables (mapped from Fortran)
    # therm_diff = 5.37e-3, temp_freq = 7.272e-5
    # therm_amp = sqrt(temp_freq/2/therm_diff)
    import math

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    n = len(prev_temp_profile)

    if n == 0:
        return []

    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    # vexp(z) = exp(-z * therm_amp) for z = 1..n (Fortran 1-based indexing)
    vexp = [math.exp(-(z) * therm_amp) for z in range(1, n + 1)]

    # temp_profile(z) = prev_temp_profile(z)
    #                   - vexp(z) * (prev_canopy_temp - min_air_temp)
    #                   + 0.1 * (prev_canopy_temp - prev_temp_profile(z))
    #                   + (temp_amp * vexp(z)) / 2
    temp_profile = [0.0] * n
    for z in range(1, n + 1):
        idx = z - 1  # convert to 0-based index
        temp_profile[idx] = (
            prev_temp_profile[idx]
            - vexp[idx] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[idx])
            + (temp_amp * vexp[idx]) / 2.0
        )

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: list[int],
) -> tuple[list[float], float]:
    """
    Initialize previous soil temperature profile (1-cm resolution) and previous canopy temperature.

    Inputs:
    - air_temp_day1: float
    - layer_thick: list[int]

    Returns:
    - prev_temp_profile: list[float]
    - prev_canopy_temp: float
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * max(soil_depth, 0)
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp