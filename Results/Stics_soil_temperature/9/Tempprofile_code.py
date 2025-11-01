def model_temp_profile(temp_amp: float,
                       prev_temp_profile: list[float],
                       prev_canopy_temp: float,
                       min_air_temp: float) -> list[float]:
    """
    Compute soil temperature profile for 1-cm layers.

    Inputs:
    - temp_amp: float
        Temperature amplitude (degC).
    - prev_temp_profile: list[float]
        Previous soil temperature profile for 1-cm layers (degC), length n.
    - prev_canopy_temp: float
        Previous canopy temperature (degC).
    - min_air_temp: float
        Current minimum air temperature (degC).

    Returns:
    - temp_profile: list[float]
        Updated soil temperature profile for 1-cm layers (degC), length n.
    """
    import math

    n = len(prev_temp_profile)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    vexp = [math.exp(-(z) * therm_amp) for z in range(1, n + 1)]
    temp_profile = [0.0] * n

    for idx in range(n):
        temp_profile[idx] = (
            prev_temp_profile[idx]
            - vexp[idx] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[idx])
            + (temp_amp * vexp[idx]) / 2.0
        )

    return temp_profile


def init_temp_profile(air_temp_day1: float,
                      layer_thick: list[int]) -> tuple[list[float], float]:
    """
    Initialize soil temperature profile and previous canopy temperature.

    Inputs:
    - air_temp_day1: float
        Air temperature of day 1 used for initialization (degC).
    - layer_thick: list[int]
        Soil layer thicknesses (cm). Sum determines profile depth in cm.

    Returns:
    - prev_temp_profile: list[float]
        Initialized soil temperature profile for 1-cm layers (degC), length sum(layer_thick).
    - prev_canopy_temp: float
        Initialized previous canopy temperature (degC).
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp