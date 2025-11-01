def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute soil temperature profile for 1-cm layers.

    Inputs:
      - temp_amp (float): current canopy temperature amplitude (degC)
      - prev_temp_profile (list[float]): previous soil temperature profile (degC) for 1-cm layers
      - prev_canopy_temp (float): previous canopy temperature (degC)
      - min_air_temp (float): current minimum air temperature (degC)

    Returns:
      - temp_profile (list[float]): current soil temperature profile (degC) for 1-cm layers
    """
    import math

    n = len(prev_temp_profile)
    temp_profile = [0.0] * n
    vexp = [0.0] * n

    # Physical constants
    therm_diff = 5.37e-3  # thermal diffusivity (cm^2/s)
    temp_freq = 7.272e-5  # diurnal frequency (1/s)

    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    # Exponential attenuation with depth (Fortran z in [1..n])
    for z in range(1, n + 1):
        vexp[z - 1] = math.exp(-z * therm_amp)

    # Update profile
    for z in range(1, n + 1):
        idx = z - 1
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
    Initialize previous soil temperature profile and previous canopy temperature.

    Inputs:
      - air_temp_day1 (float): initial air temperature (degC)
      - layer_thick (list[int]): thickness of soil layers (cm)

    Returns:
      - prev_temp_profile (list[float]): initialized soil temperature profile (degC) for 1-cm layers
      - prev_canopy_temp (float): initialized canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * max(soil_depth, 0)
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp