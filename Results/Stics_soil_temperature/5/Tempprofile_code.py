def init_temp_profile(air_temp_day1, layer_thick):
    """
    Pure initialization function translated from Fortran SUBROUTINE init_temp_profile.

    Inputs:
      - air_temp_day1: float, Mean temperature on first day
      - layer_thick: sequence of int, layers thickness (cm)

    Returns:
      - prev_temp_profile: list of float, initialized soil temperature profile (1 cm layers)
      - prev_canopy_temp: float, initialized previous crop temperature
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1 for _ in range(soil_depth)]
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Pure main biophysical process function translated from Fortran SUBROUTINE model_temp_profile.

    Inputs:
      - temp_amp: float, current temperature amplitude
      - prev_temp_profile: sequence of float, previous soil temperature profile (1 cm layers)
      - prev_canopy_temp: float, previous crop temperature
      - min_air_temp: float, current minimum air temperature

    Returns:
      - temp_profile: list of float, current soil profile temperature (1 cm layers)
    """
    import math

    n = len(prev_temp_profile)

    # Physical constants (as in original Fortran code)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5

    # Derived parameter
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    # Exponential damping with depth (1-based depth index in Fortran => use (i+1) here)
    vexp = [math.exp(-(i + 1) * therm_amp) for i in range(n)]

    # Compute updated temperature profile
    temp_profile = []
    for i in range(n):
        term1 = prev_temp_profile[i]
        term2 = -vexp[i] * (prev_canopy_temp - min_air_temp)
        term3 = 0.1 * (prev_canopy_temp - prev_temp_profile[i])
        term4 = (temp_amp * vexp[i]) / 2.0
        temp_profile.append(term1 + term2 + term3 + term4)

    return temp_profile