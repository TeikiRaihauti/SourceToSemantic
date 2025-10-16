def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    tempprofile model
    Calculates soil temperature profile for 1 cm layers.

    Inputs:
    - temp_amp: current temperature amplitude (degC)
    - prev_temp_profile: previous soil temperature profile (list of degC for 1 cm layers)
    - prev_canopy_temp: previous crop temperature (degC)
    - min_air_temp: current minimum air temperature (degC)

    Returns:
    - temp_profile: current soil profile temperature (list of degC for 1 cm layers)
    """
    from math import sqrt, exp

    n = len(prev_temp_profile)

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = sqrt(temp_freq / 2.0 / therm_diff)

    # Exponential decay with depth (1 cm indexed)
    vexp = [exp(-(z) * therm_amp) for z in range(1, n + 1)]

    temp_profile = []
    for i in range(n):
        tp = (
            prev_temp_profile[i]
            - vexp[i] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[i])
            + (temp_amp * vexp[i]) / 2.0
        )
        temp_profile.append(tp)

    return temp_profile


def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialization for soil temperature profile.

    Inputs:
    - air_temp_day1: Mean temperature on first day (degC)
    - layer_thick: layers thickness (list of cm)

    Returns:
    - prev_temp_profile: initialized previous soil temperature profile (list of degC for 1 cm layers)
    - prev_canopy_temp: initialized previous canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp