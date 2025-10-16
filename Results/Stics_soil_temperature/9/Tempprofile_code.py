import math

def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile for 1 cm layers.

    Parameters
    - temp_amp: current temperature amplitude (degC)
    - prev_temp_profile: list/sequence of previous soil temperature profile (degC) for 1 cm layers
    - prev_canopy_temp: previous crop temperature (degC)
    - min_air_temp: current minimum air temperature (degC)

    Returns
    - temp_profile: list of current soil profile temperature (degC) for 1 cm layers
    """
    n = len(prev_temp_profile)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5

    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    temp_profile = [0.0] * n
    for z in range(1, n + 1):
        vexp = math.exp(-z * therm_amp)
        prev_temp = prev_temp_profile[z - 1]
        temp_profile[z - 1] = (
            prev_temp
            - vexp * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp)
            + (temp_amp * vexp) / 2.0
        )

    return temp_profile


def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialize soil temperature profile and canopy temperature.

    Parameters
    - air_temp_day1: Mean air temperature on first day (degC)
    - layer_thick: sequence of layer thicknesses (cm)

    Returns
    - prev_temp_profile: list with soil temperature profile (degC) at 1 cm resolution
    - prev_canopy_temp: previous canopy temperature (degC), initialized to air_temp_day1
    """
    soil_depth = int(sum(layer_thick))
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp