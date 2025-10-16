from math import sqrt, exp

def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialize soil temperature profile and canopy temperature.

    Parameters:
    - air_temp_day1: Mean temperature on first day (degC)
    - layer_thick: list/sequence of integers representing layer thicknesses (cm)

    Returns:
    - prev_temp_profile: list of floats, initial soil temperature profile (1 cm layers)
    - prev_canopy_temp: float, initial canopy temperature
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp


def model_soil_temp(prev_temp_profile,
                    prev_canopy_temp,
                    min_canopy_temp,
                    max_canopy_temp,
                    min_air_temp,
                    layer_thick):
    """
    Main soil temperature process for one time step.

    Parameters:
    - prev_temp_profile: list of floats, previous soil temperature profile (degC) at 1 cm resolution
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)
    - layer_thick: list/sequence of integers, soil layer thicknesses (cm)

    Returns:
    - prev_temp_profile: updated previous soil temperature profile for next step (degC)
    - prev_canopy_temp: updated previous canopy temperature for next step (degC)
    - temp_profile: list of floats, current soil profile temperature (degC) at 1 cm resolution
    - layer_temp: list of floats, current layers mean temperature (degC)
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)

    temp_profile = model_temp_profile(temp_amp,
                                      prev_temp_profile,
                                      prev_canopy_temp,
                                      min_air_temp)

    layer_temp = model_layers_temp(temp_profile, layer_thick)

    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)

    prev_canopy_temp, prev_temp_profile = model_update(canopy_temp_avg, temp_profile)

    return prev_temp_profile, prev_canopy_temp, temp_profile, layer_temp


def model_temp_profile(temp_amp,
                       prev_temp_profile,
                       prev_canopy_temp,
                       min_air_temp):
    """
    Calculates soil temperature profile at 1 cm resolution.

    Parameters:
    - temp_amp: float, current temperature amplitude (degC)
    - prev_temp_profile: list of floats, previous soil temperature profile (degC)
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)

    Returns:
    - temp_profile: list of floats, current soil temperature profile (degC)
    """
    n = len(prev_temp_profile)
    temp_profile = [0.0] * n

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = sqrt(temp_freq / (2.0 * therm_diff))

    # Compute vexp for depths 1..n cm
    vexp = [exp(-(z + 1) * therm_amp) for z in range(n)]

    for i in range(n):
        temp_profile[i] = (prev_temp_profile[i]
                           - vexp[i] * (prev_canopy_temp - min_air_temp)
                           + 0.1 * (prev_canopy_temp - prev_temp_profile[i])
                           + (temp_amp * vexp[i]) / 2.0)

    return temp_profile


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates mean temperature per soil layer.

    Parameters:
    - temp_profile: list of floats, soil temperature profile (degC) at 1 cm resolution
    - layer_thick: list/sequence of integers, layer thicknesses (cm)

    Returns:
    - layer_temp: list of floats, mean temperature for each non-zero thickness layer (degC)
    """
    layers_nb = get_layers_number(layer_thick)

    # Up depths (upper boundaries) array, with first value 0 (surface)
    up_depth = [0] * (layers_nb + 1)

    # Bottom depths per layer (cumulative)
    layer_depth = layer_thickness2depth(layer_thick)

    # Fill up_depth using the first layers_nb layers
    for z in range(layers_nb):
        depth_value = layer_depth[z]
        up_depth[z + 1] = depth_value

    # Compute average temperature per layer
    layer_temp = [0.0] * layers_nb
    for z in range(layers_nb):
        start = up_depth[z]
        end = up_depth[z + 1]
        # Sum temperatures from start to end-1 (Python slices are end-exclusive)
        if end > start:
            avg = sum(temp_profile[start:end]) / float(layer_thick[z])
        else:
            avg = 0.0
        layer_temp[z] = avg

    return layer_temp


def layer_thickness2depth(layer_thick):
    """
    Convert layer thicknesses to cumulative depths (bottom depth per layer).

    Parameters:
    - layer_thick: list/sequence of integers, layer thicknesses (cm)

    Returns:
    - layer_depth: list of integers, cumulative depth at bottom of each layer (cm)
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            layer_depth[z] = sum(layer_thick[:(z + 1)])
        else:
            layer_depth[z] = 0
    return layer_depth


def get_soil_depth(layer_thick):
    """
    Calculate total soil depth from layer thicknesses.

    Parameters:
    - layer_thick: list/sequence of integers

    Returns:
    - soil_depth: integer total depth (cm)
    """
    return sum(layer_thick)


def get_layers_number(layer_thick_or_depth):
    """
    Count number of non-zero layers.

    Parameters:
    - layer_thick_or_depth: list/sequence of integers

    Returns:
    - layers_number: integer count of non-zero entries
    """
    layers_number = 0
    for v in layer_thick_or_depth:
        if v != 0:
            layers_number += 1
    return layers_number


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Compute canopy average temperature.

    Parameters:
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)

    Returns:
    - canopy_temp_avg: float, average canopy temperature (degC)
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def model_temp_amp(min_temp, max_temp):
    """
    Compute temperature amplitude.

    Parameters:
    - min_temp: float, current minimum temperature (degC)
    - max_temp: float, current maximum temperature (degC)

    Returns:
    - temp_amp: float, temperature amplitude (degC)
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_update(canopy_temp_avg, temp_profile):
    """
    Update previous state variables.

    Parameters:
    - canopy_temp_avg: float, current canopy mean temperature (degC)
    - temp_profile: list of floats, current soil profile temperature (degC)

    Returns:
    - prev_canopy_temp: float, updated previous canopy temperature (degC)
    - prev_temp_profile: list of floats, updated previous soil temperature profile (degC)
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile