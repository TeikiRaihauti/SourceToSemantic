import math

def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialize previous soil temperature profile and canopy temperature.
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp


def model_temp_amp(min_temp, max_temp):
    """
    Calculates temperature amplitude.
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculates canopy temperature average.
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def layer_thickness2depth(layer_thick):
    """
    Convert layer thicknesses to cumulative depths.
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            layer_depth[z] = sum(layer_thick[: z + 1])
        else:
            layer_depth[z] = 0
    return layer_depth


def get_layers_number(layer_thick_or_depth):
    """
    Count the number of non-zero layers.
    """
    layers_number = 0
    for z in range(len(layer_thick_or_depth)):
        if layer_thick_or_depth[z] != 0:
            layers_number += 1
    return layers_number


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates soil layers mean temperature.
    """
    layers_nb = get_layers_number(layer_thick)
    layer_temp = [0.0] * layers_nb
    up_depth = [0] * (layers_nb + 1)
    layer_depth = layer_thickness2depth(layer_thick)

    for z in range(layers_nb):
        depth_value = layer_depth[z]
        up_depth[z + 1] = depth_value

    for z in range(layers_nb):
        start = up_depth[z]
        end_ = up_depth[z + 1]
        if end_ > start and layer_thick[z] != 0:
            layer_temp[z] = sum(temp_profile[start:end_]) / float(layer_thick[z])
        else:
            layer_temp[z] = 0.0

    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile at 1 cm resolution.
    """
    n = len(prev_temp_profile)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / (2.0 * therm_diff))

    temp_profile = [0.0] * n

    for i in range(n):
        z = i + 1  # Fortran indices start at 1; z represents depth in cm
        vexp = math.exp(-z * therm_amp)
        temp_profile[i] = (
            prev_temp_profile[i]
            - vexp * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[i])
            + (temp_amp * vexp) / 2.0
        )

    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Update previous state variables.
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Soil temperature main process function. Computes temperature amplitude,
    soil temperature profile, layer mean temperatures, canopy temperature average,
    and updates previous state variables.
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)
    temp_profile = model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp)
    layer_temp = model_layers_temp(temp_profile, layer_thick)
    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)
    prev_canopy_temp, prev_temp_profile = model_update(canopy_temp_avg, temp_profile)

    return prev_temp_profile, prev_canopy_temp, temp_profile, layer_temp