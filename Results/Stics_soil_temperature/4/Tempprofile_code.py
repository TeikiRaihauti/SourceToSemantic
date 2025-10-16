def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialize previous soil temperature profile and canopy temperature.

    Parameters
    - air_temp_day1: float, Mean temperature on first day (degC)
    - layer_thick: list of int, layers thickness (cm)

    Returns
    - prev_temp_profile: list of float, initialized soil temperature profile (1 cm layers)
    - prev_canopy_temp: float, initialized canopy temperature
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp


def model_temp_amp(min_temp, max_temp):
    """
    Calculates temperature amplitude.

    Parameters
    - min_temp: float, current minimum temperature (degC)
    - max_temp: float, current maximum temperature (degC)

    Returns
    - temp_amp: float, current temperature amplitude (degC)
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculates canopy mean temperature.

    Parameters
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)

    Returns
    - canopy_temp_avg: float, current canopy mean temperature (degC)
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def get_layers_number(layer_thick_or_depth):
    """
    Counts the number of non-zero soil layers.

    Parameters
    - layer_thick_or_depth: list of int

    Returns
    - layers_number: int
    """
    layers_number = 0
    for z in range(len(layer_thick_or_depth)):
        if layer_thick_or_depth[z] != 0:
            layers_number += 1
    return layers_number


def layer_thickness2depth(layer_thick):
    """
    Converts layer thickness to cumulative depth.

    Parameters
    - layer_thick: list of int

    Returns
    - layer_depth: list of int, cumulative depth at bottom of each layer (same length as layer_thick)
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    cum = 0
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            cum += layer_thick[z]
            layer_depth[z] = cum
        else:
            layer_depth[z] = 0
    return layer_depth


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates soil layers' mean temperatures from the 1-cm profile.

    Parameters
    - temp_profile: list of float, soil temperature profile per 1 cm layer (degC)
    - layer_thick: list of int, layer thicknesses (cm)

    Returns
    - layer_temp: list of float, mean temperature of each soil layer (degC)
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
        end = up_depth[z + 1]
        if end > start and layer_thick[z] != 0:
            layer_temp[z] = sum(temp_profile[start:end]) / float(layer_thick[z])
        else:
            layer_temp[z] = 0.0

    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile.

    Parameters
    - temp_amp: float, current temperature amplitude (degC)
    - prev_temp_profile: list of float, previous soil temperature profile per 1 cm (degC)
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)

    Returns
    - temp_profile: list of float, current soil temperature profile per 1 cm (degC)
    """
    from math import exp, sqrt

    n = len(prev_temp_profile)
    vexp = [0.0] * n
    temp_profile = [0.0] * n

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = sqrt(temp_freq / 2.0 / therm_diff)

    for z in range(n):
        vexp[z] = exp(-(z + 1) * therm_amp)

    for z in range(n):
        temp_profile[z] = (
            prev_temp_profile[z]
            - vexp[z] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[z])
            + (temp_amp * vexp[z]) / 2.0
        )

    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Updates previous state variables.

    Parameters
    - canopy_temp_avg: float, canopy mean temperature (degC)
    - temp_profile: list of float, current soil temperature profile per 1 cm (degC)

    Returns
    - prev_canopy_temp: float, updated previous canopy temperature (degC)
    - prev_temp_profile: list of float, updated previous soil temperature profile per 1 cm (degC)
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Main soil temperature model process for one time step.

    Parameters
    - prev_temp_profile: list of float, previous soil temperature profile per 1 cm (degC) [INOUT]
    - prev_canopy_temp: float, previous canopy temperature (degC) [INOUT]
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)
    - layer_thick: list of int, layer thicknesses (cm)

    Returns
    - prev_temp_profile: list of float, updated previous soil temperature profile per 1 cm (degC)
    - prev_canopy_temp: float, updated previous canopy temperature (degC)
    - temp_profile: list of float, current soil temperature profile per 1 cm (degC)
    - layer_temp: list of float, current layers mean temperature (degC)
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)

    temp_profile = model_temp_profile(
        temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp
    )

    layer_temp = model_layers_temp(temp_profile, layer_thick)

    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)

    prev_canopy_temp, prev_temp_profile = model_update(canopy_temp_avg, temp_profile)

    return prev_temp_profile, prev_canopy_temp, temp_profile, layer_temp