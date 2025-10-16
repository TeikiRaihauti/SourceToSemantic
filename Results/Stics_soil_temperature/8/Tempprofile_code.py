def model_temp_amp(min_temp, max_temp):
    """
    Calculates temperature amplitude.
    Inputs:
      - min_temp: current minimum temperature (degC)
      - max_temp: current maximum temperature (degC)
    Output:
      - temp_amp: current temperature amplitude (degC)
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculates canopy temperature average.
    Inputs:
      - min_canopy_temp: current minimum canopy temperature (degC)
      - max_canopy_temp: current maximum canopy temperature (degC)
    Output:
      - canopy_temp_avg: current canopy mean temperature (degC)
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def layer_thickness2depth(layer_thick):
    """
    Converts layer thicknesses to cumulative depths.
    Input:
      - layer_thick: list of layer thicknesses (cm)
    Output:
      - layer_depth: list of cumulative bottom depths for layers (cm)
                     For zero-thickness layers, depth is 0.
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    running = 0
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            running += layer_thick[z]
            layer_depth[z] = running
        else:
            layer_depth[z] = 0
    return layer_depth


def get_soil_depth(layer_thick):
    """
    Computes total soil depth (sum of layer thicknesses).
    Input:
      - layer_thick: list of layer thicknesses (cm)
    Output:
      - soil_depth: total depth (cm)
    """
    soil_depth = sum(layer_thick)
    return soil_depth


def get_layers_number(layer_thick_or_depth):
    """
    Counts the number of non-zero layers.
    Input:
      - layer_thick_or_depth: list of layer thicknesses or depths
    Output:
      - layers_number: number of non-zero entries
    """
    layers_number = 0
    for z in range(len(layer_thick_or_depth)):
        if layer_thick_or_depth[z] != 0:
            layers_number += 1
    return layers_number


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates soil layers mean temperature from the 1-cm temperature profile.
    Inputs:
      - temp_profile: soil temperature profile per 1-cm layer (degC), length equals sum(layer_thick)
      - layer_thick: macro-layer thicknesses (cm)
    Output:
      - layer_temp: mean temperature per macro-layer (degC)
    """
    layers_nb = get_layers_number(layer_thick)

    layer_temp = [0.0] * layers_nb
    up_depth = [0] * (layers_nb + 1)

    layer_depth = layer_thickness2depth(layer_thick)

    for z in range(layers_nb):
        depth_value = layer_depth[z]
        up_depth[z + 1] = depth_value

    for z in range(layers_nb):
        depth_start = up_depth[z]
        depth_end = up_depth[z + 1]
        thickness = layer_thick[z] if z < len(layer_thick) else 0
        if thickness > 0 and depth_end > depth_start:
            segment = temp_profile[depth_start:depth_end]
            layer_temp[z] = sum(segment) / float(thickness)
        else:
            layer_temp[z] = 0.0

    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile for 1-cm layers.
    Inputs:
      - temp_amp: current temperature amplitude (degC)
      - prev_temp_profile: previous soil temperature profile per 1-cm layer (degC)
      - prev_canopy_temp: previous crop canopy temperature (degC)
      - min_air_temp: current minimum air temperature (degC)
    Output:
      - temp_profile: current soil temperature profile per 1-cm layer (degC)
    """
    from math import sqrt, exp

    n = len(prev_temp_profile)
    vexp = [0.0] * n
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = sqrt(temp_freq / 2.0 / therm_diff)

    for z in range(1, n + 1):
        vexp[z - 1] = exp(-z * therm_amp)

    temp_profile = [0.0] * n
    for z in range(1, n + 1):
        idx = z - 1
        temp_profile[idx] = (
            prev_temp_profile[idx]
            - vexp[idx] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[idx])
            + (temp_amp * vexp[idx]) / 2.0
        )

    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Updates previous state variables.
    Inputs:
      - canopy_temp_avg: current canopy mean temperature (degC)
      - temp_profile: current soil temperature profile per 1-cm layer (degC)
    Outputs:
      - prev_canopy_temp: updated previous canopy temperature (degC)
      - prev_temp_profile: updated previous soil temperature profile per 1-cm layer (degC)
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initializes soil temperature profile and previous canopy temperature.
    Inputs:
      - air_temp_day1: Mean temperature on first day (degC)
      - layer_thick: list of layer thicknesses (cm)
    Outputs:
      - prev_temp_profile: initialized soil temperature profile per 1-cm layer (degC)
      - prev_canopy_temp: initialized previous canopy temperature (degC)
    """
    soil_depth = get_soil_depth(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Main soil temperature model process.
    Inputs:
      - prev_temp_profile: previous soil temperature profile per 1-cm layer (degC)
      - prev_canopy_temp: previous canopy temperature (degC)
      - min_canopy_temp: current minimum canopy temperature (degC)
      - max_canopy_temp: current maximum canopy temperature (degC)
      - min_air_temp: current minimum air temperature (degC)
      - layer_thick: list of layer thicknesses (cm)
    Outputs:
      - prev_temp_profile: updated previous soil temperature profile per 1-cm layer (degC)
      - prev_canopy_temp: updated previous canopy temperature (degC)
      - temp_profile: current soil temperature profile per 1-cm layer (degC)
      - layer_temp: current layers mean temperature (degC)
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)

    temp_profile = model_temp_profile(
        temp_amp,
        prev_temp_profile,
        prev_canopy_temp,
        min_air_temp
    )

    layer_temp = model_layers_temp(temp_profile, layer_thick)

    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)

    prev_canopy_temp, prev_temp_profile = model_update(canopy_temp_avg, temp_profile)

    return prev_temp_profile, prev_canopy_temp, temp_profile, layer_temp