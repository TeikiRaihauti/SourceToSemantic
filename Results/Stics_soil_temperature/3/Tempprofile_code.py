import math

def model_temp_amp(min_temp, max_temp):
    """
    Calculates temperature amplitude.
    Inputs:
      - min_temp: current minimum temperature
      - max_temp: current maximum temperature
    Returns:
      - temp_amp: current temperature amplitude
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculates canopy average temperature.
    Inputs:
      - min_canopy_temp: current minimum canopy temperature
      - max_canopy_temp: current maximum canopy temperature
    Returns:
      - canopy_temp_avg: canopy mean temperature
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def layer_thickness2depth(layer_thick):
    """
    Convert layer thickness to cumulative depth per position.
    Inputs:
      - layer_thick: list of layer thicknesses (cm), may include trailing zeros
    Returns:
      - layer_depth: list where each non-zero thickness index contains cumulative depth (cm), zeros elsewhere
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    running_sum = 0
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            running_sum += layer_thick[z]
            layer_depth[z] = running_sum
        else:
            layer_depth[z] = 0
    return layer_depth


def get_layers_number(layer_thick_or_depth):
    """
    Count non-zero entries in a list (used to determine number of layers).
    Inputs:
      - layer_thick_or_depth: list of integers
    Returns:
      - layers_number: count of non-zero entries
    """
    layers_number = 0
    for val in layer_thick_or_depth:
        if val != 0:
            layers_number += 1
    return layers_number


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates soil layers mean temperature.
    Inputs:
      - temp_profile: soil temperature profile (1 cm layers)
      - layer_thick: layers thicknesses (cm), trailing zeros allowed
    Returns:
      - layer_temp: soil layers mean temperatures
    """
    layers_nb = get_layers_number(layer_thick)

    # Prepare cumulative depths using the original convention:
    # Consider the first 'layers_nb' entries as effective layers (typical STICS usage: trailing zeros).
    layer_depth_full = layer_thickness2depth(layer_thick)
    up_depth = [0] * (layers_nb + 1)
    for z in range(1, layers_nb + 1):
        depth_value = layer_depth_full[z - 1]
        up_depth[z] = depth_value

    # Calculate mean temperature for each layer using the first 'layers_nb' thickness entries
    layer_temp = [0.0] * layers_nb
    for z in range(layers_nb):
        start = up_depth[z]
        end = up_depth[z + 1]
        thickness = layer_thick[z]
        # Avoid division by zero; in correct inputs thickness > 0
        if thickness > 0 and end > start:
            segment_sum = sum(temp_profile[start:end])
            layer_temp[z] = segment_sum / float(thickness)
        else:
            layer_temp[z] = 0.0
    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile.
    Inputs:
      - temp_amp: current temperature amplitude
      - prev_temp_profile: previous soil temperature profile (1 cm layers)
      - prev_canopy_temp: previous crop temperature
      - min_air_temp: current minimum air temperature
    Returns:
      - temp_profile: current soil profile temperature (1 cm layers)
    """
    n = len(prev_temp_profile)

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    vexp = [0.0] * n
    for i in range(n):
        z = i + 1  # Fortran 1-based depth index
        vexp[i] = math.exp(-z * therm_amp)

    temp_profile = [0.0] * n
    for i in range(n):
        temp_profile[i] = (
            prev_temp_profile[i]
            - vexp[i] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[i])
            + (temp_amp * vexp[i]) / 2.0
        )
    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Update previous state variables.
    Inputs:
      - canopy_temp_avg: current canopy mean temperature
      - temp_profile: current soil profile temperature (1 cm layers)
    Returns:
      - prev_canopy_temp: updated previous canopy temperature
      - prev_temp_profile: updated previous soil temperature profile
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialization function.
    Inputs:
      - air_temp_day1: Mean temperature on first day
      - layer_thick: layers thickness (cm)
    Returns:
      - prev_temp_profile: initialized soil temperature profile (1 cm layers)
      - prev_canopy_temp: initialized previous canopy temperature
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * soil_depth
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Main biophysical process function: soil temperature model.
    Inputs:
      - prev_temp_profile: previous soil temperature profile (1 cm layers)
      - prev_canopy_temp: previous crop temperature
      - min_canopy_temp: current minimum canopy temperature
      - max_canopy_temp: current maximum canopy temperature
      - min_air_temp: current minimum air temperature
      - layer_thick: layers thickness (cm)
    Returns:
      - prev_temp_profile: updated soil temperature profile (1 cm layers)
      - prev_canopy_temp: updated crop temperature
      - temp_profile: current soil profile temperature (1 cm layers)
      - layer_temp: current layers mean temperature
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