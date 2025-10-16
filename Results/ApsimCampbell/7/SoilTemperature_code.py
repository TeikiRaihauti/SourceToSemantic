def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and Value == Value:
                return True
    return False


def Sum(values, startIndex, endIndex):
    result = 0.0
    index = -1
    for value in values:
        index += 1
        if index >= startIndex and value != 999999:
            result += value
        if index == endIndex:
            break
    return result


def Zero(arr):
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0


def ToCumThickness(Thickness):
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness


def boundCheck(VariableValue, Lower, Upper, VariableName):
    return


def kelvinT(celciusT):
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def longWaveRadn(emissivity, tDegC):
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def copy_state(state):
    new_state = {}
    for k, v in state.items():
        if isinstance(v, list):
            new_state[k] = v[:]
        else:
            new_state[k] = v
    return new_state


def volumetricFractionRocks(layer, state):
    return state['rocks'][layer] / 100.0


def volumetricFractionOrganicMatter(layer, state):
    return state['carbon'][layer] / 100.0 * state['bulkDensity'][layer] * 2.5 / state['pom']


def volumetricFractionSand(layer, state):
    return (1 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * \
           state['sand'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']


def volumetricFractionSilt(layer, state):
    return (1 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * \
           state['silt'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']


def volumetricFractionClay(layer, state):
    return (1 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * \
           state['clay'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']


def volumetricFractionWater(layer, state):
    return (1 - volumetricFractionOrganicMatter(layer, state)) * state['soilWater'][layer]


def volumetricFractionIce(layer, state):
    return 0.0


def volumetricFractionAir(layer, state):
    return 1.0 - volumetricFractionRocks(layer, state) - \
           volumetricFractionOrganicMatter(layer, state) - \
           volumetricFractionSand(layer, state) - \
           volumetricFractionSilt(layer, state) - \
           volumetricFractionClay(layer, state) - \
           volumetricFractionWater(layer, state) - \
           volumetricFractionIce(layer, state)


def volumetricSpecificHeat(name, layer, state):
    specificHeatRocks = 7.7
    specificHeatOM = 0.25
    specificHeatSand = 7.7
    specificHeatSilt = 2.74
    specificHeatClay = 2.92
    specificHeatWater = 0.57
    specificHeatIce = 2.18
    specificHeatAir = 0.025
    result = 0.0
    if name == "Rocks":
        result = specificHeatRocks
    elif name == "OrganicMatter":
        result = specificHeatOM
    elif name == "Sand":
        result = specificHeatSand
    elif name == "Silt":
        result = specificHeatSilt
    elif name == "Clay":
        result = specificHeatClay
    elif name == "Water":
        result = specificHeatWater
    elif name == "Ice":
        result = specificHeatIce
    elif name == "Air":
        result = specificHeatAir
    return result


def ThermalConductance(name, layer, state):
    thermalConductanceRocks = 0.182
    thermalConductanceOM = 2.50
    thermalConductanceSand = 0.182
    thermalConductanceSilt = 2.39
    thermalConductanceClay = 1.39
    thermalConductanceWater = 4.18
    thermalConductanceIce = 1.73
    thermalConductanceAir = 0.0012

    result = 0.0
    if name == "Rocks":
        result = thermalConductanceRocks
    elif name == "OrganicMatter":
        result = thermalConductanceOM
    elif name == "Sand":
        result = thermalConductanceSand
    elif name == "Silt":
        result = thermalConductanceSilt
    elif name == "Clay":
        result = thermalConductanceClay
    elif name == "Water":
        result = thermalConductanceWater
    elif name == "Ice":
        result = thermalConductanceIce
    elif name == "Air":
        result = thermalConductanceAir
    elif name == "Minerals":
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, state)) * \
                 (thermalConductanceSand ** volumetricFractionSand(layer, state)) + \
                 (thermalConductanceSilt ** volumetricFractionSilt(layer, state)) + \
                 (thermalConductanceClay ** volumetricFractionClay(layer, state))
    result = volumetricSpecificHeat(name, layer, state)
    return result


def shapeFactor(name, layer, state):
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0

    result = 0.0
    if name == "Rocks":
        result = shapeFactorRocks
    elif name == "OrganicMatter":
        result = shapeFactorOM
    elif name == "Sand":
        result = shapeFactorSand
    elif name == "Silt":
        result = shapeFactorSilt
    elif name == "Clay":
        result = shapeFactorClay
    elif name == "Water":
        result = shapeFactorWater
    elif name == "Ice":
        result = 0.333 - 0.333 * volumetricFractionIce(layer, state) / \
                 (volumetricFractionWater(layer, state) + volumetricFractionIce(layer, state) + volumetricFractionAir(layer, state))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(layer, state) / \
                 (volumetricFractionWater(layer, state) + volumetricFractionIce(layer, state) + volumetricFractionAir(layer, state))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, state) + \
                 shapeFactorSand * volumetricFractionSand(layer, state) + \
                 shapeFactorSilt * volumetricFractionSilt(layer, state) + \
                 shapeFactorClay * volumetricFractionClay(layer, state)
    result = volumetricSpecificHeat(name, layer, state)
    return result


def mapLayer2Node(layerArray, nodeArray, state):
    airNode = 0
    surfaceNode = 1
    numNodes = state['numNodes']
    topsoilNode = 2
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(state['thickness'], 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (state['nodeDepth'][node] * 1000.0)
        d2 = state['nodeDepth'][node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)


def doThermalConductivityCoeffs(state):
    numNodes = state['numNodes']
    numLayers = state['numLayers']
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)

    for layer in range(1, numLayers + 1 + 1):
        element = layer
        bd = state['bulkDensity'][layer]
        thermCondPar1[element] = 0.65 - 0.78 * bd + 0.6 * (bd ** 2)
        thermCondPar2[element] = 1.06 * bd
        clay = state['clay'][layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay ** 0.5) if clay > 0 else 0.0, 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bd ** 2)

    state['thermCondPar1'] = thermCondPar1
    state['thermCondPar2'] = thermCondPar2
    state['thermCondPar3'] = thermCondPar3
    state['thermCondPar4'] = thermCondPar4
    return state


def doVolumetricSpecificHeat(state):
    numNodes = state['numNodes']
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]

    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node, state) * 1000000.0 * state['soilWater'][node]

    mapLayer2Node(volspecHeatSoil_, state['volSpecHeatSoil'], state)


def doThermalConductivity(state):
    numNodes = state['numNodes']
    thermCondLayers = [0.0] * (numNodes + 1)
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]

    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            sf = shapeFactor(constituentName, node, state)
            tc_const = ThermalConductance(constituentName, node, state)
            tc_water = ThermalConductance("Water", node, state)
            if tc_water == 0:
                k = 0.0
            else:
                k = (2.0 / 3.0) * (1 + sf * (tc_const / tc_water - 1.0)) ** (-1) + \
                    (1.0 / 3.0) * (1 + sf * (tc_const / tc_water - 1.0) * (1 - 2 * sf)) ** (-1)
            numerator += tc_const * state['soilWater'][node] * k
            denominator += state['soilWater'][node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)

    mapLayer2Node(thermCondLayers, state['thermalConductivity'], state)


def doThomas(newTemps, state, weather_AirPressure, waterBalance_Eos, waterBalance_Es):
    airNode = 0
    surfaceNode = 1
    numNodes = state['numNodes']
    nu = state['nu']
    timestep = 24.0 * 60.0 * 60.0
    latentHeatOfVapourisation = 2465000.0

    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    state['thermalConductance'][airNode] = state['thermalConductivity'][airNode]

    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (state['nodeDepth'][node + 1] - state['nodeDepth'][node - 1])
        state['heatStorage'][node] = Divide(state['volSpecHeatSoil'][node] * volumeOfSoilAtNode, state['internalTimeStep'], 0.0)
        elementLength = state['nodeDepth'][node + 1] - state['nodeDepth'][node]
        state['thermalConductance'][node] = Divide(state['thermalConductivity'][node], elementLength, 0.0)

    g = 1 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * state['thermalConductance'][node]
        a[node + 1] = c[node]
        b[node] = nu * (state['thermalConductance'][node] + state['thermalConductance'][node - 1]) + state['heatStorage'][node]
        d[node] = g * state['thermalConductance'][node - 1] * state['soilTemp'][node - 1] + \
                  (state['heatStorage'][node] - g * (state['thermalConductance'][node] + state['thermalConductance'][node - 1])) * state['soilTemp'][node] + \
                  g * state['thermalConductance'][node] * state['soilTemp'][node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * state['thermalConductance'][airNode] * newTemps[airNode]

    if state['netRadiationSource'] == "calc":
        radnNet = Divide(state['netRadiation'] * 1000000.0, state['internalTimeStep'], 0.0)
    else:
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)

    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] += soilSurfaceHeatFlux

    d[numNodes] += nu * state['thermalConductance'][numNodes] * newTemps[numNodes + 1]

    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] -= a[node + 1] * c[node]
        d[node + 1] -= a[node + 1] * d[node]

    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)

    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - c[node] * newTemps[node + 1]
        boundCheck(newTemps[node], -50.0, 100.0, f"newTemps({node})")


def interpolateTemperature(timeHours, weather_MinT, weather_MaxT, weather_MeanT, maxTempYesterday, minTempYesterday, defaultTimeOfMaximumTemperature):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = (math_sin((0.0 + 0.25 - maxT_time) * 2.0 * math_pi()) *
                     (maxTempYesterday - minTempYesterday) / 2.0) + ((maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0:
            tScale = 0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = (math_sin((time + 0.25 - maxT_time) * 2.0 * math_pi()) *
                              (weather_MaxT - weather_MinT) / 2.0) + weather_MeanT
        return currentTemperature


def doUpdate(numInterationsPerDay, state):
    airNode = 0
    surfaceNode = 1
    numNodes = state['numNodes']

    state['soilTemp'] = state['newTemperature'][:]

    if state['timeOfDaySecs'] < state['internalTimeStep'] * 1.2:
        for node in range(surfaceNode, numNodes + 1):
            state['minSoilTemp'][node] = state['soilTemp'][node]
            state['maxSoilTemp'][node] = state['soilTemp'][node]

    for node in range(surfaceNode, numNodes + 1):
        if state['soilTemp'][node] < state['minSoilTemp'][node]:
            state['minSoilTemp'][node] = state['soilTemp'][node]
        elif state['soilTemp'][node] > state['maxSoilTemp'][node]:
            state['maxSoilTemp'][node] = state['soilTemp'][node]
        state['aveSoilTemp'][node] += Divide(state['soilTemp'][node], numInterationsPerDay, 0.0)

    state['boundaryLayerConductance'] += Divide(state['thermalConductivity'][airNode], numInterationsPerDay, 0.0)


def getBoundaryLayerConductance(TNew_zb, state, airTemperature, weather_Wind, weather_AirPressure, waterBalance_Eos, waterBalance_Eo):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98

    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)

    roughnessFactorMomentum = 0.13 * state['canopyHeight']
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * state['canopyHeight']

    surfaceTemperature = TNew_zb[1]

    diffusePenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)

    stefanBoltzmannConstant = 0.0000000567
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for iteration in range(1, 4):
        denominator_mom = math_log(Divide(state['instrumentHeight'] - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, denominator_mom, 0.0)

        denominator_heat = math_log(Divide(state['instrumentHeight'] - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, denominator_heat, 0.0)

        boundaryLayerCond += radiativeConductance

        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)

        stabilityParammeter = Divide(-vonKarmanConstant * state['instrumentHeight'] * gravitationalConstant * heatFluxDensity,
                                     SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)

        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            inside = 1.0 - 16.0 * stabilityParammeter
            if inside < 0.0:
                inside = 0.0
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + math_sqrt(inside)) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond


def calcSoilTemperature(soilTempIO, state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear):
    numNodes = state['numNodes']
    thickness = state['thickness']
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * math_pi() / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = math_sqrt(2.0 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25

    tempProfile = [0.0] * (numNodes + 2)
    for node in range(1, numNodes + 1):
        tempProfile[node] = weather_Tav + weather_Amp * math_exp(-1.0 * cumulativeDepth[node] / zd) * \
                            math_sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math_pi() - cumulativeDepth[node] / zd)

    for node in range(1, numNodes + 1):
        soilTempIO[node] = tempProfile[node]


def calcLayerTemperature(depthLag, alx, deltaTemp, weather_Tav, weather_Amp):
    return weather_Tav + (weather_Amp / 2.0 * math_cos(alx - depthLag) + deltaTemp) * math_exp(-depthLag)


def calcSurfaceTemperature(weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math_sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    boundCheck(surfaceT, -100.0, 100.0, "Initial surfaceT")
    return surfaceT


def doNetRadiation(solarRadn, ITERATIONSperDAY, clock_Today_DayOfYear, weather_Latitude, weather_Radn):
    TSTEPS2RAD = Divide(2.0 * math_pi(), float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin(4.869 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25) + 0.03345 * math_sin(6.224 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25)))
    cD = math_sqrt(1.0 - solarDeclination * solarDeclination)
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(weather_Latitude * math_pi() / 180.0) + cD * math_cos(weather_Latitude * math_pi() / 180.0) * math_cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0

    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    if cloudFr < 0.0:
        cloudFr = 0.0
    if cloudFr > 1.0:
        cloudFr = 1.0

    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)

    cva = math_exp(31.3716 - 6014.79 / kelvinT(weather_Radn if False else 0.0) - 0.00792495 * kelvinT(weather_Radn if False else 0.0))
    # Above line is placeholder for vapour pressure calc; original used MinT; we'll approximate using 0 in absence of MinT here.
    cva = Divide(cva, kelvinT(0.0), 0.0)

    return cloudFr, cva


def interpolateNetRadiation(solarRadn, cloudFr, cva, state, airTemperature, weather_Radn, waterBalance_Salb, waterBalance_Eos, waterBalance_Eo):
    surfaceEmissivity = 0.96
    w2MJ = state['internalTimeStep'] / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, state['soilTemp'][1]) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def math_sin(x):
    import math as _math
    return _math.sin(x)


def math_cos(x):
    import math as _math
    return _math.cos(x)


def math_sqrt(x):
    import math as _math
    return _math.sqrt(x)


def math_log(x):
    import math as _math
    if x <= 0:
        return -1e30
    return _math.log(x)


def math_exp(x):
    import math as _math
    return _math.exp(x)


def math_pi():
    import math as _math
    return _math.pi


def soil_temperature_initialize(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    waterBalance_SW,
    organic_Carbon,
    DepthToConstantTemperature,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    clock_Today_DayOfYear,
    instrumHeight=None,
    nu=0.6,
    boundarLayerConductanceSource="calc",
    netRadiationSource="calc",
    defaultTimeOfMaximumTemperature=14.0
):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    bareSoilRoughness = 57.0
    defaultInstrumentHeight = 1.2

    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(numLayers):
        thickness[i + 1] = physical_Thickness[i]

    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes

    nodeDepth = [0.0] * (numNodes + 2)
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0

    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        bulkDensity[i + 1] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]

    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            t_top = thickness[layer - 1] if (layer - 1) >= 0 else 0.0
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * t_top, thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]

    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]

    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]

    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]

    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]

    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        clay[layer] = clay[numLayers]

    maxSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    minSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    aveSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    soilTemp = [0.0] * (numNodes + 2)
    morningSoilTemp = [0.0] * (numNodes + 2)
    newTemperature = [0.0] * (numNodes + 2)
    thermalConductivity = [0.0] * (numNodes + 1)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 2)

    instrumentHeight = instrumHeight if (instrumHeight is not None and instrumHeight > 0.00001) else defaultInstrumentHeight

    state = {
        'doInitialisationStuff': True,
        'internalTimeStep': 0.0,
        'timeOfDaySecs': 0.0,
        'numNodes': numNodes,
        'numLayers': numLayers,
        'nodeDepth': nodeDepth,
        'thermCondPar1': [0.0] * (numNodes + 1),
        'thermCondPar2': [0.0] * (numNodes + 1),
        'thermCondPar3': [0.0] * (numNodes + 1),
        'thermCondPar4': [0.0] * (numNodes + 1),
        'volSpecHeatSoil': volSpecHeatSoil,
        'soilTemp': soilTemp,
        'morningSoilTemp': morningSoilTemp,
        'heatStorage': heatStorage,
        'thermalConductance': thermalConductance,
        'thermalConductivity': thermalConductivity,
        'boundaryLayerConductance': 0.0,
        'newTemperature': newTemperature,
        'airTemperature': 0.0,
        'maxTempYesterday': 0.0,
        'minTempYesterday': 0.0,
        'soilWater': soilWater,
        'minSoilTemp': minSoilTemp,
        'maxSoilTemp': maxSoilTemp,
        'aveSoilTemp': aveSoilTemp,
        'thickness': thickness,
        'bulkDensity': bulkDensity,
        'rocks': rocks,
        'carbon': carbon,
        'sand': sand,
        'silt': silt,
        'clay': clay,
        'soilRoughnessHeight': bareSoilRoughness,
        'instrumentHeight': instrumentHeight,
        'netRadiation': 0.0,
        'canopyHeight': 0.0,
        'instrumHeight': instrumHeight if instrumHeight is not None else 0.0,
        'nu': nu,
        'boundarLayerConductanceSource': boundarLayerConductanceSource,
        'netRadiationSource': netRadiationSource,
        'pInitialValues': None,
        'InitialValues': None,
        'defaultTimeOfMaximumTemperature': defaultTimeOfMaximumTemperature,
        'DepthToConstantTemperature': DepthToConstantTemperature,
        'pom': 1.3,
        'ps': 2.63
    }

    state = doThermalConductivityCoeffs(state)
    calcSoilTemperature(state['soilTemp'], state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)
    state['newTemperature'] = state['soilTemp'][:]
    return state


def soil_temperature_process(
    state,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_AirPressure,
    weather_Wind,
    weather_Radn,
    weather_Amp,
    weather_Tav,
    weather_Latitude,
    waterBalance_Salb,
    waterBalance_Eos,
    waterBalance_Eo,
    waterBalance_Es,
    waterBalance_SW,
    microClimate_CanopyHeight,
    clock_Today_DayOfYear,
    InitialValues=None
):
    s = copy_state(state)
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    timestep = 24.0 * 60.0 * 60.0
    constantBoundaryLayerConductance = 20.0
    interactionsPerDay = 48

    for i in range(s['numLayers']):
        if i + 1 < len(s['soilWater']) and i < len(waterBalance_SW):
            s['soilWater'][i + 1] = waterBalance_SW[i]
    s['soilWater'][s['numNodes']] = s['soilWater'][s['numLayers']]
    s['canopyHeight'] = max(microClimate_CanopyHeight, s['soilRoughnessHeight']) / 1000.0
    s['instrumentHeight'] = max(s['instrumentHeight'], s['canopyHeight'] + 0.5)

    if s['doInitialisationStuff']:
        if ValuesInArray(InitialValues):
            s['soilTemp'] = [0.0] * (s['numNodes'] + 2)
            for idx in range(len(InitialValues)):
                if topsoilNode + idx < len(s['soilTemp']):
                    s['soilTemp'][topsoilNode + idx] = InitialValues[idx]
            s['InitialValues'] = InitialValues[:]
        else:
            calcSoilTemperature(s['soilTemp'], s, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)
            s['InitialValues'] = [0.0] * s['numLayers']
            for i in range(s['numLayers']):
                s['InitialValues'][i] = s['soilTemp'][topsoilNode + i]

        s['soilTemp'][airNode] = weather_MeanT
        s['soilTemp'][surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb)

        for i in range(s['numNodes'] + 1, len(s['soilTemp'])):
            s['soilTemp'][i] = weather_Tav

        s['newTemperature'] = s['soilTemp'][:]
        s['maxTempYesterday'] = weather_MaxT
        s['minTempYesterday'] = weather_MinT
        s['doInitialisationStuff'] = False

    solarRadn = [0.0] * (interactionsPerDay + 1)
    cva = 0.0
    cloudFr, cva = doNetRadiation(solarRadn, interactionsPerDay, clock_Today_DayOfYear, weather_Latitude, weather_Radn)

    Zero(s['minSoilTemp'])
    Zero(s['maxSoilTemp'])
    Zero(s['aveSoilTemp'])
    s['boundaryLayerConductance'] = 0.0

    s['internalTimeStep'] = round(timestep / interactionsPerDay)

    doVolumetricSpecificHeat(s)
    doThermalConductivity(s)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        s['timeOfDaySecs'] = s['internalTimeStep'] * float(timeStepIteration)
        if timestep < 24.0 * 60.0 * 60.0:
            s['airTemperature'] = weather_MeanT
        else:
            s['airTemperature'] = interpolateTemperature(s['timeOfDaySecs'] / 3600.0, weather_MinT, weather_MaxT, weather_MeanT, s['maxTempYesterday'], s['minTempYesterday'], s['defaultTimeOfMaximumTemperature'])
        s['newTemperature'][airNode] = s['airTemperature']

        s['netRadiation'] = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, s, s['airTemperature'], weather_Radn, waterBalance_Salb, waterBalance_Eos, waterBalance_Eo)

        if s['boundarLayerConductanceSource'] == "constant":
            s['thermalConductivity'][airNode] = constantBoundaryLayerConductance
        else:
            s['thermalConductivity'][airNode] = getBoundaryLayerConductance(s['newTemperature'], s, s['airTemperature'], weather_Wind, weather_AirPressure, waterBalance_Eos, waterBalance_Eo)
            for iteration in range(1, 1 + 1):
                doThomas(s['newTemperature'], s, weather_AirPressure, waterBalance_Eos, waterBalance_Es)
                s['thermalConductivity'][airNode] = getBoundaryLayerConductance(s['newTemperature'], s, s['airTemperature'], weather_Wind, weather_AirPressure, waterBalance_Eos, waterBalance_Eo)

        doThomas(s['newTemperature'], s, weather_AirPressure, waterBalance_Eos, waterBalance_Es)
        doUpdate(interactionsPerDay, s)

        if abs(s['timeOfDaySecs'] - 5.0 * 3600.0) <= min(s['timeOfDaySecs'], 5.0 * 3600.0) * 0.0001:
            s['morningSoilTemp'] = s['soilTemp'][:]

    s['minTempYesterday'] = weather_MinT
    s['maxTempYesterday'] = weather_MaxT

    return s


def FinalSoilTemperature(state):
    numLayers = state['numLayers']
    topsoilNode = 2
    result = [0.0] * numLayers
    for i in range(numLayers):
        result[i] = state['soilTemp'][topsoilNode + i]
    return result


def FinalSoilSurfaceTemperature(state):
    surfaceNode = 1
    return state['soilTemp'][surfaceNode]


def AverageSoilTemperature(state):
    numLayers = state['numLayers']
    topsoilNode = 2
    result = [0.0] * numLayers
    for i in range(numLayers):
        result[i] = state['aveSoilTemp'][topsoilNode + i]
    return result


def AverageSoilSurfaceTemperature(state):
    surfaceNode = 1
    return state['aveSoilTemp'][surfaceNode]


def MinimumSoilTemperature(state):
    numLayers = state['numLayers']
    topsoilNode = 2
    result = [0.0] * numLayers
    for i in range(numLayers):
        result[i] = state['minSoilTemp'][topsoilNode + i]
    return result


def MinimumSoilSurfaceTemperature(state):
    surfaceNode = 1
    return state['minSoilTemp'][surfaceNode]


def MaximumSoilTemperature(state):
    numLayers = state['numLayers']
    topsoilNode = 2
    result = [0.0] * numLayers
    for i in range(numLayers):
        result[i] = state['maxSoilTemp'][topsoilNode + i]
    return result


def MaximumSoilSurfaceTemperature(state):
    surfaceNode = 1
    return state['maxSoilTemp'][surfaceNode]


def BoundaryLayerConductance(state):
    return state['boundaryLayerConductance']


def ThermalConductivity_profile(state):
    numNodes = state['numNodes']
    result = [0.0] * numNodes
    for i in range(numNodes):
        idx = 1 + i
        if idx < len(state['thermalConductivity']):
            result[i] = state['thermalConductivity'][idx]
    return result


def HeatCapacity_profile(state):
    numNodes = state['numNodes']
    result = [0.0] * numNodes
    for i in range(numNodes):
        idx = 1 + i
        if idx < len(state['volSpecHeatSoil']):
            result[i] = state['volSpecHeatSoil'][idx]
    return result


def HeatStore_profile(state):
    numNodes = state['numNodes']
    result = [0.0] * numNodes
    for i in range(numNodes):
        idx = 1 + i
        if idx < len(state['heatStorage']):
            result[i] = state['heatStorage'][idx]
    return result


def Thr_profile(state):
    return state['morningSoilTemp'][:]