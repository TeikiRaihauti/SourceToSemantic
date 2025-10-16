import math

def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal

def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            try:
                if Value != 999999 and not math.isnan(Value):
                    return True
            except Exception:
                if Value is not None and Value != 999999:
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
    CumThickness = [0.0 for _ in range(len(Thickness))]
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness

def kelvinT(celciusT):
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin

def longWaveRadn(emissivity, tDegC):
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)

def boundCheck(VariableValue, Lower, Upper, VariableName):
    precisionMargin = 0.00001
    lowerBound = Lower - precisionMargin
    upperBound = Upper + precisionMargin
    # No exception raising per original commented code
    return

def volumetricFractionRocks(layer, state):
    return state['rocks'][layer] / 100.0

def volumetricFractionOrganicMatter(layer, state):
    return state['carbon'][layer] / 100.0 * state['bulkDensity'][layer] * 2.5 / state['pom']

def volumetricFractionSand(layer, state):
    return (1.0 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * state['sand'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']

def volumetricFractionSilt(layer, state):
    return (1.0 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * state['silt'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']

def volumetricFractionClay(layer, state):
    return (1.0 - volumetricFractionOrganicMatter(layer, state) - volumetricFractionRocks(layer, state)) * state['clay'][layer] / 100.0 * state['bulkDensity'][layer] / state['ps']

def volumetricFractionWater(layer, state):
    return (1.0 - volumetricFractionOrganicMatter(layer, state)) * state['soilWater'][layer]

def volumetricFractionIce(layer, state):
    return 0.0

def volumetricFractionAir(layer, state):
    return 1.0 - volumetricFractionRocks(layer, state) - volumetricFractionOrganicMatter(layer, state) - volumetricFractionSand(layer, state) - volumetricFractionSilt(layer, state) - volumetricFractionClay(layer, state) - volumetricFractionWater(layer, state) - volumetricFractionIce(layer, state)

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
        # Per original odd formulation comment
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, state)) * (thermalConductanceSand ** volumetricFractionSand(layer, state)) + (thermalConductanceSilt ** volumetricFractionSilt(layer, state)) + (thermalConductanceClay ** volumetricFractionClay(layer, state))
    # Original code erroneously set result = volumetricSpecificHeat, replicate
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
        denom = (volumetricFractionWater(layer, state) + volumetricFractionIce(layer, state) + volumetricFractionAir(layer, state))
        result = 0.333 - 0.333 * Divide(volumetricFractionIce(layer, state), denom, 0.0)
        return result
    elif name == "Air":
        denom = (volumetricFractionWater(layer, state) + volumetricFractionIce(layer, state) + volumetricFractionAir(layer, state))
        result = 0.333 - 0.333 * Divide(volumetricFractionAir(layer, state), denom, 0.0)
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, state) + shapeFactorSand * volumetricFractionSand(layer, state) + shapeFactorSilt * volumetricFractionSilt(layer, state) + shapeFactorClay * volumetricFractionClay(layer, state)
    # Original code erroneously set result = volumetricSpecificHeat, replicate
    result = volumetricSpecificHeat(name, layer, state)
    return result

def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)

def mapLayer2Node(layerArray, state, nodeArray):
    surfaceNode = state['surfaceNode']
    numNodes = state['numNodes']
    thickness = state['thickness']
    nodeDepth = state['nodeDepth']
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)
    return nodeArray

def doThermalConductivityCoeffs(state):
    numLayers = state['numLayers']
    numNodes = state['numNodes']
    clay = state['clay']
    bulkDensity = state['bulkDensity']

    thermCondPar1 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar2 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar3 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar4 = [0.0 for _ in range(numNodes + 1)]

    for layer in range(1, numLayers + 1 + 1):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, math.sqrt(clay[layer]), 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)

    state['thermCondPar1'] = thermCondPar1
    state['thermCondPar2'] = thermCondPar2
    state['thermCondPar3'] = thermCondPar3
    state['thermCondPar4'] = thermCondPar4
    return state

def doVolumetricSpecificHeat(state):
    numNodes = state['numNodes']
    volSpecHeatSoil_ = [0.0 for _ in range(numNodes + 1)]
    for node in range(1, numNodes + 1):
        volSpecHeatSoil_[node] = 0.0
        for constituentName in [n for n in state['soilConstituentNames'] if n != "Minerals"]:
            volSpecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node, state) * 1000000.0 * state['soilWater'][node]
    volSpecHeatSoil = [0.0 for _ in range(numNodes + 1)]
    volSpecHeatSoil = mapLayer2Node(volSpecHeatSoil_, state, volSpecHeatSoil)
    state['volSpecHeatSoil'] = volSpecHeatSoil
    return state

def doThermalConductivity(state):
    numNodes = state['numNodes']
    thermCondLayers = [0.0 for _ in range(numNodes + 1)]
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in state['soilConstituentNames']:
            shapeFactorConstituent = shapeFactor(constituentName, node, state)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, state)
            thermalConductanceWater = ThermalConductance("Water", node, state)
            denom1 = 1.0 + shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0)
            denom2 = 1.0 + shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0) * (1.0 - 2.0 * shapeFactorConstituent)
            k = (2.0 / 3.0) * (denom1 ** -1.0) + (1.0 / 3.0) * (denom2 ** -1.0)
            numerator += thermalConductanceConstituent * state['soilWater'][node] * k
            denominator += state['soilWater'][node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = [0.0 for _ in range(numNodes + 1)]
    thermalConductivity = mapLayer2Node(thermCondLayers, state, thermalConductivity)
    state['thermalConductivity'] = thermalConductivity
    return state

def doThomas(state, newTemps, weather_AirPressure, weather_Radn, waterBalance_Eos, waterBalance_Es):
    airNode = state['airNode']
    surfaceNode = state['surfaceNode']
    numNodes = state['numNodes']
    thermalConductivity = state['thermalConductivity']
    volSpecHeatSoil = state['volSpecHeatSoil']
    soilTemp = state['soilTemp']
    nu = state['nu']
    timestep = state['timestep']
    netRadiationSource = state['netRadiationSource']
    netRadiation = state['netRadiation']
    internalTimeStep = state['internalTimeStep']
    thermalConductance = state['thermalConductance']
    heatStorage = state['heatStorage']
    nodeDepth = state['nodeDepth']
    waterBalance_Salb = state['waterBalance_Salb']

    a = [0.0 for _ in range(numNodes + 2)]
    b = [0.0 for _ in range(numNodes + 1)]
    c = [0.0 for _ in range(numNodes + 1)]
    d = [0.0 for _ in range(numNodes + 1)]

    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (state['nodeDepth'][node + 1] - state['nodeDepth'][node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1.0 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + g * thermalConductance[node] * soilTemp[node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * state['latentHeatOfVapourisation'], timestep, 0.0)
    else:
        radnNet = 0.0

    latentHeatFlux = Divide(waterBalance_Es * state['latentHeatOfVapourisation'], timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] += soilSurfaceHeatFlux

    d[numNodes] += nu * thermalConductance[numNodes] * newTemps[numNodes + 1]

    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] -= a[node + 1] * c[node]
        d[node + 1] -= a[node + 1] * d[node]

    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)
    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - c[node] * newTemps[node + 1]
        boundCheck(newTemps[node], -50.0, 100.0, f"newTemps({node})")

    state['thermalConductance'] = thermalConductance
    state['heatStorage'] = heatStorage
    return newTemps, state

def interpolateTemperature(state, timeHours, weather_MeanT, weather_MaxT, weather_MinT):
    defaultTimeOfMaximumTemperature = state['defaultTimeOfMaximumTemperature']
    maxTempYesterday = state['maxTempYesterday']
    minTempYesterday = state['minTempYesterday']

    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5

    if time < minT_time:
        midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (maxTempYesterday - minTempYesterday) / 2.0 + (maxTempYesterday + minTempYesterday) / 2.0
        tScale = Divide((minT_time - time), minT_time, 0.0)
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature

def doUpdate(state, numInterationsPerDay):
    surfaceNode = state['surfaceNode']
    numNodes = state['numNodes']
    internalTimeStep = state['internalTimeStep']

    state['soilTemp'] = state['newTemperature'][:]

    if state['timeOfDaySecs'] < internalTimeStep * 1.2:
        for node in range(surfaceNode, numNodes + 1):
            state['minSoilTemp'][node] = state['soilTemp'][node]
            state['maxSoilTemp'][node] = state['soilTemp'][node]

    for node in range(surfaceNode, numNodes + 1):
        if state['soilTemp'][node] < state['minSoilTemp'][node]:
            state['minSoilTemp'][node] = state['soilTemp'][node]
        elif state['soilTemp'][node] > state['maxSoilTemp'][node]:
            state['maxSoilTemp'][node] = state['soilTemp'][node]
        state['aveSoilTemp'][node] += Divide(state['soilTemp'][node], numInterationsPerDay, 0.0)

    state['boundaryLayerConductance'] += Divide(state['thermalConductivity'][state['airNode']], numInterationsPerDay, 0.0)
    return state

def getBoundaryLayerConductance(state, TNew_zb, weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98

    SpecificHeatAir = specificHeatOfAir * airDensity(state['airTemperature'], weather_AirPressure)
    roughnessFactorMomentum = 0.13 * state['canopyHeight']
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * state['canopyHeight']

    surfaceTemperature = TNew_zb[state['surfaceNode']]

    diffusePenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    radiativeConductance = 4.0 * state['stefanBoltzmannConstant'] * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(state['airTemperature']) ** 3)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for _ in range(1, 3 + 1):
        log_mom = math.log(Divide((state['instrumentHeight'] - d + roughnessFactorMomentum), roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, log_mom, 0.0)
        log_heat = math.log(Divide((state['instrumentHeight'] - d + roughnessFactorHeat), roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, log_heat, 0.0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - state['airTemperature'])
        denom = SpecificHeatAir * kelvinT(state['airTemperature']) * (frictionVelocity ** 3.0)
        stabilityParammeter = Divide(-vonKarmanConstant * state['instrumentHeight'] * gravitationalConstant * heatFluxDensity, denom, 0.0)

        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            term = (1.0 + math.sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0
            stabilityCorrectionHeat = -2.0 * math.log(term) if term > 0 else 0.0
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond

def calcSoilTemperature(state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, soilTempIO):
    cumulativeDepth = ToCumThickness(state['thickness'])
    w = 2.0 * math.pi / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = math.sqrt(2.0 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTemp_local = [0.0 for _ in range(state['numNodes'] + 2)]
    for nodes in range(1, state['numNodes'] + 1):
        soilTemp_local[nodes] = weather_Tav + weather_Amp * math.exp(-1.0 * cumulativeDepth[nodes] / zd) * math.sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)
    for idx in range(state['surfaceNode'], state['surfaceNode'] + state['numNodes']):
        if idx < len(soilTempIO) and idx < len(soilTemp_local):
            soilTempIO[idx] = soilTemp_local[idx - (state['surfaceNode'] - 0)]
    return soilTempIO

def calcLayerTemperature(state, depthLag, alx, deltaTemp, weather_Tav, weather_Amp):
    return weather_Tav + (weather_Amp / 2.0 * math.cos(alx - depthLag) + deltaTemp) * math.exp(-depthLag)

def calcSurfaceTemperature(state, weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    boundCheck(surfaceT, -100.0, 100.0, "Initial surfaceT")
    return surfaceT

def doNetRadiation(state, ITERATIONSperDAY, weather_Latitude, weather_Radn, weather_MinT, clock_Today_DayOfYear):
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin(4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + 0.03345 * math.sin(6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))
    cD = math.sqrt(1.0 - solarDeclination * solarDeclination)
    m1 = [0.0 for _ in range(ITERATIONSperDAY + 1)]
    m1Tot = 0.0
    lat_rad = weather_Latitude * math.pi / 180.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(lat_rad) + cD * math.cos(lat_rad) * math.cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0

    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)

    solarRadn = [0.0 for _ in range(ITERATIONSperDAY + 1)]
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)

    cva = math.exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva

def interpolateNetRadiation(state, solarRadn, cloudFr, cva, airTemperature, waterBalance_Eo, waterBalance_Eos, waterBalance_Salb):
    surfaceEmissivity = 0.96
    internalTimeStep = state['internalTimeStep']
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, state['soilTemp'][state['surfaceNode']]) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil

def getIniVariables(state, weather_Tav, instrumHeight, defaultInstrumentHeight):
    boundCheck(weather_Tav, -30.0, 50.0, "tav (oC)")
    if instrumHeight > 0.00001:
        state['instrumentHeight'] = instrumHeight
    else:
        state['instrumentHeight'] = defaultInstrumentHeight
    return state

def getProfileVariables(state, physical_Thickness, physical_BD, waterBalance_SW, organic_Carbon, physical_Rocks, physical_ParticleSizeSand, physical_ParticleSizeSilt, physical_ParticleSizeClay, DepthToConstantTemperature):
    airNode = state['airNode']
    surfaceNode = state['surfaceNode']
    topsoilNode = state['topsoilNode']
    numPhantomNodes = state['numPhantomNodes']

    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes
    state['numLayers'] = numLayers
    state['numNodes'] = numNodes

    thickness = [0.0 for _ in range(numLayers + numPhantomNodes + 1)]
    for i in range(numLayers):
        thickness[i + 1] = physical_Thickness[i]
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes
    state['thickness'] = thickness

    nodeDepth = [0.0 for _ in range(numNodes + 2)]
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0
    state['nodeDepth'] = nodeDepth

    bulkDensity = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for i in range(numLayers):
        bulkDensity[i + 1] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]
    state['bulkDensity'] = bulkDensity

    soilWater = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            # replicate original odd indexing conversion
            th_prev = thickness[layer - 1] if (layer - 1) < len(thickness) else 0.0
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * th_prev, thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]
    state['soilWater'] = soilWater

    carbon = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]
    state['carbon'] = carbon

    rocks = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]
    state['rocks'] = rocks

    sand = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]
    state['sand'] = sand

    silt = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]
    state['silt'] = silt

    clay = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        clay[layer] = clay[numLayers]
    state['clay'] = clay

    state['maxSoilTemp'] = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    state['minSoilTemp'] = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    state['aveSoilTemp'] = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    state['volSpecHeatSoil'] = [0.0 for _ in range(numNodes + 1)]
    state['soilTemp'] = [0.0 for _ in range(numNodes + 2)]
    state['morningSoilTemp'] = [0.0 for _ in range(numNodes + 2)]
    state['newTemperature'] = [0.0 for _ in range(numNodes + 2)]
    state['thermalConductivity'] = [0.0 for _ in range(numNodes + 1)]
    state['heatStorage'] = [0.0 for _ in range(numNodes + 1)]
    state['thermalConductance'] = [0.0 for _ in range(numNodes + 2)]

    return state

def readParam(state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear):
    state = doThermalConductivityCoeffs(state)
    state['soilTemp'] = calcSoilTemperature(state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, state['soilTemp'])
    state['newTemperature'] = state['soilTemp'][:]
    state['soilRoughnessHeight'] = state['bareSoilRoughness']
    return state

def getOtherVariables(state, waterBalance_SW, microClimate_CanopyHeight):
    numLayers = state['numLayers']
    numNodes = state['numNodes']
    if waterBalance_SW is not None:
        for i in range(numLayers):
            state['soilWater'][i + 1] = waterBalance_SW[i]
    state['soilWater'][numNodes] = state['soilWater'][state['numLayers']]
    canopyHeight_mm = microClimate_CanopyHeight if microClimate_CanopyHeight is not None else 0.0
    state['canopyHeight'] = max(canopyHeight_mm, state['soilRoughnessHeight']) / 1000.0
    state['instrumentHeight'] = max(state['instrumentHeight'], state['canopyHeight'] + 0.5)
    return state

def doProcess(state, weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp, weather_AirPressure, weather_Radn, weather_Wind, weather_Latitude, waterBalance_Eo, waterBalance_Eos, waterBalance_Es, waterBalance_Salb, microClimate_CanopyHeight, clock_Today_DayOfYear):
    interactionsPerDay = 48
    solarRadn, cloudFr, cva = doNetRadiation(state, interactionsPerDAY=interactionsPerDay, weather_Latitude=weather_Latitude, weather_Radn=weather_Radn, weather_MinT=weather_MinT, clock_Today_DayOfYear=clock_Today_DayOfYear)

    Zero(state['minSoilTemp'])
    Zero(state['maxSoilTemp'])
    Zero(state['aveSoilTemp'])
    state['boundaryLayerConductance'] = 0.0

    state['internalTimeStep'] = round(state['timestep'] / interactionsPerDay)

    state = doVolumetricSpecificHeat(state)
    state = doThermalConductivity(state)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        state['timeOfDaySecs'] = state['internalTimeStep'] * float(timeStepIteration)
        if state['timestep'] < 24.0 * 60.0 * 60.0:
            state['airTemperature'] = weather_MeanT
        else:
            state['airTemperature'] = interpolateTemperature(state, (state['timeOfDaySecs'] / 3600.0), weather_MeanT, weather_MaxT, weather_MinT)
        state['newTemperature'][state['airNode']] = state['airTemperature']

        state['netRadiation'] = interpolateNetRadiation(state, solarRadn[timeStepIteration], cloudFr, cva, state['airTemperature'], waterBalance_Eo, waterBalance_Eos, waterBalance_Salb)

        if state['boundarLayerConductanceSource'] == "constant":
            state['thermalConductivity'][state['airNode']] = state['constantBoundaryLayerConductance']
        elif state['boundarLayerConductanceSource'] == "calc":
            state['thermalConductivity'][state['airNode']] = getBoundaryLayerConductance(state, state['newTemperature'], weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo)
            for _ in range(1, state['numIterationsForBoundaryLayerConductance'] + 1):
                state['newTemperature'], state = doThomas(state, state['newTemperature'], weather_AirPressure, weather_Radn, waterBalance_Eos, waterBalance_Es)
                state['thermalConductivity'][state['airNode']] = getBoundaryLayerConductance(state, state['newTemperature'], weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo)
        state['newTemperature'], state = doThomas(state, state['newTemperature'], weather_AirPressure, weather_Radn, waterBalance_Eos, waterBalance_Es)
        state = doUpdate(state, interactionsPerDay)

        if abs(state['timeOfDaySecs'] - 5.0 * 3600.0) <= min(state['timeOfDaySecs'], 5.0 * 3600.0) * 0.0001:
            state['morningSoilTemp'] = state['soilTemp'][:]

    state['minTempYesterday'] = weather_MinT
    state['maxTempYesterday'] = weather_MaxT
    return state

def SoilTemperature_initialize(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    waterBalance_Eo,
    waterBalance_Eos,
    waterBalance_Es,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_AirPressure,
    weather_Radn,
    weather_Wind,
    DepthToConstantTemperature=10000.0,
    instrumHeight=0.0,
    defaultTimeOfMaximumTemperature=14.0,
    InitialValues=None,
    clock_Today_DayOfYear=1
):
    state = {}
    state['pom'] = 1.3
    state['ps'] = 2.63
    state['soilConstituentNames'] = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    state['timestep'] = 24.0 * 60.0 * 60.0
    state['latentHeatOfVapourisation'] = 2465000.0
    state['stefanBoltzmannConstant'] = 0.0000000567
    state['airNode'] = 0
    state['surfaceNode'] = 1
    state['topsoilNode'] = 2
    state['numPhantomNodes'] = 5
    state['constantBoundaryLayerConductance'] = 20.0
    state['numIterationsForBoundaryLayerConductance'] = 1
    state['defaultTimeOfMaximumTemperature'] = defaultTimeOfMaximumTemperature
    state['defaultInstrumentHeight'] = 1.2
    state['bareSoilRoughness'] = 57.0
    state['doInitialisationStuff'] = True
    state['internalTimeStep'] = 0.0
    state['timeOfDaySecs'] = 0.0
    state['numNodes'] = 0
    state['numLayers'] = 0
    state['nodeDepth'] = []
    state['thermCondPar1'] = []
    state['thermCondPar2'] = []
    state['thermCondPar3'] = []
    state['thermCondPar4'] = []
    state['volSpecHeatSoil'] = []
    state['soilTemp'] = []
    state['morningSoilTemp'] = []
    state['heatStorage'] = []
    state['thermalConductance'] = []
    state['thermalConductivity'] = []
    state['boundaryLayerConductance'] = 0.0
    state['newTemperature'] = []
    state['airTemperature'] = 0.0
    state['maxTempYesterday'] = 0.0
    state['minTempYesterday'] = 0.0
    state['soilWater'] = []
    state['minSoilTemp'] = []
    state['maxSoilTemp'] = []
    state['aveSoilTemp'] = []
    state['thickness'] = []
    state['bulkDensity'] = []
    state['rocks'] = []
    state['carbon'] = []
    state['sand'] = []
    state['silt'] = []
    state['clay'] = []
    state['soilRoughnessHeight'] = 0.0
    state['instrumentHeight'] = 0.0
    state['netRadiation'] = 0.0
    state['canopyHeight'] = 0.0
    state['instrumHeight'] = instrumHeight
    state['nu'] = 0.6
    state['boundarLayerConductanceSource'] = "calc"
    state['netRadiationSource'] = "calc"
    state['InitialValues'] = InitialValues[:] if InitialValues is not None else None
    state['waterBalance_Salb'] = waterBalance_Salb

    state = getIniVariables(state, weather_Tav, instrumHeight, state['defaultInstrumentHeight'])
    state = getProfileVariables(
        state,
        physical_Thickness=physical_Thickness,
        physical_BD=physical_BD,
        waterBalance_SW=waterBalance_SW,
        organic_Carbon=organic_Carbon,
        physical_Rocks=physical_Rocks,
        physical_ParticleSizeSand=physical_ParticleSizeSand,
        physical_ParticleSizeSilt=physical_ParticleSizeSilt,
        physical_ParticleSizeClay=physical_ParticleSizeClay,
        DepthToConstantTemperature=DepthToConstantTemperature
    )
    state = readParam(state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)

    return state

def SoilTemperature_process(
    state,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Tav,
    weather_Amp,
    weather_AirPressure,
    weather_Radn,
    weather_Wind,
    weather_Latitude,
    waterBalance_SW,
    waterBalance_Eo,
    waterBalance_Eos,
    waterBalance_Es,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    clock_Today_DayOfYear
):
    state = getOtherVariables(state, waterBalance_SW, microClimate_CanopyHeight)

    if state['doInitialisationStuff']:
        if ValuesInArray(state['InitialValues']):
            state['soilTemp'] = [0.0 for _ in range(state['numNodes'] + 2)]
            for i in range(len(state['InitialValues'])):
                idx = state['topsoilNode'] + i
                if idx < len(state['soilTemp']):
                    state['soilTemp'][idx] = state['InitialValues'][i]
        else:
            state['soilTemp'] = calcSoilTemperature(state, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, state['soilTemp'])
            state['InitialValues'] = [0.0 for _ in range(state['numLayers'])]
            for i in range(state['numLayers']):
                state['InitialValues'][i] = state['soilTemp'][state['topsoilNode'] + i]

        state['soilTemp'][state['airNode']] = weather_MeanT
        state['soilTemp'][state['surfaceNode']] = calcSurfaceTemperature(state, weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb)

        for i in range(state['numNodes'] + 1, len(state['soilTemp'])):
            state['soilTemp'][i] = weather_Tav

        state['newTemperature'] = state['soilTemp'][:]
        state['maxTempYesterday'] = weather_MaxT
        state['minTempYesterday'] = weather_MinT
        state['doInitialisationStuff'] = False

    state['waterBalance_Salb'] = waterBalance_Salb
    state = doProcess(
        state,
        weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp,
        weather_AirPressure, weather_Radn, weather_Wind, weather_Latitude,
        waterBalance_Eo, waterBalance_Eos, waterBalance_Es, waterBalance_Salb,
        microClimate_CanopyHeight, clock_Today_DayOfYear
    )
    return state

def FinalSoilTemperature(state):
    numLayers = state['numLayers']
    result = [0.0 for _ in range(numLayers)]
    for i in range(numLayers):
        result[i] = state['soilTemp'][state['topsoilNode'] + i]
    return result

def FinalSoilSurfaceTemperature(state):
    return state['soilTemp'][state['surfaceNode']]

def AverageSoilTemperature(state):
    numLayers = state['numLayers']
    result = [0.0 for _ in range(numLayers)]
    for i in range(numLayers):
        result[i] = state['aveSoilTemp'][state['topsoilNode'] + i]
    return result

def AverageSoilSurfaceTemperature(state):
    return state['aveSoilTemp'][state['surfaceNode']]

def MinimumSoilTemperature(state):
    numLayers = state['numLayers']
    result = [0.0 for _ in range(numLayers)]
    for i in range(numLayers):
        result[i] = state['minSoilTemp'][state['topsoilNode'] + i]
    return result

def MinimumSoilSurfaceTemperature(state):
    return state['minSoilTemp'][state['surfaceNode']]

def MaximumSoilTemperature(state):
    numLayers = state['numLayers']
    result = [0.0 for _ in range(numLayers)]
    for i in range(numLayers):
        result[i] = state['maxSoilTemp'][state['topsoilNode'] + i]
    return result

def MaximumSoilSurfaceTemperature(state):
    return state['maxSoilTemp'][state['surfaceNode']]

def BoundaryLayerConductance(state):
    return state['boundaryLayerConductance']

def ThermalConductivity_profile(state):
    numNodes = state['numNodes']
    result = [0.0 for _ in range(numNodes)]
    # Original returned from index 1
    for i in range(numNodes):
        result[i] = state['thermalConductivity'][i + 1] if (i + 1) < len(state['thermalConductivity']) else 0.0
    return result

def HeatCapacity_profile(state):
    numNodes = state['numNodes']
    result = [0.0 for _ in range(numNodes)]
    for i in range(numNodes):
        idx = state['surfaceNode'] + i
        result[i] = state['volSpecHeatSoil'][idx] if idx < len(state['volSpecHeatSoil']) else 0.0
    return result

def HeatStore_profile(state):
    numNodes = state['numNodes']
    result = [0.0 for _ in range(numNodes)]
    for i in range(numNodes):
        idx = state['surfaceNode'] + i
        result[i] = state['heatStorage'][idx] if idx < len(state['heatStorage']) else 0.0
    return result

def Thr_profile(state):
    return state['morningSoilTemp'][:]