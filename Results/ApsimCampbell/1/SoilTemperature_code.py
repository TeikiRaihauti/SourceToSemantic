def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and not (Value != Value):
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


def volumetricSpecificHeat(name, layer, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity):
    specificHeatRocks = 7.7
    specificHeatOM = 0.25
    specificHeatSand = 7.7
    specificHeatSilt = 2.74
    specificHeatClay = 2.92
    specificHeatWater = 0.57
    specificHeatIce = 2.18
    specificHeatAir = 0.025

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
    else:
        result = 0.0

    return result


def ThermalConductance(name, layer, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity):
    thermalConductanceRocks = 0.182
    thermalConductanceOM = 2.50
    thermalConductanceSand = 0.182
    thermalConductanceSilt = 2.39
    thermalConductanceClay = 1.39
    thermalConductanceWater = 4.18
    thermalConductanceIce = 1.73
    thermalConductanceAir = 0.0012

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
        # Placeholder as in original (odd formula commented there); here we don't compute this branch further
        result = thermalConductanceRocks
    else:
        result = 0.0

    # Preserve original erroneous assignment
    result = volumetricSpecificHeat(name, layer, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity)
    return result


def volumetricFractionRocks(layer, rocks):
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom):
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(layer, carbon, rocks, sand, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer, carbon, rocks, silt, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(layer, carbon, rocks, clay, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(layer, carbon, soilWater):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, None, 1.0)) * soilWater[layer]  # bulkDensity not used in this call as per original


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
    return 1.0 - volumetricFractionRocks(layer, rocks) - \
           volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - \
           volumetricFractionSand(layer, carbon, rocks, sand, bulkDensity, pom, ps) - \
           volumetricFractionSilt(layer, carbon, rocks, silt, bulkDensity, pom, ps) - \
           volumetricFractionClay(layer, carbon, rocks, clay, bulkDensity, pom, ps) - \
           volumetricFractionWater(layer, carbon, soilWater) - \
           volumetricFractionIce(layer)


def shapeFactor(name, layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0

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
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / \
                 (volumetricFractionWater(layer, carbon, soilWater) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps) / \
                 (volumetricFractionWater(layer, carbon, soilWater) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + \
                 shapeFactorSand * volumetricFractionSand(layer, carbon, rocks, sand, bulkDensity, pom, ps) + \
                 shapeFactorSilt * volumetricFractionSilt(layer, carbon, rocks, silt, bulkDensity, pom, ps) + \
                 shapeFactorClay * volumetricFractionClay(layer, carbon, rocks, clay, bulkDensity, pom, ps)
    else:
        result = 0.0

    # Preserve original erroneous assignment
    result = volumetricSpecificHeat(name, layer, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity)
    return result


def mapLayer2Node(layerArray, nodeDepth, thickness, numNodes):
    surfaceNode = 1
    nodeArray = [0.0] * (numNodes + 1)
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0) + Divide(layerArray[layer + 1] * d2, dSum, 0)
    return nodeArray


def calcLayerTemperature(depthLag, alx, deltaTemp, weather_Tav, weather_Amp):
    return weather_Tav + (weather_Amp / 2.0 * (math.cos(alx - depthLag)) + deltaTemp) * (math.exp(-depthLag))


def calcSurfaceTemperature(weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * ((max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5)) + waterBalance_Salb * weather_MeanT
    return surfaceT


def doNetRadiation(clock_Today_DayOfYear, ITERATIONSperDAY, weather_Latitude, weather_Radn, weather_MinT):
    TSTEPS2RAD = Divide(2.0 * 3.141592653589793, float(ITERATIONSperDAY), 0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin(4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + 0.03345 * math.sin(6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))
    cD = math.sqrt(1.0 - solarDeclination * solarDeclination)
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * math.pi / 180.0) + cD * math.cos(weather_Latitude * math.pi / 180.0) * math.cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0

    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0)
    cva = math.exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(internalTimeStep, solarRadn, cloudFr, cva, airTemperature, soilSurfaceTemperature, waterBalance_Eo, waterBalance_Eos, weather_Wind, instrumentHeight, canopyHeight):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0

    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0)

    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilSurfaceTemperature) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil

    swRin = solarRadn
    # Albedo applied later in doThomas
    swRnetSoil = (swRin - 0.0) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def interpolateTemperature(timeHours, defaultTimeOfMaximumTemperature, weather_MaxT, weather_MinT, weather_MeanT, maxTempYesterday, minTempYesterday):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5

    if time < minT_time:
        midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (maxTempYesterday - minTempYesterday) / 2.0 + (maxTempYesterday + minTempYesterday) / 2.0
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0:
            tScale = 0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


def doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, nodeDepth, thickness):
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in [n for n in soilConstituentNames if n != "Minerals"]:
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, nodeDepth, thickness, numNodes)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, nodeDepth, thickness):
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity)
            thermalConductanceWater = ThermalConductance("Water", node, pom, ps, rocks, carbon, sand, silt, clay, soilWater, bulkDensity)
            k = (2.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** (-1) + \
                (1.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** (-1)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = mapLayer2Node(thermCondLayers, nodeDepth, thickness, numNodes)
    return thermalConductivity


def doThomas(newTemps, soilTemp, volSpecHeatSoil, thermalConductivity, internalTimeStep, netRadiation, waterBalance_Salb, waterBalance_Es, waterBalance_Eos, timestep_seconds, numNodes, nodeDepth, nu):
    airNode = 0
    surfaceNode = 1

    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    thermalConductance = [0.0] * (numNodes + 2)
    heatStorage = [0.0] * (numNodes + 1)

    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0)

    g = 1 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + g * thermalConductance[node] * soilTemp[node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0)

    latentHeatOfVapourisation = 2465000.0
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep_seconds, 0)

    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] += soilSurfaceHeatFlux

    d[numNodes] += nu * thermalConductance[numNodes] * newTemps[numNodes + 1]

    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0)
        d[node] = Divide(d[node], b[node], 0)
        b[node + 1] -= a[node + 1] * c[node]
        d[node + 1] -= a[node + 1] * d[node]
    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0)

    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - c[node] * newTemps[node + 1]
    return newTemps, thermalConductance, heatStorage


def getBoundaryLayerConductance(TNew_zb, weather_Wind, airTemperature, weather_AirPressure, canopyHeight_m, instrumentHeight, waterBalance_Eo, waterBalance_Eos):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98

    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)

    roughnessFactorMomentum = 0.13 * canopyHeight_m
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight_m

    surfaceNode = 1
    surfaceTemperature = TNew_zb[surfaceNode]

    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    stefanBoltzmannConstant = 0.0000000567
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for iteration in range(1, 4):
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0)) + stabilityCorrectionMomentum, 0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0)) + stabilityCorrectionHeat, 0)
        boundaryLayerCond += radiativeConductance

        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) ** 3 if False else SpecificHeatAir * (frictionVelocity ** 3), 0)

        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + math.sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond


def calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25

    soilTemp_local = [0.0] * (numNodes + 2)
    for nodes in range(1, numNodes + 1):
        soilTemp_local[nodes] = weather_Tav + weather_Amp * math.exp(-1 * cumulativeDepth[nodes] / zd) * math.sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)
    return soilTemp_local


def OnStartOfSimulation(
    weather_MeanT,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_Radn,
    weather_AirPressure,
    weather_MaxT,
    weather_MinT,
    weather_Wind,
    waterBalance_SW,
    waterBalance_Eo,
    waterBalance_Eos,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    clock_Today_DayOfYear,
    instrumHeight=0.0,
    InitialValues=None,
    DepthToConstantTemperature=10000.0
):
    import math as _math
    globals()['_math'] = _math  # ensure math is available to other functions
    global math
    math = _math

    # Constants and indices
    timestep_seconds = 24.0 * 60.0 * 60.0
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    defaultInstrumentHeight = 1.2
    bareSoilRoughness = 57
    defaultTimeOfMaximumTemperature = 14.0

    # Instrument height
    instrumentHeight = instrumHeight if instrumHeight > 0.00001 else defaultInstrumentHeight

    # Determine layers and nodes
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    # Thickness array with phantom layers and surface slot
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
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * (physical_Thickness[layer - 1] if layer - 1 < len(physical_Thickness) else 0.0), thickness[layer], 0)
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

    # Parameters for thermal conductivity, computed but not used further (preserved behavior)
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 2):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5 if clay[layer] > 0 else 1.0), 0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)

    # Initial temperatures
    if ValuesInArray(InitialValues):
        soilTemp = [0.0] * (numNodes + 2)
        for i in range(len(InitialValues)):
            soilTemp[topsoilNode + i] = InitialValues[i]
    else:
        soilTemp_calc = calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)
        # Copy into soilTemp starting at surfaceNode
        for i in range(numNodes + 2):
            if i < len(soilTemp_calc):
                if surfaceNode + i < len(soilTemp):
                    soilTemp[surfaceNode + i] = soilTemp_calc[i]
        InitialValues = [soilTemp[topsoilNode + i] for i in range(numLayers)]

    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, weather_Radn, waterBalance_Salb)

    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav

    newTemperature = list(soilTemp)
    maxTempYesterday = weather_MaxT
    minTempYesterday = weather_MinT

    # Fixed parameters and settings
    soilRoughnessHeight = bareSoilRoughness
    boundaryLayerConductance = 0.0
    nu = 0.6
    boundarLayerConductanceSource = "calc"
    netRadiationSource = "calc"

    state = {
        "timestep_seconds": timestep_seconds,
        "airNode": airNode,
        "surfaceNode": surfaceNode,
        "topsoilNode": topsoilNode,
        "numPhantomNodes": numPhantomNodes,
        "defaultTimeOfMaximumTemperature": defaultTimeOfMaximumTemperature,
        "instrumentHeight": instrumentHeight,
        "numLayers": numLayers,
        "numNodes": numNodes,
        "thickness": thickness,
        "nodeDepth": nodeDepth,
        "bulkDensity": bulkDensity,
        "soilWater": soilWater,
        "carbon": carbon,
        "rocks": rocks,
        "sand": sand,
        "silt": silt,
        "clay": clay,
        "maxSoilTemp": maxSoilTemp,
        "minSoilTemp": minSoilTemp,
        "aveSoilTemp": aveSoilTemp,
        "volSpecHeatSoil": volSpecHeatSoil,
        "soilTemp": soilTemp,
        "morningSoilTemp": morningSoilTemp,
        "newTemperature": newTemperature,
        "thermalConductivity": thermalConductivity,
        "heatStorage": heatStorage,
        "thermalConductance": thermalConductance,
        "thermCondPar1": thermCondPar1,
        "thermCondPar2": thermCondPar2,
        "thermCondPar3": thermCondPar3,
        "thermCondPar4": thermCondPar4,
        "soilRoughnessHeight": soilRoughnessHeight,
        "boundaryLayerConductance": boundaryLayerConductance,
        "airTemperature": 0.0,
        "maxTempYesterday": maxTempYesterday,
        "minTempYesterday": minTempYesterday,
        "netRadiation": 0.0,
        "canopyHeight": 0.0,
        "instrumHeight": instrumHeight,
        "nu": nu,
        "boundarLayerConductanceSource": boundarLayerConductanceSource,
        "netRadiationSource": netRadiationSource,
        "InitialValues": InitialValues,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "pom": 1.3,
        "ps": 2.63,
        "soilConstituentNames": ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"],
    }
    return state


def OnProcess(
    state,
    weather_MeanT,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_Radn,
    weather_AirPressure,
    weather_MaxT,
    weather_MinT,
    weather_Wind,
    waterBalance_SW,
    waterBalance_Eo,
    waterBalance_Eos,
    waterBalance_Salb,
    waterBalance_Es,
    microClimate_CanopyHeight,
    clock_Today_DayOfYear
):
    # Unpack state
    timestep_seconds = state["timestep_seconds"]
    airNode = state["airNode"]
    surfaceNode = state["surfaceNode"]
    topsoilNode = state["topsoilNode"]
    numPhantomNodes = state["numPhantomNodes"]
    defaultTimeOfMaximumTemperature = state["defaultTimeOfMaximumTemperature"]
    instrumentHeight = state["instrumentHeight"]
    numLayers = state["numLayers"]
    numNodes = state["numNodes"]
    thickness = state["thickness"]
    nodeDepth = state["nodeDepth"]
    bulkDensity = state["bulkDensity"]
    soilWater = state["soilWater"]
    carbon = state["carbon"]
    rocks = state["rocks"]
    sand = state["sand"]
    silt = state["silt"]
    clay = state["clay"]
    maxSoilTemp = state["maxSoilTemp"]
    minSoilTemp = state["minSoilTemp"]
    aveSoilTemp = state["aveSoilTemp"]
    volSpecHeatSoil = state["volSpecHeatSoil"]
    soilTemp = state["soilTemp"]
    morningSoilTemp = state["morningSoilTemp"]
    newTemperature = state["newTemperature"]
    thermalConductivity = state["thermalConductivity"]
    heatStorage = state["heatStorage"]
    thermalConductance = state["thermalConductance"]
    soilRoughnessHeight = state["soilRoughnessHeight"]
    boundaryLayerConductance = 0.0
    airTemperature = state["airTemperature"]
    maxTempYesterday = state["maxTempYesterday"]
    minTempYesterday = state["minTempYesterday"]
    netRadiation = state["netRadiation"]
    nu = state["nu"]
    boundarLayerConductanceSource = state["boundarLayerConductanceSource"]
    netRadiationSource = state["netRadiationSource"]
    DepthToConstantTemperature = state["DepthToConstantTemperature"]
    pom = state["pom"]
    ps = state["ps"]
    soilConstituentNames = state["soilConstituentNames"]

    # getOtherVariables
    for i in range(1, min(len(soilWater), len(waterBalance_SW) + 1)):
        soilWater[i] = waterBalance_SW[i - 1]
    soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)

    # Process a day
    interactionsPerDay = 48
    solarRadn, cloudFr, cva = doNetRadiation(clock_Today_DayOfYear, interactionsPerDay, weather_Latitude, weather_Radn, weather_MinT)

    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0

    internalTimeStep = round(timestep_seconds / interactionsPerDay)

    volSpecHeatSoil = doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, nodeDepth, thickness)
    thermalConductivity = doThermalConductivity(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, nodeDepth, thickness)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep_seconds < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0, defaultTimeOfMaximumTemperature, weather_MaxT, weather_MinT, weather_MeanT, maxTempYesterday, minTempYesterday)
        newTemperature[airNode] = airTemperature

        netRadiation = interpolateNetRadiation(internalTimeStep, solarRadn[timeStepIteration], cloudFr, cva, airTemperature, soilTemp[surfaceNode], waterBalance_Eo, waterBalance_Eos, weather_Wind, instrumentHeight, canopyHeight)

        if boundarLayerConductanceSource == "constant":
            constantBoundaryLayerConductance = 20
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_Wind, airTemperature, weather_AirPressure, canopyHeight, instrumentHeight, waterBalance_Eo, waterBalance_Eos)
            numIterationsForBoundaryLayerConductance = 1
            for iteration in range(1, numIterationsForBoundaryLayerConductance + 1):
                newTemperature, thermalConductance, heatStorage = doThomas(newTemperature, soilTemp, volSpecHeatSoil, thermalConductivity, internalTimeStep, netRadiation, waterBalance_Salb, waterBalance_Es, waterBalance_Eos, timestep_seconds, numNodes, nodeDepth, nu)
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_Wind, airTemperature, weather_AirPressure, canopyHeight, instrumentHeight, waterBalance_Eo, waterBalance_Eos)

        newTemperature, thermalConductance, heatStorage = doThomas(newTemperature, soilTemp, volSpecHeatSoil, thermalConductivity, internalTimeStep, netRadiation, waterBalance_Salb, waterBalance_Es, waterBalance_Eos, timestep_seconds, numNodes, nodeDepth, nu)

        # doUpdate
        soilTemp = list(newTemperature)
        if timeOfDaySecs < internalTimeStep * 1.2:
            for node in range(surfaceNode, numNodes + 1):
                minSoilTemp[node] = soilTemp[node]
                maxSoilTemp[node] = soilTemp[node]
        for node in range(surfaceNode, numNodes + 1):
            if soilTemp[node] < minSoilTemp[node]:
                minSoilTemp[node] = soilTemp[node]
            elif soilTemp[node] > maxSoilTemp[node]:
                maxSoilTemp[node] = soilTemp[node]
            aveSoilTemp[node] += Divide(soilTemp[node], interactionsPerDay, 0)
        boundaryLayerConductance += Divide(thermalConductivity[airNode], interactionsPerDay, 0)

        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            morningSoilTemp = list(soilTemp)

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    # Pack state
    state["instrumentHeight"] = instrumentHeight
    state["soilWater"] = soilWater
    state["canopyHeight"] = canopyHeight
    state["volSpecHeatSoil"] = volSpecHeatSoil
    state["thermalConductivity"] = thermalConductivity
    state["thermalConductance"] = thermalConductance
    state["heatStorage"] = heatStorage
    state["soilTemp"] = soilTemp
    state["newTemperature"] = newTemperature
    state["minSoilTemp"] = minSoilTemp
    state["maxSoilTemp"] = maxSoilTemp
    state["aveSoilTemp"] = aveSoilTemp
    state["boundaryLayerConductance"] = boundaryLayerConductance
    state["airTemperature"] = airTemperature
    state["netRadiation"] = netRadiation
    state["morningSoilTemp"] = morningSoilTemp
    state["maxTempYesterday"] = maxTempYesterday
    state["minTempYesterday"] = minTempYesterday

    return state

