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


def volumetricSpecificHeat(name, layer):
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


def volumetricFractionRocks(layer, rocks):
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom):
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps):
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps):
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps):
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(layer, carbon, soilWater, bulkDensity, pom):
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom)) * soilWater[layer]


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps) - volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps) - volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps) - volumetricFractionWater(layer, carbon, soilWater, bulkDensity, pom) - volumetricFractionIce(layer)


def ThermalConductance(name, layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
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
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, rocks)) * (thermalConductanceSand ** volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps)) + (thermalConductanceSilt ** volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps)) + (thermalConductanceClay ** volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps))
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(name, layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
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
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / (volumetricFractionWater(layer, carbon, soilWater, bulkDensity, pom) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps) / (volumetricFractionWater(layer, carbon, soilWater, bulkDensity, pom) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + shapeFactorSand * volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps) + shapeFactorSilt * volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps) + shapeFactorClay * volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps)
    result = volumetricSpecificHeat(name, layer)
    return result


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def kelvinT(celciusT):
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def longWaveRadn(emissivity, tDegC):
    return 0.0000000567 * emissivity * (kelvinT(tDegC) ** 4)


def mapLayer2Node(layerArray, nodeArray, surfaceNode, numNodes, nodeDepth, thickness):
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)


def doVolumetricSpecificHeat(numNodes, soilWater, surfaceNode, nodeDepth, thickness):
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName != "Minerals":
                volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, surfaceNode, numNodes, nodeDepth, thickness)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps, surfaceNode, nodeDepth, thickness):
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
            k_val = (2.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** -1) + (1.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** -1)
            numerator += thermalConductanceConstituent * soilWater[node] * k_val
            denominator += soilWater[node] * k_val
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = [0.0] * (numNodes + 1)
    mapLayer2Node(thermCondLayers, thermalConductivity, surfaceNode, numNodes, nodeDepth, thickness)
    return thermalConductivity


def doThomas(newTemps, volSpecHeatSoil, thermalConductivity, soilTemp, surfaceNode, numNodes, nodeDepth, internalTimeStep, nu, netRadiation, netRadiationSource, waterBalance_Eos, waterBalance_Es, timestep, latentHeatOfVapourisation):
    airNode = 0
    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 2)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
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
    radnNet = 0.0
    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
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
    return newTemps, heatStorage, thermalConductance


def interpolateTemperature(timeHours, defaultTimeOfMaximumTemperature, maxTempYesterday, minTempYesterday, weather_MaxT, weather_MinT, weather_MeanT):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = (math_sin((0.0 + 0.25 - maxT_time) * 2.0 * math_pi()) * (maxTempYesterday - minTempYesterday) / 2.0) + ((maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = math_sin((time + 0.25 - maxT_time) * 2.0 * math_pi()) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


def doUpdate(numInterationsPerDay, newTemperature, soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, timeOfDaySecs, internalTimeStep, boundaryLayerConductance, thermalConductivity_airNode, surfaceNode, numNodes):
    for i in range(len(soilTemp)):
        soilTemp[i] = newTemperature[i]
    if timeOfDaySecs < internalTimeStep * 1.2:
        for node in range(surfaceNode, numNodes + 1):
            minSoilTemp[node] = soilTemp[node]
            maxSoilTemp[node] = soilTemp[node]
    for node in range(surfaceNode, numNodes + 1):
        if soilTemp[node] < minSoilTemp[node]:
            minSoilTemp[node] = soilTemp[node]
        elif soilTemp[node] > maxSoilTemp[node]:
            maxSoilTemp[node] = soilTemp[node]
        aveSoilTemp[node] += Divide(soilTemp[node], numInterationsPerDay, 0.0)
    boundaryLayerConductance += Divide(thermalConductivity_airNode, numInterationsPerDay, 0.0)
    return soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def getBoundaryLayerConductance(TNew_zb, surfaceNode, canopyHeight, instrumentHeight, weather_Wind, airTemperature, weather_AirPressure, waterBalance_Eos, waterBalance_Eo):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight
    surfaceTemperature = TNew_zb[surfaceNode]
    diffusePenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    radiativeConductance = 4.0 * 0.0000000567 * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)
    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0
    for _ in range(1, 3 + 1):
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, math_log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum, 0.0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, math_log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat, 0.0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + math_sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


def calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * math_pi() / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = math_sqrt(2.0 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTemp_local = [0.0] * (numNodes + 2)
    for nodes in range(1, numNodes + 1):
        soilTemp_local[nodes] = weather_Tav + weather_Amp * math_exp(-1.0 * cumulativeDepth[nodes] / zd) * math_sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math_pi() - cumulativeDepth[nodes] / zd)
    return soilTemp_local


def calcLayerTemperature(depthLag, alx, deltaTemp, weather_Tav, weather_Amp):
    return weather_Tav + (weather_Amp / 2.0 * math_cos(alx - depthLag) + deltaTemp) * math_exp(-depthLag)


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math_sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    return surfaceT


def doNetRadiation(ITERATIONSperDAY, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT):
    TSTEPS2RAD = Divide(2.0 * math_pi(), float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin(4.869 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25) + 0.03345 * math_sin(6.224 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25)))
    cD = math_sqrt(1.0 - solarDeclination * solarDeclination)
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    latRad = weather_Latitude * math_pi() / 180.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(latRad) + cD * math_cos(latRad) * math_cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0
    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = math_exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(solarRadn_val, cloudFr, cva, internalTimeStep, waterBalance_Eos, waterBalance_Eo, waterBalance_Salb, soilSurfaceTemp, airTemperature):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilSurfaceTemp) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn_val
    swRout = waterBalance_Salb * solarRadn_val
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def soil_temperature_init(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    waterBalance_SW,
    organic_Carbon,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    clock_Today_DayOfYear,
    InitialValues=None,
    DepthToConstantTemperature=10000.0,
    instrumHeight=0.0
):
    import math as _math  # local import to avoid module-level state
    _ = _math  # keep reference
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    defaultInstrumentHeight = 1.2
    bareSoilRoughness = 57.0
    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
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
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0.0)
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
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 2):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5) if clay[layer] > 0 else 0.0, 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    soilTemp_init = calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude)
    for i in range(len(soilTemp)):
        if i < len(soilTemp_init):
            soilTemp[i] = soilTemp_init[i]
    for i in range(len(newTemperature)):
        newTemperature[i] = soilTemp[i]
    soilRoughnessHeight = bareSoilRoughness
    instrumentHeight = instrumHeight if (instrumHeight > 0.00001) else defaultInstrumentHeight
    state = {
        "airNode": airNode,
        "surfaceNode": surfaceNode,
        "topsoilNode": topsoilNode,
        "numPhantomNodes": numPhantomNodes,
        "doInitialisationStuff": doInitialisationStuff,
        "internalTimeStep": internalTimeStep,
        "timeOfDaySecs": timeOfDaySecs,
        "numLayers": numLayers,
        "numNodes": numNodes,
        "nodeDepth": nodeDepth,
        "thermCondPar1": thermCondPar1,
        "thermCondPar2": thermCondPar2,
        "thermCondPar3": thermCondPar3,
        "thermCondPar4": thermCondPar4,
        "volSpecHeatSoil": volSpecHeatSoil,
        "soilTemp": soilTemp,
        "morningSoilTemp": morningSoilTemp,
        "heatStorage": heatStorage,
        "thermalConductance": thermalConductance,
        "thermalConductivity": thermalConductivity,
        "boundaryLayerConductance": 0.0,
        "newTemperature": newTemperature,
        "airTemperature": 0.0,
        "maxTempYesterday": 0.0,
        "minTempYesterday": 0.0,
        "soilWater": soilWater,
        "minSoilTemp": minSoilTemp,
        "maxSoilTemp": maxSoilTemp,
        "aveSoilTemp": aveSoilTemp,
        "thickness": thickness,
        "bulkDensity": bulkDensity,
        "rocks": rocks,
        "carbon": carbon,
        "sand": sand,
        "silt": silt,
        "clay": clay,
        "soilRoughnessHeight": soilRoughnessHeight,
        "instrumentHeight": instrumentHeight,
        "netRadiation": 0.0,
        "canopyHeight": 0.0,
        "instrumHeight": instrumHeight,
        "nu": 0.6,
        "boundarLayerConductanceSource": "calc",
        "netRadiationSource": "calc",
        "defaultTimeOfMaximumTemperature": 14.0,
        "InitialValues": InitialValues if InitialValues is not None else None,
        "weather_Tav": weather_Tav,
        "weather_Amp": weather_Amp,
        "weather_Latitude": weather_Latitude,
        "clock_Today_DayOfYear": clock_Today_DayOfYear,
    }
    return state


def soil_temperature_process(
    state,
    waterBalance_SW,
    waterBalance_Es,
    waterBalance_Eos,
    waterBalance_Eo,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    weather_Wind,
    weather_AirPressure,
    weather_Radn,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Tav,
    weather_Latitude,
    clock_Today_DayOfYear
):
    airNode = state["airNode"]
    surfaceNode = state["surfaceNode"]
    topsoilNode = state["topsoilNode"]
    numPhantomNodes = state["numPhantomNodes"]
    doInitialisationStuff = state["doInitialisationStuff"]
    internalTimeStep = state["internalTimeStep"]
    timeOfDaySecs = state["timeOfDaySecs"]
    numLayers = state["numLayers"]
    numNodes = state["numNodes"]
    nodeDepth = state["nodeDepth"]
    volSpecHeatSoil = state["volSpecHeatSoil"]
    soilTemp = state["soilTemp"]
    morningSoilTemp = state["morningSoilTemp"]
    heatStorage = state["heatStorage"]
    thermalConductance = state["thermalConductance"]
    thermalConductivity = state["thermalConductivity"]
    boundaryLayerConductance = state["boundaryLayerConductance"]
    newTemperature = state["newTemperature"]
    airTemperature = state["airTemperature"]
    maxTempYesterday = state["maxTempYesterday"]
    minTempYesterday = state["minTempYesterday"]
    soilWater = state["soilWater"]
    minSoilTemp = state["minSoilTemp"]
    maxSoilTemp = state["maxSoilTemp"]
    aveSoilTemp = state["aveSoilTemp"]
    thickness = state["thickness"]
    bulkDensity = state["bulkDensity"]
    rocks = state["rocks"]
    carbon = state["carbon"]
    sand = state["sand"]
    silt = state["silt"]
    clay = state["clay"]
    soilRoughnessHeight = state["soilRoughnessHeight"]
    instrumentHeight = state["instrumentHeight"]
    netRadiation = state["netRadiation"]
    canopyHeight = state["canopyHeight"]
    nu = state["nu"]
    boundarLayerConductanceSource = state["boundarLayerConductanceSource"]
    netRadiationSource = state["netRadiationSource"]
    defaultTimeOfMaximumTemperature = state["defaultTimeOfMaximumTemperature"]
    InitialValues = state["InitialValues"]
    timestep = 24.0 * 60.0 * 60.0
    latentHeatOfVapourisation = 2465000.0
    if waterBalance_SW is not None:
        for i in range(numLayers):
            soilWater[i + 1] = waterBalance_SW[i]
        soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight if microClimate_CanopyHeight is not None else 0.0, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)
    if doInitialisationStuff:
        if ValuesInArray(InitialValues):
            soilTemp_local = [0.0] * (numNodes + 2)
            for i in range(len(InitialValues)):
                soilTemp_local[topsoilNode + i] = InitialValues[i]
            for i in range(len(soilTemp_local)):
                soilTemp[i] = soilTemp_local[i]
        else:
            soilTemp_init = calcSoilTemperature(numNodes, thickness, state["weather_Tav"], state["weather_Amp"], state["clock_Today_DayOfYear"], state["weather_Latitude"])
            for i in range(len(soilTemp)):
                if i < len(soilTemp_init):
                    soilTemp[i] = soilTemp_init[i]
            InitialValues = [0.0] * numLayers
            for i in range(numLayers):
                InitialValues[i] = soilTemp[topsoilNode + i]
            state["InitialValues"] = InitialValues
        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        for i in range(len(newTemperature)):
            newTemperature[i] = soilTemp[i]
        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn_list, cloudFr, cva = doNetRadiation(interactionsPerDay, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT)
    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0
    internalTimeStep = round(timestep / interactionsPerDay)
    volSpecHeatSoil = doVolumetricSpecificHeat(numNodes, soilWater, surfaceNode, nodeDepth, thickness)
    thermalConductivity = doThermalConductivity(numNodes, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, 1.3, 2.63, surfaceNode, nodeDepth, thickness)
    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0, defaultTimeOfMaximumTemperature, maxTempYesterday, minTempYesterday, weather_MaxT, weather_MinT, weather_MeanT)
        newTemperature[airNode] = airTemperature
        netRadiation = interpolateNetRadiation(solarRadn_list[timeStepIteration], cloudFr, cva, internalTimeStep, waterBalance_Eos, waterBalance_Eo, waterBalance_Salb, soilTemp[surfaceNode], airTemperature)
        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = 20.0
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, surfaceNode, canopyHeight, instrumentHeight, weather_Wind, airTemperature, weather_AirPressure, waterBalance_Eos, waterBalance_Eo)
            for _ in range(1, 1 + 1):
                newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, volSpecHeatSoil, thermalConductivity, soilTemp, surfaceNode, numNodes, nodeDepth, internalTimeStep, nu, netRadiation, netRadiationSource, waterBalance_Eos, waterBalance_Es, timestep, latentHeatOfVapourisation)
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, surfaceNode, canopyHeight, instrumentHeight, weather_Wind, airTemperature, weather_AirPressure, waterBalance_Eos, waterBalance_Eo)
        newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, volSpecHeatSoil, thermalConductivity, soilTemp, surfaceNode, numNodes, nodeDepth, internalTimeStep, nu, netRadiation, netRadiationSource, waterBalance_Eos, waterBalance_Es, timestep, latentHeatOfVapourisation)
        soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = doUpdate(interactionsPerDay, newTemperature, soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, timeOfDaySecs, internalTimeStep, boundaryLayerConductance, thermalConductivity[airNode], surfaceNode, numNodes)
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            for i in range(len(morningSoilTemp)):
                morningSoilTemp[i] = soilTemp[i]
    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT
    state["doInitialisationStuff"] = doInitialisationStuff
    state["internalTimeStep"] = internalTimeStep
    state["timeOfDaySecs"] = timeOfDaySecs
    state["volSpecHeatSoil"] = volSpecHeatSoil
    state["soilTemp"] = soilTemp
    state["morningSoilTemp"] = morningSoilTemp
    state["heatStorage"] = heatStorage
    state["thermalConductance"] = thermalConductance
    state["thermalConductivity"] = thermalConductivity
    state["boundaryLayerConductance"] = boundaryLayerConductance
    state["newTemperature"] = newTemperature
    state["airTemperature"] = airTemperature
    state["maxTempYesterday"] = maxTempYesterday
    state["minTempYesterday"] = minTempYesterday
    state["soilWater"] = soilWater
    state["minSoilTemp"] = minSoilTemp
    state["maxSoilTemp"] = maxSoilTemp
    state["aveSoilTemp"] = aveSoilTemp
    state["instrumentHeight"] = instrumentHeight
    state["netRadiation"] = netRadiation
    state["canopyHeight"] = canopyHeight
    FinalSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        FinalSoilTemperature[i] = soilTemp[topsoilNode + i]
    FinalSoilSurfaceTemperature = soilTemp[surfaceNode]
    AverageSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        AverageSoilTemperature[i] = aveSoilTemp[topsoilNode + i]
    AverageSoilSurfaceTemperature = aveSoilTemp[surfaceNode]
    MinimumSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        MinimumSoilTemperature[i] = minSoilTemp[topsoilNode + i]
    MinimumSoilSurfaceTemperature = minSoilTemp[surfaceNode]
    MaximumSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        MaximumSoilTemperature[i] = maxSoilTemp[topsoilNode + i]
    MaximumSoilSurfaceTemperature = maxSoilTemp[surfaceNode]
    BoundaryLayerConductance = boundaryLayerConductance
    ThermalConductivity_result = [0.0] * numNodes
    for i in range(numNodes):
        ThermalConductivity_result[i] = thermalConductivity[i + 1]
    HeatCapacity = [0.0] * numNodes
    for i in range(numNodes):
        HeatCapacity[i] = volSpecHeatSoil[surfaceNode + i]
    HeatStore = [0.0] * numNodes
    for i in range(numNodes):
        HeatStore[i] = heatStorage[surfaceNode + i]
    Thr_profile = [0.0] * (numNodes + 2)
    for i in range(len(Thr_profile)):
        Thr_profile[i] = morningSoilTemp[i]
    outputs = {
        "FinalSoilTemperature": FinalSoilTemperature,
        "FinalSoilSurfaceTemperature": FinalSoilSurfaceTemperature,
        "AverageSoilTemperature": AverageSoilTemperature,
        "AverageSoilSurfaceTemperature": AverageSoilSurfaceTemperature,
        "MinimumSoilTemperature": MinimumSoilTemperature,
        "MinimumSoilSurfaceTemperature": MinimumSoilSurfaceTemperature,
        "MaximumSoilTemperature": MaximumSoilTemperature,
        "MaximumSoilSurfaceTemperature": MaximumSoilSurfaceTemperature,
        "BoundaryLayerConductance": BoundaryLayerConductance,
        "ThermalConductivity": ThermalConductivity_result,
        "HeatCapacity": HeatCapacity,
        "HeatStore": HeatStore,
        "Thr_profile": Thr_profile
    }
    return state, outputs


def math_pi():
    return 3.141592653589793


def math_sin(x):
    import math
    return math.sin(x)


def math_cos(x):
    import math
    return math.cos(x)


def math_sqrt(x):
    import math
    if x < 0:
        return 0.0
    return math.sqrt(x)


def math_log(x):
    import math
    if x <= 0:
        return 0.0
    return math.log(x)


def math_exp(x):
    import math
    return math.exp(x)


def test_initialize_and_process_simple():
    physical_Thickness = [150.0, 150.0, 300.0, 300.0, 300.0]
    physical_BD = [1.3, 1.3, 1.35, 1.35, 1.4]
    physical_Rocks = [0.0] * 5
    physical_ParticleSizeSand = [40.0] * 5
    physical_ParticleSizeSilt = [40.0] * 5
    physical_ParticleSizeClay = [20.0] * 5
    waterBalance_SW = [0.25, 0.25, 0.26, 0.27, 0.28]
    organic_Carbon = [1.0, 0.8, 0.6, 0.4, 0.3]
    weather_Tav = 15.0
    weather_Amp = 10.0
    weather_Latitude = -27.0
    clock_Today_DayOfYear = 200
    state = soil_temperature_init(
        physical_Thickness,
        physical_BD,
        physical_Rocks,
        physical_ParticleSizeSand,
        physical_ParticleSizeSilt,
        physical_ParticleSizeClay,
        waterBalance_SW,
        organic_Carbon,
        weather_Tav,
        weather_Amp,
        weather_Latitude,
        clock_Today_DayOfYear,
        InitialValues=None,
        DepthToConstantTemperature=10000.0,
        instrumHeight=1.2
    )
    waterBalance_Es = 2.0
    waterBalance_Eos = 5.0
    waterBalance_Eo = 6.0
    waterBalance_Salb = 0.2
    microClimate_CanopyHeight = 0.0
    weather_Wind = 3.0
    weather_AirPressure = 1010.0
    weather_Radn = 18.0
    weather_MeanT = 18.0
    weather_MaxT = 25.0
    weather_MinT = 10.0
    weather_Tav_day = weather_Tav
    weather_Latitude_day = weather_Latitude
    clock_Today_DayOfYear_day = clock_Today_DayOfYear
    state2, outputs = soil_temperature_process(
        state,
        waterBalance_SW,
        waterBalance_Es,
        waterBalance_Eos,
        waterBalance_Eo,
        waterBalance_Salb,
        microClimate_CanopyHeight,
        weather_Wind,
        weather_AirPressure,
        weather_Radn,
        weather_MeanT,
        weather_MaxT,
        weather_MinT,
        weather_Tav_day,
        weather_Latitude_day,
        clock_Today_DayOfYear_day
    )
    assert isinstance(outputs["FinalSoilTemperature"], list)
    assert len(outputs["FinalSoilTemperature"]) == len(physical_Thickness)