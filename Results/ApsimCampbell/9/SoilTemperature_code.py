def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and not (Value != Value):  # NaN check
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
    precisionMargin = 0.00001
    lowerBound = Lower - precisionMargin
    upperBound = Upper + precisionMargin
    return


def kelvinT(celciusT):
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def longWaveRadn(emissivity, tDegC):
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


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


def volumetricFractionRocks(rocks, layer):
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(carbon, bulkDensity, pom, layer):
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(sand, bulkDensity, ps, carbon, rocks, layer):
    return (1 - volumetricFractionOrganicMatter(carbon, bulkDensity, 1.3, layer) - volumetricFractionRocks(rocks, layer)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(silt, bulkDensity, ps, carbon, rocks, layer):
    return (1 - volumetricFractionOrganicMatter(carbon, bulkDensity, 1.3, layer) - volumetricFractionRocks(rocks, layer)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(clay, bulkDensity, ps, carbon, rocks, layer):
    return (1 - volumetricFractionOrganicMatter(carbon, bulkDensity, 1.3, layer) - volumetricFractionRocks(rocks, layer)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(soilWater, carbon, bulkDensity, layer):
    return (1 - volumetricFractionOrganicMatter(carbon, bulkDensity, 1.3, layer)) * soilWater[layer]


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, layer):
    return 1.0 - volumetricFractionRocks(rocks, layer) - \
           volumetricFractionOrganicMatter(carbon, bulkDensity, 1.3, layer) - \
           volumetricFractionSand(sand, bulkDensity, ps, carbon, rocks, layer) - \
           volumetricFractionSilt(silt, bulkDensity, ps, carbon, rocks, layer) - \
           volumetricFractionClay(clay, bulkDensity, ps, carbon, rocks, layer) - \
           volumetricFractionWater(soilWater, carbon, bulkDensity, layer) - \
           volumetricFractionIce(layer)


def ThermalConductance(name, layer, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps):
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
        result = (thermalConductanceRocks ** volumetricFractionRocks(rocks, layer)) * \
                 (thermalConductanceSand ** volumetricFractionSand(sand, bulkDensity, ps, carbon, rocks, layer)) + \
                 (thermalConductanceSilt ** volumetricFractionSilt(silt, bulkDensity, ps, carbon, rocks, layer)) + \
                 (thermalConductanceClay ** volumetricFractionClay(clay, bulkDensity, ps, carbon, rocks, layer))
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(name, layer, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps):
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
        den = volumetricFractionWater(soilWater, carbon, bulkDensity, layer) + volumetricFractionIce(layer) + volumetricFractionAir(rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, layer)
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / den if den != 0 else 0.333
        return result
    elif name == "Air":
        den = volumetricFractionWater(soilWater, carbon, bulkDensity, layer) + volumetricFractionIce(layer) + volumetricFractionAir(rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, layer)
        result = 0.333 - 0.333 * volumetricFractionAir(rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, layer) / den if den != 0 else 0.333
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(rocks, layer) + \
                 shapeFactorSand * volumetricFractionSand(sand, bulkDensity, ps, carbon, rocks, layer) + \
                 shapeFactorSilt * volumetricFractionSilt(silt, bulkDensity, ps, carbon, rocks, layer) + \
                 shapeFactorClay * volumetricFractionClay(clay, bulkDensity, ps, carbon, rocks, layer)
    result = volumetricSpecificHeat(name, layer)
    return result


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def mapLayer2Node(layerArray, nodeArray, nodeDepth, thickness, numNodes):
    surfaceNode = 1
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0) + Divide(layerArray[layer + 1] * d2, dSum, 0)


def doThermalConductivityCoeffs(numNodes, numLayers, bulkDensity):
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 1 + 1):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clamp_sqrt(bulkDensity, layer)), 0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def clamp_sqrt(bulkDensity, layer):
    import math
    val = bulkDensity[layer]
    if val < 0:
        return 0.0
    return math.sqrt(val) if val > 0 else 0.0


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * ((max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5)) + waterBalance_Salb * weather_MeanT
    boundCheck(surfaceT, -100.0, 100.0, "Initial surfaceT")
    return surfaceT


def calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, soilTempIO):
    import math
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTemp = [0.0] * (numNodes + 2)
    for nodes in range(1, numNodes + 1):
        soilTemp[nodes] = weather_Tav + weather_Amp * math.exp(-1 * cumulativeDepth[nodes] / zd) * math.sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)
    for i in range(0, numNodes):
        if 1 + i < len(soilTempIO) and i < len(soilTemp):
            soilTempIO[1 + i] = soilTemp[i]


def interpolateTemperature(timeHours, defaultTimeOfMaximumTemperature, maxTempYesterday, minTempYesterday, weather_MinT, weather_MaxT, weather_MeanT):
    import math
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


def doVolumetricSpecificHeat(numNodes, soilWater, nodeDepth, thickness):
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName == "Minerals":
                continue
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, nodeDepth, thickness, numNodes)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps, nodeDepth, thickness):
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, ps)
            k = (2.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** -1) + \
                (1.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** -1)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = numerator / denominator if denominator != 0 else 0.0
    thermalConductivity = [0.0] * (numNodes + 1)
    mapLayer2Node(thermCondLayers, thermalConductivity, nodeDepth, thickness, numNodes)
    return thermalConductivity


def doThomas(newTemps, numNodes, nodeDepth, volSpecHeatSoil, internalTimeStep, thermalConductivity, soilTemp, netRadiation, netRadiationSource, waterBalance_Eos, waterBalance_Es, timestep):
    surfaceNode = 1
    airNode = 0
    latentHeatOfVapourisation = 2465000.0
    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 2)
    heatStorage = [0.0] * (numNodes + 1)
    nu = 0.6
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
    radnNet = 0.0
    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0)
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0)
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
        boundCheck(newTemps[node], -50.0, 100.0, "newTemps(" + str(node) + ")")
    return newTemps, thermalConductance, heatStorage


def doUpdate(numInterationsPerDay, numNodes, soilTemp, newTemperature, minSoilTemp, maxSoilTemp, aveSoilTemp, thermalConductivity, timeOfDaySecs, internalTimeStep):
    surfaceNode = 1
    airNode = 0
    for i in range(len(newTemperature)):
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
        aveSoilTemp[node] += Divide(soilTemp[node], numInterationsPerDay, 0)
    boundaryLayerConductanceIncrement = Divide(thermalConductivity[airNode], numInterationsPerDay, 0)
    return boundaryLayerConductanceIncrement


def getBoundaryLayerConductance(TNew_zb, canopyHeight, instrumentHeight, airTemperature, weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo):
    import math
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight
    surfaceTemperature = TNew_zb[1]
    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    stefanBoltzmannConstant = 0.0000000567
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)
    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0
    for iteration in range(1, 3 + 1):
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant,
                                  math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0)) + stabilityCorrectionMomentum,
                                  0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity,
                                   math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0)) + stabilityCorrectionHeat,
                                   0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity,
                                     SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + math.sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0) if (1.0 - 16.0 * stabilityParammeter) > 0 else 0.0
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


def doNetRadiation(ITERATIONSperDAY, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT):
    import math
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin(4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + 0.03345 * math.sin(6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))
    cD = math.sqrt(max(0.0, 1.0 - solarDeclination * solarDeclination))
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


def interpolateNetRadiation(solarRadn, cloudFr, cva, internalTimeStep, soilTemp_surface, airTemperature, waterBalance_Eos, waterBalance_Eo):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilTemp_surface) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = 0.0  # actual shortwave outgoing handled via albedo in soilwat; here considered zero for net to soil through cover for this function
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def SoilTemperature_initialize(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    waterBalance_Salb,
    weather_Tav,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Amp,
    weather_Radn,
    weather_AirPressure,
    weather_Latitude,
    weather_Wind,
    clock_Today_DayOfYear,
    microClimate_CanopyHeight,
    InitialValues=None,
    DepthToConstantTemperature=10000.0,
    instrumHeight=None
):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    defaultInstrumentHeight = 1.2
    bareSoilRoughness = 57.0
    boundCheck(weather_Tav, -30.0, 50.0, "tav (oC)")
    instrumentHeight = instrumHeight if (instrumHeight is not None and instrumHeight > 0.00001) else defaultInstrumentHeight
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(numLayers):
        thickness[1 + i] = physical_Thickness[i]
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
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0)
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
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numNodes, numLayers, bulkDensity)
    calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, soilTemp)
    for i in range(len(soilTemp)):
        newTemperature[i] = soilTemp[i]
    soilRoughnessHeight = bareSoilRoughness
    if ValuesInArray(InitialValues):
        soilTemp = [0.0] * (numNodes + 2)
        for i in range(len(InitialValues)):
            soilTemp[topsoilNode + i] = InitialValues[i]
    else:
        calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, soilTemp)
        InitialValuesComputed = [0.0] * numLayers
        for i in range(numLayers):
            InitialValuesComputed[i] = soilTemp[topsoilNode + i]
    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav
    for i in range(len(soilTemp)):
        newTemperature[i] = soilTemp[i]
    maxTempYesterday = weather_MaxT
    minTempYesterday = weather_MinT
    state = {
        "numNodes": numNodes,
        "numLayers": numLayers,
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
        "maxTempYesterday": maxTempYesterday,
        "minTempYesterday": minTempYesterday,
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
        "instrumHeight": instrumHeight if instrumHeight is not None else 0.0,
        "nu": 0.6,
        "boundarLayerConductanceSource": "calc",
        "netRadiationSource": "calc",
        "defaultTimeOfMaximumTemperature": 14.0,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "pom": 1.3,
        "ps": 2.63
    }
    return state


def SoilTemperature_process(
    state,
    waterBalance_SW,
    waterBalance_Eos,
    waterBalance_Es,
    waterBalance_Eo,
    waterBalance_Salb,
    weather_Tav,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Amp,
    weather_Radn,
    weather_AirPressure,
    weather_Latitude,
    weather_Wind,
    clock_Today_DayOfYear,
    microClimate_CanopyHeight
):
    airNode = 0
    surfaceNode = 1
    timestep = 24.0 * 60.0 * 60.0
    interactionsPerDay = 48
    soilWater = state["soilWater"]
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    for i in range(numLayers):
        soilWater[1 + i] = waterBalance_SW[i]
    soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, state["soilRoughnessHeight"]) / 1000.0
    instrumentHeight = max(state["instrumentHeight"], canopyHeight + 0.5)
    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT)
    Zero(state["minSoilTemp"])
    Zero(state["maxSoilTemp"])
    Zero(state["aveSoilTemp"])
    state["boundaryLayerConductance"] = 0.0
    internalTimeStep = round(timestep / interactionsPerDay)
    state["internalTimeStep"] = internalTimeStep
    state["volSpecHeatSoil"] = doVolumetricSpecificHeat(numNodes, soilWater, state["nodeDepth"], state["thickness"])
    state["thermalConductivity"] = doThermalConductivity(numNodes, state["rocks"], state["sand"], state["silt"], state["clay"], soilWater, state["carbon"], state["bulkDensity"], state["ps"], state["nodeDepth"], state["thickness"])
    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        state["timeOfDaySecs"] = timeOfDaySecs
        if timestep < 24.0 * 60.0 * 60.0:
            state["airTemperature"] = weather_MeanT
        else:
            state["airTemperature"] = interpolateTemperature(timeOfDaySecs / 3600.0, state["defaultTimeOfMaximumTemperature"], state["maxTempYesterday"], state["minTempYesterday"], weather_MinT, weather_MaxT, weather_MeanT)
        state["newTemperature"][airNode] = state["airTemperature"]
        state["netRadiation"] = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, internalTimeStep, state["soilTemp"][surfaceNode], state["airTemperature"], waterBalance_Eos, waterBalance_Eo)
        if state["boundarLayerConductanceSource"] == "constant":
            state["thermalConductivity"][airNode] = 20.0
        else:
            state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"], canopyHeight, instrumentHeight, state["airTemperature"], weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo)
            for iteration in range(1, 1 + 1):
                newTemps, thermCond, heatStor = doThomas(
                    state["newTemperature"][:],
                    numNodes,
                    state["nodeDepth"],
                    state["volSpecHeatSoil"],
                    internalTimeStep,
                    state["thermalConductivity"],
                    state["soilTemp"],
                    state["netRadiation"],
                    state["netRadiationSource"],
                    waterBalance_Eos,
                    waterBalance_Es,
                    timestep
                )
                state["newTemperature"] = newTemps
                state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"], canopyHeight, instrumentHeight, state["airTemperature"], weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo)
        newTemps, thermCond, heatStor = doThomas(
            state["newTemperature"][:],
            numNodes,
            state["nodeDepth"],
            state["volSpecHeatSoil"],
            internalTimeStep,
            state["thermalConductivity"],
            state["soilTemp"],
            state["netRadiation"],
            state["netRadiationSource"],
            waterBalance_Eos,
            waterBalance_Es,
            timestep
        )
        state["newTemperature"] = newTemps
        state["thermalConductance"] = thermCond
        state["heatStorage"] = heatStor
        inc = doUpdate(interactionsPerDay, numNodes, state["soilTemp"], state["newTemperature"], state["minSoilTemp"], state["maxSoilTemp"], state["aveSoilTemp"], state["thermalConductivity"], timeOfDaySecs, internalTimeStep)
        state["boundaryLayerConductance"] += inc
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            for i in range(len(state["soilTemp"])):
                state["morningSoilTemp"][i] = state["soilTemp"][i]
    state["minTempYesterday"] = weather_MinT
    state["maxTempYesterday"] = weather_MaxT
    return state