def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and not (Value is None):
                # In C# they also checked for NaN. Here we assume None is not a number.
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


def volumetricFractionSand(sand, carbon, rocks, bulkDensity, ps, pom, layer):
    return (1.0 - volumetricFractionOrganicMatter(carbon, bulkDensity, pom, layer) - volumetricFractionRocks(rocks, layer)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(silt, carbon, rocks, bulkDensity, ps, pom, layer):
    return (1.0 - volumetricFractionOrganicMatter(carbon, bulkDensity, pom, layer) - volumetricFractionRocks(rocks, layer)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(clay, carbon, rocks, bulkDensity, ps, pom, layer):
    return (1.0 - volumetricFractionOrganicMatter(carbon, bulkDensity, pom, layer) - volumetricFractionRocks(rocks, layer)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(soilWater, carbon, layer):
    return (1 - volumetricFractionOrganicMatter(carbon, None, 1.0, layer)) * soilWater[layer]


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, pom, layer):
    return 1.0 - volumetricFractionRocks(rocks, layer) - volumetricFractionOrganicMatter(carbon, bulkDensity, pom, layer) - volumetricFractionSand(sand, carbon, rocks, bulkDensity, ps, pom, layer) - volumetricFractionSilt(silt, carbon, rocks, bulkDensity, ps, pom, layer) - volumetricFractionClay(clay, carbon, rocks, bulkDensity, ps, pom, layer) - soilWater[layer] - volumetricFractionIce(layer)


def ThermalConductance(name, layer, rocks, sand, silt, clay, bulkDensity, ps, pom):
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
        # Replicate odd formula then overwrite as in original code
        _ = (thermalConductanceRocks ** volumetricFractionRocks(rocks, layer)) * (thermalConductanceSand ** volumetricFractionSand(sand, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer))
        _ = _ + (thermalConductanceSilt ** volumetricFractionSilt(silt, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer))
        _ = _ + (thermalConductanceClay ** volumetricFractionClay(clay, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer))
        result = _

    # As in original code, overwrite with volumetricSpecificHeat (wrong but preserved)
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(name, layer, rocks, sand, silt, clay, soilWater, bulkDensity, ps, pom):
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
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / (soilWater[layer] + volumetricFractionIce(layer) + volumetricFractionAir(rocks, carbon=None, sand=sand, silt=silt, clay=clay, soilWater=soilWater, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(rocks, carbon=None, sand=sand, silt=silt, clay=clay, soilWater=soilWater, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer) / (soilWater[layer] + volumetricFractionIce(layer) + volumetricFractionAir(rocks, carbon=None, sand=sand, silt=silt, clay=clay, soilWater=soilWater, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(rocks, layer) + shapeFactorSand * volumetricFractionSand(sand, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer) + shapeFactorSilt * volumetricFractionSilt(silt, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer) + shapeFactorClay * volumetricFractionClay(clay, carbon=None, rocks=rocks, bulkDensity=bulkDensity, ps=ps, pom=pom, layer=layer)

    # As in original code, overwrite with volumetricSpecificHeat (wrong but preserved)
    result = volumetricSpecificHeat(name, layer)
    return result


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def mapLayer2Node(layerArray, nodeArray, surfaceNode, numNodes, nodeDepth, thickness):
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)


def doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay):
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 1 + 0):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5) if clay[layer] > 0 else 0, 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, surfaceNode):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * 3.141592653589793 / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = (2.0 * dh / w) ** 0.5
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    temp_profile = [0.0] * (numNodes + 1 + 1)
    for node in range(1, numNodes + 1):
        temp_profile[node] = weather_Tav + weather_Amp * (2.718281828459045 ** (-1.0 * cumulativeDepth[node] / zd)) * (math_sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * 3.141592653589793 - cumulativeDepth[node] / zd))
    # As in original: copy from index 0 into soilTemp starting at surfaceNode for numNodes elements (this places 0 at index 1)
    soilTempIO = [0.0] * (numNodes + 1 + 1)
    for i in range(0, numNodes):
        dst_index = surfaceNode + i
        if dst_index < len(soilTempIO) and i < len(temp_profile):
            soilTempIO[dst_index] = temp_profile[i]
    return soilTempIO


def interpolateTemperature(timeHours, weather_MaxT, weather_MinT, weather_MeanT, maxTempYesterday, minTempYesterday, defaultTimeOfMaximumTemperature):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math_sin((0.0 + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (maxTempYesterday - minTempYesterday) / 2.0 + (maxTempYesterday + minTempYesterday) / 2.0
        tScale = Divide((minT_time - time), minT_time, 0.0)
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = math_sin((time + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


def longWaveRadn(emissivity, tDegC):
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def doNetRadiation(ITERATIONSperDAY, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT):
    TSTEPS2RAD = Divide(2.0 * 3.141592653589793, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin(4.869 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25) + 0.03345 * math_sin(6.224 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25)))
    cD = (1.0 - solarDeclination * solarDeclination) ** 0.5
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(weather_Latitude * 3.141592653589793 / 180.0) + cD * math_cos(weather_Latitude * 3.141592653589793 / 180.0) * math_cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
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
    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = math_exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(solarRadn, cloudFr, cva, internalTimeStep, airTemperature, surfaceTemperature, waterBalance_Salb, waterBalance_Eo, waterBalance_Eos):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, surfaceTemperature) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * (max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5) + waterBalance_Salb * weather_MeanT
    return surfaceT


def doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, surfaceNode, nodeDepth, thickness):
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in [n for n in soilConstituentNames if n != "Minerals"]:
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, surfaceNode, numNodes, nodeDepth, thickness)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, soilConstituentNames, rocks, sand, silt, clay, soilWater, bulkDensity, ps, pom, surfaceNode, nodeDepth, thickness):
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, sand, silt, clay, soilWater, bulkDensity, ps, pom)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, sand, silt, clay, bulkDensity, ps, pom)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, sand, silt, clay, bulkDensity, ps, pom)
            k = (2.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** -1 + (1.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** -1
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = [0.0] * (numNodes + 1)
    mapLayer2Node(thermCondLayers, thermalConductivity, surfaceNode, numNodes, nodeDepth, thickness)
    return thermalConductivity


def doThomas(newTemps, soilTemp, numNodes, surfaceNode, airNode, volSpecHeatSoil, thermalConductivity, nodeDepth, internalTimeStep, nu, netRadiationSource, netRadiation_MJ, waterBalance_Eos, waterBalance_Es, timestep, thermalConductance):
    a = [0.0] * (numNodes + 1 + 1)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    heatStorage = [0.0] * (numNodes + 1)

    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + g * thermalConductance[node] * soilTemp[node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation_MJ * 1000000.0, internalTimeStep, 0.0)
    else:
        radnNet = Divide(waterBalance_Eos * 2465000.0, timestep, 0.0)

    latentHeatFlux = Divide(waterBalance_Es * 2465000.0, timestep, 0.0)

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


def getBoundaryLayerConductance(TNew_zb_surface, airTemperature, weather_AirPressure, weather_Wind, canopyHeight_m, instrumentHeight_m, waterBalance_Eo, waterBalance_Eos):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)

    roughnessFactorMomentum = 0.13 * canopyHeight_m
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight_m
    surfaceTemperature = TNew_zb_surface

    diffusePenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    stefanBoltzmannConstant = 0.0000000567
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for _ in range(1, 3 + 1):
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, math_log(Divide(instrumentHeight_m - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum, 0.0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, math_log(Divide(instrumentHeight_m - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat, 0.0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight_m * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + (1.0 - 16.0 * stabilityParammeter) ** 0.5) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


def math_sin(x):
    import math
    return math.sin(x)


def math_cos(x):
    import math
    return math.cos(x)


def math_exp(x):
    import math
    return math.exp(x)


def math_log(x):
    import math
    return math.log(x)


def soiltemperature_initialize(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    waterBalance_Eos,
    waterBalance_Eo,
    waterBalance_Salb,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Radn,
    weather_AirPressure,
    weather_Wind,
    clock_Today_DayOfYear,
    microClimate_CanopyHeight,
    instrumHeight=None,
    DepthToConstantTemperature=10000.0,
    InitialValues=None,
    defaultTimeOfMaximumTemperature=14.0,
    timestep_seconds=24.0 * 60.0 * 60.0
):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5

    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0

    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    # thickness with phantom layers
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(numLayers):
        thickness[i + 1] = physical_Thickness[i]
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes

    nodeDepth = [0.0] * (numNodes + 1 + 1)
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
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * (thickness[layer - 1] if (layer - 1) >= 0 else 0.0), thickness[layer], 0.0)
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
    soilTemp = [0.0] * (numNodes + 1 + 1)
    morningSoilTemp = [0.0] * (numNodes + 1 + 1)
    newTemperature = [0.0] * (numNodes + 1 + 1)
    thermalConductivity = [0.0] * (numNodes + 1)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 1 + 1)

    # instrument height setup
    defaultInstrumentHeight = 1.2
    instrumentHeight = instrumHeight if (instrumHeight is not None and instrumHeight > 0.00001) else defaultInstrumentHeight

    # Thermal conductivity coefficients (not used elsewhere but kept to preserve behaviour)
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay)

    # initial soil temperature profile (seasonal)
    initial_profile = calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, surfaceNode)
    # Copy result to soilTemp as in original readParam behaviour
    for i in range(len(initial_profile)):
        if i < len(soilTemp):
            soilTemp[i] = initial_profile[i]
    for i in range(len(soilTemp)):
        if i < len(newTemperature):
            newTemperature[i] = soilTemp[i]

    soilRoughnessHeight = 57.0

    state = {
        "doInitialisationStuff": doInitialisationStuff,
        "internalTimeStep": internalTimeStep,
        "timeOfDaySecs": timeOfDaySecs,
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
        "maxTempYesterday": weather_MaxT,
        "minTempYesterday": weather_MinT,
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
        "InitialValues": InitialValues,
        "defaultTimeOfMaximumTemperature": defaultTimeOfMaximumTemperature,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "timestep": timestep_seconds,
        "pom": 1.3,
        "ps": 2.63,
        "soilConstituentNames": ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"],
    }
    return state


def soiltemperature_process(
    state,
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    waterBalance_Es,
    waterBalance_Eos,
    waterBalance_Eo,
    waterBalance_Salb,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Radn,
    weather_AirPressure,
    weather_Wind,
    clock_Today_DayOfYear,
    microClimate_CanopyHeight
):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5

    # Unpack state
    doInitialisationStuff = state["doInitialisationStuff"]
    internalTimeStep = state["internalTimeStep"]
    timeOfDaySecs = state["timeOfDaySecs"]
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    nodeDepth = state["nodeDepth"]
    thermCondPar1 = state["thermCondPar1"]
    thermCondPar2 = state["thermCondPar2"]
    thermCondPar3 = state["thermCondPar3"]
    thermCondPar4 = state["thermCondPar4"]
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
    instrumHeight = state["instrumHeight"]
    nu = state["nu"]
    boundarLayerConductanceSource = state["boundarLayerConductanceSource"]
    netRadiationSource = state["netRadiationSource"]
    InitialValues = state["InitialValues"]
    defaultTimeOfMaximumTemperature = state["defaultTimeOfMaximumTemperature"]
    DepthToConstantTemperature = state["DepthToConstantTemperature"]
    timestep = state["timestep"]
    pom = state["pom"]
    ps = state["ps"]
    soilConstituentNames = state["soilConstituentNames"]

    # Update other variables
    if waterBalance_SW is not None:
        for i in range(numLayers):
            if (1 + i) < len(soilWater):
                soilWater[1 + i] = waterBalance_SW[i]
    soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)

    # Initialisation block on first process call
    if doInitialisationStuff:
        if ValuesInArray(InitialValues):
            soilTemp = [0.0] * (numNodes + 1 + 1)
            for i in range(len(InitialValues)):
                if topsoilNode + i < len(soilTemp):
                    soilTemp[topsoilNode + i] = InitialValues[i]
        else:
            profile = calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude, surfaceNode)
            soilTemp = [x for x in soilTemp]
            for i in range(len(profile)):
                if i < len(soilTemp):
                    soilTemp[i] = profile[i]
            InitialValues = [0.0] * numLayers
            for i in range(numLayers):
                if topsoilNode + i < len(soilTemp):
                    InitialValues[i] = soilTemp[topsoilNode + i]
        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        newTemperature = [x for x in soilTemp]
        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False

    # Daily process
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, clock_Today_DayOfYear, weather_Latitude, weather_Radn, weather_MinT)

    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0

    internalTimeStep = round(timestep / interactionsPerDay)

    volSpecHeatSoil = doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, surfaceNode, nodeDepth, thickness)
    thermalConductivity = doThermalConductivity(numNodes, soilConstituentNames, rocks, sand, silt, clay, soilWater, bulkDensity, ps, pom, surfaceNode, nodeDepth, thickness)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0, weather_MaxT, weather_MinT, weather_MeanT, maxTempYesterday, minTempYesterday, defaultTimeOfMaximumTemperature)
        newTemperature[airNode] = airTemperature

        netRadiation = interpolateNetRadiation(
            solarRadn[timeStepIteration],
            cloudFr,
            cva,
            internalTimeStep,
            airTemperature,
            soilTemp[surfaceNode],
            waterBalance_Salb,
            waterBalance_Eo,
            waterBalance_Eos
        )

        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = 20.0
        else:
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature[surfaceNode], airTemperature, weather_AirPressure, weather_Wind, canopyHeight, instrumentHeight, waterBalance_Eo, waterBalance_Eos)
            for _ in range(1, 1 + 1):
                newTemperature, heatStorage, thermalConductance = doThomas(
                    newTemperature,
                    soilTemp,
                    numNodes,
                    surfaceNode,
                    airNode,
                    volSpecHeatSoil,
                    thermalConductivity,
                    nodeDepth,
                    internalTimeStep,
                    nu,
                    netRadiationSource,
                    netRadiation,
                    waterBalance_Eos,
                    waterBalance_Es,
                    timestep,
                    thermalConductance
                )
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature[surfaceNode], airTemperature, weather_AirPressure, weather_Wind, canopyHeight, instrumentHeight, waterBalance_Eo, waterBalance_Eos)

        newTemperature, heatStorage, thermalConductance = doThomas(
            newTemperature,
            soilTemp,
            numNodes,
            surfaceNode,
            airNode,
            volSpecHeatSoil,
            thermalConductivity,
            nodeDepth,
            internalTimeStep,
            nu,
            netRadiationSource,
            netRadiation,
            waterBalance_Eos,
            waterBalance_Es,
            timestep,
            thermalConductance
        )

        # doUpdate
        soilTemp = [x for x in newTemperature]
        if timeOfDaySecs < internalTimeStep * 1.2:
            for node in range(surfaceNode, numNodes + 1):
                minSoilTemp[node] = soilTemp[node]
                maxSoilTemp[node] = soilTemp[node]
        for node in range(surfaceNode, numNodes + 1):
            if soilTemp[node] < minSoilTemp[node]:
                minSoilTemp[node] = soilTemp[node]
            elif soilTemp[node] > maxSoilTemp[node]:
                maxSoilTemp[node] = soilTemp[node]
            aveSoilTemp[node] += Divide(soilTemp[node], interactionsPerDay, 0.0)
        boundaryLayerConductance += Divide(thermalConductivity[airNode], interactionsPerDay, 0.0)

        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            morningSoilTemp = [x for x in soilTemp]

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    # Prepare outputs
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

    ThermalConductivity_profile = [0.0] * (numNodes)
    for i in range(numNodes):
        ThermalConductivity_profile[i] = thermalConductivity[1 + i]

    HeatCapacity_profile = [0.0] * (numNodes)
    for i in range(numNodes):
        HeatCapacity_profile[i] = volSpecHeatSoil[surfaceNode + i]

    HeatStore_profile = [0.0] * (numNodes)
    for i in range(numNodes):
        HeatStore_profile[i] = heatStorage[surfaceNode + i]

    outputs = {
        "FinalSoilTemperature": FinalSoilTemperature,
        "FinalSoilSurfaceTemperature": FinalSoilSurfaceTemperature,
        "AverageSoilTemperature": AverageSoilTemperature,
        "AverageSoilSurfaceTemperature": AverageSoilSurfaceTemperature,
        "MinimumSoilTemperature": MinimumSoilTemperature,
        "MinimumSoilSurfaceTemperature": MinimumSoilSurfaceTemperature,
        "MaximumSoilTemperature": MaximumSoilTemperature,
        "MaximumSoilSurfaceTemperature": MaximumSoilSurfaceTemperature,
        "BoundaryLayerConductance": boundaryLayerConductance,
        "ThermalConductivity": ThermalConductivity_profile,
        "HeatCapacity": HeatCapacity_profile,
        "HeatStore": HeatStore_profile,
        "Thr_profile": [x for x in morningSoilTemp],
    }

    # Pack state back
    state["doInitialisationStuff"] = doInitialisationStuff
    state["internalTimeStep"] = internalTimeStep
    state["timeOfDaySecs"] = timeOfDaySecs
    state["nodeDepth"] = nodeDepth
    state["thermCondPar1"] = thermCondPar1
    state["thermCondPar2"] = thermCondPar2
    state["thermCondPar3"] = thermCondPar3
    state["thermCondPar4"] = thermCondPar4
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
    state["thickness"] = thickness
    state["bulkDensity"] = bulkDensity
    state["rocks"] = rocks
    state["carbon"] = carbon
    state["sand"] = sand
    state["silt"] = silt
    state["clay"] = clay
    state["soilRoughnessHeight"] = soilRoughnessHeight
    state["instrumentHeight"] = instrumentHeight
    state["netRadiation"] = netRadiation
    state["canopyHeight"] = canopyHeight
    state["instrumHeight"] = instrumHeight
    state["nu"] = nu
    state["boundarLayerConductanceSource"] = boundarLayerConductanceSource
    state["netRadiationSource"] = netRadiationSource
    state["InitialValues"] = InitialValues
    state["defaultTimeOfMaximumTemperature"] = defaultTimeOfMaximumTemperature
    state["DepthToConstantTemperature"] = DepthToConstantTemperature
    state["timestep"] = timestep
    state["pom"] = pom
    state["ps"] = ps
    state["soilConstituentNames"] = soilConstituentNames

    return state, outputs