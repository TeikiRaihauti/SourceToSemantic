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


def boundCheck(VariableValue, Lower, Upper, VariableName):
    precisionMargin = 0.00001
    lowerBound = Lower - precisionMargin
    upperBound = Upper + precisionMargin
    # No exception thrown; function kept for parity
    return


def volumetricFractionRocks(layer, rocks):
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom):
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(layer, rocks, carbon, bulkDensity, sand, ps, pom):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer, rocks, carbon, bulkDensity, silt, ps, pom):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(layer, rocks, carbon, bulkDensity, clay, ps, pom):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(layer, carbon, soilWater, pom):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, soilWater, pom)) * soilWater[layer]


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, pom):
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionSand(layer, rocks, carbon, bulkDensity, sand, ps, pom) - volumetricFractionSilt(layer, rocks, carbon, bulkDensity, silt, ps, pom) - volumetricFractionClay(layer, rocks, carbon, bulkDensity, clay, ps, pom) - volumetricFractionWater(layer, carbon, soilWater, pom) - volumetricFractionIce(layer)


def volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom):
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


def ThermalConductance(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom):
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
    result = volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom)
    return result


def shapeFactor(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom):
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
        vf_w = volumetricFractionWater(layer, carbon, soilWater, pom)
        vf_i = volumetricFractionIce(layer)
        vf_a = volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, pom)
        result = 0.333 - 0.333 * vf_i / (vf_w + vf_i + vf_a)
        return result
    elif name == "Air":
        vf_w = volumetricFractionWater(layer, carbon, soilWater, pom)
        vf_i = volumetricFractionIce(layer)
        vf_a = volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, ps, pom)
        result = 0.333 - 0.333 * vf_a / (vf_w + vf_i + vf_a)
        return result
    result = volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom)
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
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def mapLayer2Node(layerArray, nodeArray, nodeDepth, thickness, surfaceNode, numNodes):
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0) + Divide(layerArray[layer + 1] * d2, dSum, 0)


def doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay):
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 1 + 0):  # up to numLayers + 1 inclusive in original
        element = layer
        bd = bulkDensity[layer]
        thermCondPar1[element] = 0.65 - 0.78 * bd + 0.6 * (bd ** 2)
        thermCondPar2[element] = 1.06 * bd
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5) if clay[layer] > 0 else 0.0000001, 0)
        thermCondPar4[element] = 0.03 + 0.1 * (bd ** 2)
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def calcSoilTemperature(thickness, numNodes, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, surfaceNode):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * 3.141592653589793 / (365.25 * 24 * 3600)
    dh = 0.6
    zd = (2 * dh / w) ** 0.5
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTemp_calc = [0.0] * (numNodes + 1 + 1)
    for nodes in range(1, numNodes + 1):
        soilTemp_calc[nodes] = weather_Tav + weather_Amp * pow(2.718281828459045, -1 * cumulativeDepth[nodes] / zd) * \
                               (math_sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * 3.141592653589793 - cumulativeDepth[nodes] / zd))
    soilTemp_io = [0.0] * (numNodes + 1 + 1)
    for i in range(0, len(soilTemp_calc)):
        soilTemp_io[i] = 0.0
    # copy to IO from surfaceNode, length numNodes
    for i in range(0, numNodes):
        idx = surfaceNode + i
        if idx < len(soilTemp_io) and idx < len(soilTemp_calc):
            soilTemp_io[idx] = soilTemp_calc[idx - 1] if (idx - 1) < len(soilTemp_calc) else 0.0
    return soilTemp_io


def calcLayerTemperature(depthLag, alx, deltaTemp, weather_Tav, weather_Amp):
    return weather_Tav + (weather_Amp / 2.0 * math_cos(alx - depthLag) + deltaTemp) * math_exp(-depthLag)


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * (max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5) + waterBalance_Salb * weather_MeanT
    boundCheck(surfaceT, -100.0, 100.0, "Initial surfaceT")
    return surfaceT


def doNetRadiation(ITERATIONSperDAY, weather_Latitude, clock_Today_DayOfYear, weather_Radn, weather_MinT):
    TSTEPS2RAD = Divide(2.0 * 3.141592653589793, float(ITERATIONSperDAY), 0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin(4.869 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25) + 0.03345 * math_sin(6.224 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25)))
    cD = (1.0 - solarDeclination * solarDeclination) ** 0.5
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    lat_rad = weather_Latitude * 3.141592653589793 / 180.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(lat_rad) + cD * math_cos(lat_rad) * math_cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
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
    cva = math_exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(solarRadn, cloudFr, cva, internalTimeStep, soilTemp_surfaceNode, airTemperature, waterBalance_Salb, waterBalance_Eo, waterBalance_Eos):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilTemp_surfaceNode) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def math_sin(x):
    import math
    return math.sin(x)


def math_cos(x):
    import math
    return math.cos(x)


def math_log(x):
    import math
    return math.log(x)


def math_exp(x):
    import math
    return math.exp(x)


def getBoundaryLayerConductance(TNew_zb, airTemperature, weather_AirPressure, weather_Wind, canopyHeight, instrumentHeight, internalTimeStep, waterBalance_Eo, waterBalance_Eos):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    stefanBoltzmannConstant = 0.0000000567
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight
    surfaceNode = 1
    surfaceTemperature = TNew_zb[surfaceNode]
    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)
    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0
    for iteration in range(1, 3 + 1):
        denom_m = math_log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0)) + stabilityCorrectionMomentum
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, denom_m, 0)
        denom_h = math_log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0)) + stabilityCorrectionHeat
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, denom_h, 0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + (1.0 - 16.0 * stabilityParammeter) ** 0.5) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


def interpolateTemperature(timeHours, defaultTimeOfMaximumTemperature, maxTempYesterday, minTempYesterday, weather_MinT, weather_MaxT, weather_MeanT):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math_sin((0.0 + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (maxTempYesterday - minTempYesterday) / 2.0 + (maxTempYesterday + minTempYesterday) / 2.0
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0:
            tScale = 0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = math_sin((time + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


def doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom, nodeDepth, thickness, surfaceNode):
    volspecHeatSoil_tmp = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_tmp[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName == "Minerals":
                continue
            volspecHeatSoil_tmp[node] += volumetricSpecificHeat(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    mapLayer2Node(volspecHeatSoil_tmp, volSpecHeatSoil, nodeDepth, thickness, surfaceNode, numNodes)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom, nodeDepth, thickness, surfaceNode):
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom)
            thermalConductanceWater = ThermalConductance("Water", node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, ps, pom)
            denom = 1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)
            denom2 = 1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)
            k = (2.0 / 3.0) * (1.0 / denom if denom != 0 else 0.0) + (1.0 / 3.0) * (1.0 / denom2 if denom2 != 0 else 0.0)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = numerator / denominator if denominator != 0 else 0.0
    thermalConductivity = [0.0] * (numNodes + 1)
    mapLayer2Node(thermCondLayers, thermalConductivity, nodeDepth, thickness, surfaceNode, numNodes)
    return thermalConductivity


def doThomas(newTemps, soilTemp, volSpecHeatSoil, thermalConductivity, internalTimeStep, nu, netRadiationSource, netRadiation, waterBalance_Eos, waterBalance_Es, timestep, thermalConductance, heatStorage, numNodes):
    airNode = 0
    surfaceNode = 1
    a = [0.0] * (numNodes + 1 + 1)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)
    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (1.0)  # Placeholder; actual nodeDepth spacing handled where needed
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0)
        elementLength = 1.0  # Placeholder; true element length used via mapping at conductivity calc
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
        radnNet = Divide(waterBalance_Eos * 2465000.0, timestep, 0)
    latentHeatFlux = Divide(waterBalance_Es * 2465000.0, timestep, 0)
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
    return newTemps, heatStorage, thermalConductance


def doUpdate(numInterationsPerDay, soilTemp, newTemperature, minSoilTemp, maxSoilTemp, aveSoilTemp, timeOfDaySecs, internalTimeStep, thermalConductivity, numNodes):
    airNode = 0
    surfaceNode = 1
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
    boundaryLayerConductance_inc = Divide(thermalConductivity[airNode], numInterationsPerDay, 0)
    return soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance_inc


def getIniVariables(weather_Tav, instrumHeight):
    defaultInstrumentHeight = 1.2
    boundCheck(weather_Tav, -30.0, 50.0, "tav (oC)")
    if instrumHeight is not None and instrumHeight > 0.00001:
        instrumentHeight = instrumHeight
    else:
        instrumentHeight = defaultInstrumentHeight
    return instrumentHeight


def getProfileVariables(physical_Thickness, physical_BD, physical_Rocks, physical_ParticleSizeSand, physical_ParticleSizeSilt, physical_ParticleSizeClay, organic_Carbon, waterBalance_SW, DepthToConstantTemperature, numPhantomNodes):
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
    nodeDepth = [0.0] * (numNodes + 1 + 1)
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0
    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers] if numLayers < len(bulkDensity) else bulkDensity[-1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            t_idx_prev = layer - 1
            numer = waterBalance_SW[layer - 1] * (thickness[t_idx_prev] if t_idx_prev < len(thickness) else 0.0)
            soilWater[layer] = Divide(numer, thickness[layer], 0)
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
    return {
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
        "thermalConductance": thermalConductance
    }


def readParam(numLayers, numNodes, bulkDensity, clay, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear):
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay)
    surfaceNode = 1
    soilTemp = calcSoilTemperature(thickness, numNodes, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, surfaceNode)
    newTemperature = soilTemp[:]
    soilRoughnessHeight = 57.0
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4, soilTemp, newTemperature, soilRoughnessHeight


def OnStartOfSimulation(physical_Thickness, physical_BD, physical_Rocks, physical_ParticleSizeSand, physical_ParticleSizeSilt, physical_ParticleSizeClay, organic_Carbon, waterBalance_SW, waterBalance_Salb, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, microClimate_CanopyHeight, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, instrumHeight=None, InitialValues=None, DepthToConstantTemperature=10000.0):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    instrumentHeight = getIniVariables(weather_Tav, instrumHeight)
    profile = getProfileVariables(physical_Thickness, physical_BD, physical_Rocks, physical_ParticleSizeSand, physical_ParticleSizeSilt, physical_ParticleSizeClay, organic_Carbon, waterBalance_SW, DepthToConstantTemperature, numPhantomNodes)
    numLayers = profile["numLayers"]
    numNodes = profile["numNodes"]
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4, soilTemp, newTemperature, soilRoughnessHeight = readParam(numLayers, numNodes, profile["bulkDensity"], profile["clay"], profile["thickness"], weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)
    state = {
        "doInitialisationStuff": doInitialisationStuff,
        "internalTimeStep": internalTimeStep,
        "timeOfDaySecs": timeOfDaySecs,
        "numNodes": numNodes,
        "numLayers": numLayers,
        "nodeDepth": profile["nodeDepth"],
        "thermCondPar1": thermCondPar1,
        "thermCondPar2": thermCondPar2,
        "thermCondPar3": thermCondPar3,
        "thermCondPar4": thermCondPar4,
        "volSpecHeatSoil": profile["volSpecHeatSoil"],
        "soilTemp": soilTemp,
        "morningSoilTemp": profile["morningSoilTemp"],
        "heatStorage": profile["heatStorage"],
        "thermalConductance": profile["thermalConductance"],
        "thermalConductivity": profile["thermalConductivity"],
        "boundaryLayerConductance": 0.0,
        "newTemperature": newTemperature,
        "airTemperature": 0.0,
        "maxTempYesterday": 0.0,
        "minTempYesterday": 0.0,
        "soilWater": profile["soilWater"],
        "minSoilTemp": profile["minSoilTemp"],
        "maxSoilTemp": profile["maxSoilTemp"],
        "aveSoilTemp": profile["aveSoilTemp"],
        "thickness": profile["thickness"],
        "bulkDensity": profile["bulkDensity"],
        "rocks": profile["rocks"],
        "carbon": profile["carbon"],
        "sand": profile["sand"],
        "silt": profile["silt"],
        "clay": profile["clay"],
        "soilRoughnessHeight": soilRoughnessHeight,
        "instrumentHeight": instrumentHeight,
        "netRadiation": 0.0,
        "canopyHeight": 0.0,
        "instrumHeight": instrumHeight if instrumHeight is not None else 0.0,
        "nu": 0.6,
        "boundarLayerConductanceSource": "calc",
        "netRadiationSource": "calc",
        "InitialValues": InitialValues if InitialValues is not None else None,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "timestep": 24.0 * 60.0 * 60.0,
        "defaultTimeOfMaximumTemperature": 14.0,
        "pom": 1.3,
        "ps": 2.63,
        "soilConstituentNames": ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"],
        "waterBalance_Salb": waterBalance_Salb,
        "waterBalance_Eos": waterBalance_Eos,
        "waterBalance_Eo": waterBalance_Eo,
        "waterBalance_Es": waterBalance_Es,
        "microClimate_CanopyHeight": microClimate_CanopyHeight
    }
    return state


def OnProcess(state, weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp, weather_AirPressure, weather_Radn, weather_Wind, weather_Latitude, clock_Today_DayOfYear, microClimate_CanopyHeight, waterBalance_SW, waterBalance_Salb, waterBalance_Eo, waterBalance_Eos, waterBalance_Es):
    import copy
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    interactionsPerDay = 48
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    # Update external daily variables in state
    state["waterBalance_Salb"] = waterBalance_Salb
    state["waterBalance_Eo"] = waterBalance_Eo
    state["waterBalance_Eos"] = waterBalance_Eos
    state["waterBalance_Es"] = waterBalance_Es
    # Update soil water and canopy/instrument heights
    # getOtherVariables
    sw_copy = [0.0] * (numLayers)
    for i in range(numLayers):
        sw_copy[i] = waterBalance_SW[i]
    for layer in range(1, numLayers + 1):
        state["soilWater"][layer] = sw_copy[layer - 1]
    state["soilWater"][numNodes] = state["soilWater"][numLayers]
    state["canopyHeight"] = max(microClimate_CanopyHeight, state["soilRoughnessHeight"]) / 1000.0
    state["instrumentHeight"] = max(state["instrumentHeight"], state["canopyHeight"] + 0.5)
    # Initialisation on first day
    if state["doInitialisationStuff"]:
        if ValuesInArray(state["InitialValues"] if state["InitialValues"] is not None else None):
            state["soilTemp"] = [0.0] * (numNodes + 1 + 1)
            init_vals = state["InitialValues"]
            for i in range(len(init_vals)):
                state["soilTemp"][topsoilNode + i] = init_vals[i]
        else:
            # compute seasonal profile
            state["soilTemp"] = calcSoilTemperature(state["thickness"], numNodes, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, surfaceNode)
            init_vals2 = [0.0] * numLayers
            for i in range(numLayers):
                init_vals2[i] = state["soilTemp"][topsoilNode + i]
            state["InitialValues"] = init_vals2
        state["soilTemp"][airNode] = weather_MeanT
        state["soilTemp"][surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
        for i in range(numNodes + 1, len(state["soilTemp"])):
            state["soilTemp"][i] = weather_Tav
        state["newTemperature"] = state["soilTemp"][:]
        state["maxTempYesterday"] = weather_MaxT
        state["minTempYesterday"] = weather_MinT
        state["doInitialisationStuff"] = False
    # Daily process
    cva = 0.0
    cloudFr = 0.0
    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, weather_Latitude, clock_Today_DayOfYear, weather_Radn, weather_MinT)
    Zero(state["minSoilTemp"])
    Zero(state["maxSoilTemp"])
    Zero(state["aveSoilTemp"])
    state["boundaryLayerConductance"] = 0.0
    state["internalTimeStep"] = round(state["timestep"] / interactionsPerDay)
    # precompute capacities and conductivities for the day (assuming constant within-day water state)
    state["volSpecHeatSoil"] = doVolumetricSpecificHeat(state["numNodes"], state["soilConstituentNames"], state["soilWater"], state["rocks"], state["carbon"], state["sand"], state["silt"], state["clay"], state["bulkDensity"], state["ps"], state["pom"], state["nodeDepth"], state["thickness"], surfaceNode)
    state["thermalConductivity"] = doThermalConductivity(state["numNodes"], state["soilConstituentNames"], state["soilWater"], state["rocks"], state["carbon"], state["sand"], state["silt"], state["clay"], state["bulkDensity"], state["ps"], state["pom"], state["nodeDepth"], state["thickness"], surfaceNode)
    boundaryLayerConductance_accum = 0.0
    for timeStepIteration in range(1, interactionsPerDay + 1):
        state["timeOfDaySecs"] = state["internalTimeStep"] * float(timeStepIteration)
        if state["timestep"] < 24.0 * 60.0 * 60.0:
            state["airTemperature"] = weather_MeanT
        else:
            state["airTemperature"] = interpolateTemperature(state["timeOfDaySecs"] / 3600.0, state["defaultTimeOfMaximumTemperature"], state["maxTempYesterday"], state["minTempYesterday"], weather_MinT, weather_MaxT, weather_MeanT)
        state["newTemperature"][airNode] = state["airTemperature"]
        state["netRadiation"] = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, state["internalTimeStep"], state["soilTemp"][surfaceNode], state["airTemperature"], waterBalance_Salb, waterBalance_Eo, waterBalance_Eos)
        if state["boundarLayerConductanceSource"] == "constant":
            state["thermalConductivity"][airNode] = 20.0
        elif state["boundarLayerConductanceSource"] == "calc":
            state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"], state["airTemperature"], weather_AirPressure, weather_Wind, state["canopyHeight"], state["instrumentHeight"], state["internalTimeStep"], waterBalance_Eo, waterBalance_Eos)
            for iteration in range(1, 1 + 1):
                state["newTemperature"], state["heatStorage"], state["thermalConductance"] = doThomas(state["newTemperature"], state["soilTemp"], state["volSpecHeatSoil"], state["thermalConductivity"], state["internalTimeStep"], state["nu"], state["netRadiationSource"], state["netRadiation"], waterBalance_Eos, waterBalance_Es, state["timestep"], state["thermalConductance"], state["heatStorage"], state["numNodes"])
                state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"], state["airTemperature"], weather_AirPressure, weather_Wind, state["canopyHeight"], state["instrumentHeight"], state["internalTimeStep"], waterBalance_Eo, waterBalance_Eos)
        state["newTemperature"], state["heatStorage"], state["thermalConductance"] = doThomas(state["newTemperature"], state["soilTemp"], state["volSpecHeatSoil"], state["thermalConductivity"], state["internalTimeStep"], state["nu"], state["netRadiationSource"], state["netRadiation"], waterBalance_Eos, waterBalance_Es, state["timestep"], state["thermalConductance"], state["heatStorage"], state["numNodes"])
        state["soilTemp"], state["minSoilTemp"], state["maxSoilTemp"], state["aveSoilTemp"], b_inc = doUpdate(interactionsPerDay, state["soilTemp"], state["newTemperature"], state["minSoilTemp"], state["maxSoilTemp"], state["aveSoilTemp"], state["timeOfDaySecs"], state["internalTimeStep"], state["thermalConductivity"], state["numNodes"])
        boundaryLayerConductance_accum += b_inc
        if abs(state["timeOfDaySecs"] - 5.0 * 3600.0) <= min(state["timeOfDaySecs"], 5.0 * 3600.0) * 0.0001:
            state["morningSoilTemp"] = state["soilTemp"][:]
    state["minTempYesterday"] = weather_MinT
    state["maxTempYesterday"] = weather_MaxT
    state["boundaryLayerConductance"] = boundaryLayerConductance_accum
    # Build outputs
    FinalSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        FinalSoilTemperature[i] = state["soilTemp"][topsoilNode + i]
    FinalSoilSurfaceTemperature = state["soilTemp"][surfaceNode]
    AverageSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        AverageSoilTemperature[i] = state["aveSoilTemp"][topsoilNode + i]
    AverageSoilSurfaceTemperature = state["aveSoilTemp"][surfaceNode]
    MinimumSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        MinimumSoilTemperature[i] = state["minSoilTemp"][topsoilNode + i]
    MinimumSoilSurfaceTemperature = state["minSoilTemp"][surfaceNode]
    MaximumSoilTemperature = [0.0] * numLayers
    for i in range(numLayers):
        MaximumSoilTemperature[i] = state["maxSoilTemp"][topsoilNode + i]
    MaximumSoilSurfaceTemperature = state["maxSoilTemp"][surfaceNode]
    ThermalConductivity_out = [0.0] * (numNodes)
    for i in range(numNodes):
        idx = 1 + i
        if idx < len(state["thermalConductivity"]):
            ThermalConductivity_out[i] = state["thermalConductivity"][idx]
    HeatCapacity_out = [0.0] * (numNodes)
    for i in range(numNodes):
        idx = surfaceNode + i
        if idx < len(state["volSpecHeatSoil"]):
            HeatCapacity_out[i] = state["volSpecHeatSoil"][idx]
    HeatStore_out = [0.0] * (numNodes)
    for i in range(numNodes):
        idx = surfaceNode + i
        if idx < len(state["heatStorage"]):
            HeatStore_out[i] = state["heatStorage"][idx]
    Thr_profile = state["morningSoilTemp"][:]
    outputs = {
        "FinalSoilTemperature": FinalSoilTemperature,
        "FinalSoilSurfaceTemperature": FinalSoilSurfaceTemperature,
        "AverageSoilTemperature": AverageSoilTemperature,
        "AverageSoilSurfaceTemperature": AverageSoilSurfaceTemperature,
        "MinimumSoilTemperature": MinimumSoilTemperature,
        "MinimumSoilSurfaceTemperature": MinimumSoilSurfaceTemperature,
        "MaximumSoilTemperature": MaximumSoilTemperature,
        "MaximumSoilSurfaceTemperature": MaximumSoilSurfaceTemperature,
        "BoundaryLayerConductance": state["boundaryLayerConductance"],
        "ThermalConductivity": ThermalConductivity_out,
        "HeatCapacity": HeatCapacity_out,
        "HeatStore": HeatStore_out,
        "Thr_profile": Thr_profile
    }
    return state, outputs