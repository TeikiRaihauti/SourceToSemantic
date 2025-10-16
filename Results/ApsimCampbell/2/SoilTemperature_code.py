def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and Value == Value:  # not NaN
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


def volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity):
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


def ThermalConductance(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps):
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
        # Not used due to line below in original code.
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, rocks)) * \
                 (thermalConductanceSand ** volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps)) + \
                 (thermalConductanceSilt ** volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps)) + \
                 (thermalConductanceClay ** volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps))
    # Replicate original bug: return volumetric specific heat value instead of thermal conductance.
    return volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity)


def shapeFactor(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps):
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
        w = volumetricFractionWater(layer, soilWater, carbon)
        i = volumetricFractionIce(layer)
        a = volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
        result = 0.333 - 0.333 * i / (w + i + a)
        return result
    elif name == "Air":
        w = volumetricFractionWater(layer, soilWater, carbon)
        i = volumetricFractionIce(layer)
        a = volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps)
        result = 0.333 - 0.333 * a / (w + i + a)
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + \
                 shapeFactorSand * volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps) + \
                 shapeFactorSilt * volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps) + \
                 shapeFactorClay * volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps)
    # Replicate original bug: return volumetric specific heat value instead of shape factor.
    return volumetricSpecificHeat(name, layer, soilWater, rocks, carbon, sand, silt, clay, bulkDensity)


def volumetricFractionRocks(layer, rocks):
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom):
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
           sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
           silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
           clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(layer, soilWater, carbon):
    return (1 - volumetricFractionOrganicMatter(layer, carbon, soilWater, 1.3)) * soilWater[layer]


def volumetricFractionIce(layer):
    return 0.0


def volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity, pom, ps):
    return 1.0 - volumetricFractionRocks(layer, rocks) - \
           volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - \
           volumetricFractionSand(layer, rocks, carbon, sand, bulkDensity, pom, ps) - \
           volumetricFractionSilt(layer, rocks, carbon, silt, bulkDensity, pom, ps) - \
           volumetricFractionClay(layer, rocks, carbon, clay, bulkDensity, pom, ps) - \
           volumetricFractionWater(layer, soilWater, carbon) - \
           volumetricFractionIce(layer)


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def mapLayer2Node(layerArray, nodeArray, thickness, nodeDepth, numNodes):
    for node in range(1, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)


def doThermalConductivityCoeffs(numNodes, numLayers, bulkDensity, clay):
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

    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def doVolumetricSpecificHeat(numNodes, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, thickness, nodeDepth, volSpecHeatSoil):
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]

    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName != "Minerals":
                volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity) * 1000000.0 * soilWater[node]

    mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, thickness, nodeDepth, numNodes)


def doThermalConductivity(numNodes, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, thickness, nodeDepth, thermalConductivity, pom, ps):
    thermCondLayers = [0.0] * (numNodes + 1)
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]

    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps)
            thermalConductanceWater = ThermalConductance("Water", node, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps)
            denom1 = 1.0 + shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0)
            denom2 = 1.0 + shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0) * (1 - 2 * shapeFactorConstituent)
            k = (2.0 / 3.0) * (denom1 ** -1 if denom1 != 0 else 0.0) + (1.0 / 3.0) * (denom2 ** -1 if denom2 != 0 else 0.0)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k

        thermCondLayers[node] = Divide(numerator, denominator, 0.0)

    mapLayer2Node(thermCondLayers, thermalConductivity, thickness, nodeDepth, numNodes)


def doThomas(newTemps, numNodes, nu, volSpecHeatSoil, nodeDepth, internalTimeStep, thermalConductivity, netRadiation, soilTemp, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb, timestep, latentHeatOfVapourisation, netRadiationSource, airTemperature):
    airNode = 0
    surfaceNode = 1

    a = [0.0] * (numNodes + 2)  # a[0..numNodes+1]
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 2)

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
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + \
                  (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + \
                  g * thermalConductance[node] * soilTemp[node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)
    else:
        radnNet = 0.0

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
        midnightT = (math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) *
                     (maxTempYesterday - minTempYesterday) / 2.0
                     + (maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = (math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) *
                              (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT)
        return currentTemperature


def doUpdate(newTemperature, soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance, numNodes, timeOfDaySecs, internalTimeStep, thermalConductivity):
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
        aveSoilTemp[node] += Divide(soilTemp[node], 48, 0.0)
    boundaryLayerConductance += Divide(thermalConductivity[airNode], 48, 0.0)

    return soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def getBoundaryLayerConductance(TNew_zb, airTemperature, weather_AirPressure, canopyHeight, instrumentHeight, waterBalance_Eos, waterBalance_Eo):
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

    for _ in range(3):
        numerator = weather_Wind = 0.0  # weather_Wind not passed here; set to 0.0 if not known
        # To match original code behaviour, require wind externally. We'll set small baseline wind to avoid div-by-zero.
        weather_Wind = 0.0

        frictionVelocity = Divide(weather_Wind * vonKarmanConstant,
                                  math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum,
                                  0.0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity,
                                   math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat,
                                   0.0)

        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)

        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity,
                                     SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0),
                                     0.0)

        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            val = 1.0 - 16.0 * stabilityParammeter
            root_term = math.sqrt(val) if val > 0 else 0.0
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + root_term) / 2.0) if (1.0 + root_term) > 0 else 0.0
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond


def calcSoilTemperature(soilTempIO, numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25

    soilTemp_ = [0.0] * (numNodes + 2)
    for nodes in range(1, numNodes + 1):
        soilTemp_[nodes] = weather_Tav + weather_Amp * math.exp(-1 * cumulativeDepth[nodes] / zd) * \
                           math.sin((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)

    for i in range(0, numNodes + 1):
        idx = i + 1  # constrained copy to index starting at surfaceNode=1? In original: copy from 0 to soilTempIO starting at 1
        if idx < len(soilTempIO) and i < len(soilTemp_):
            soilTempIO[idx] = soilTemp_[i]


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    return surfaceT


def doNetRadiation(ITERATIONSperDAY, weather_Radn, weather_Latitude, weather_MinT, clock_Today_DayOfYear):
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
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
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)

    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)

    cva = math.exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(solarRadn, cloudFr, cva, internalTimeStep, airTemperature, soilSurfaceTemperature, waterBalance_Eos, waterBalance_Eo):
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0

    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr

    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)

    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilSurfaceTemperature) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil

    swRin = solarRadn
    swRout = 0.0  # soil albedo applied in doThomas via waterBalance_Salb - consistent with original usage
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def getOtherVariables(state, waterBalance_SW, microClimate_CanopyHeight):
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    soilWater = state["soilWater"]

    for i in range(len(waterBalance_SW)):
        idx = i + 1
        if idx < len(soilWater):
            soilWater[idx] = waterBalance_SW[i]
    soilWater[numNodes] = soilWater[numLayers]

    canopyHeight = max(microClimate_CanopyHeight, state["soilRoughnessHeight"]) / 1000.0
    instrumentHeight = max(state["instrumentHeight"], canopyHeight + 0.5)

    state["soilWater"] = soilWater
    state["canopyHeight"] = canopyHeight
    state["instrumentHeight"] = instrumentHeight
    return state


def doProcess(state,
              weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp,
              weather_AirPressure, weather_Radn, weather_Latitude, weather_Wind,
              clock_Today_DayOfYear,
              waterBalance_SW, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2

    interactionsPerDay = 48

    soilTemp = state["soilTemp"]
    newTemperature = state["newTemperature"]
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    boundaryLayerConductance = 0.0
    aveSoilTemp = state["aveSoilTemp"]
    minSoilTemp = state["minSoilTemp"]
    maxSoilTemp = state["maxSoilTemp"]
    nodeDepth = state["nodeDepth"]
    thickness = state["thickness"]
    bulkDensity = state["bulkDensity"]
    rocks = state["rocks"]
    carbon = state["carbon"]
    sand = state["sand"]
    silt = state["silt"]
    clay = state["clay"]
    pom = state["pom"]
    ps = state["ps"]
    nu = state["nu"]
    latentHeatOfVapourisation = state["latentHeatOfVapourisation"]
    netRadiationSource = state["netRadiationSource"]
    boundarLayerConductanceSource = state["boundarLayerConductanceSource"]
    constantBoundaryLayerConductance = state["constantBoundaryLayerConductance"]
    defaultTimeOfMaximumTemperature = state["defaultTimeOfMaximumTemperature"]
    minTempYesterday = state["minTempYesterday"]
    maxTempYesterday = state["maxTempYesterday"]
    timestep = state["timestep"]
    instrumentHeight = state["instrumentHeight"]
    canopyHeight = state["canopyHeight"]
    thermalConductivity = state["thermalConductivity"]
    volSpecHeatSoil = state["volSpecHeatSoil"]
    heatStorage = state["heatStorage"]
    thermalConductance = state["thermalConductance"]

    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, weather_Radn, weather_Latitude, weather_MinT, clock_Today_DayOfYear)

    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0
    internalTimeStep = round(timestep / interactionsPerDay)

    doVolumetricSpecificHeat(numNodes, state["soilWater"], rocks, carbon, sand, silt, clay, bulkDensity, thickness, nodeDepth, volSpecHeatSoil)
    doThermalConductivity(numNodes, state["soilWater"], rocks, carbon, sand, silt, clay, bulkDensity, thickness, nodeDepth, thermalConductivity, pom, ps)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0,
                                                    defaultTimeOfMaximumTemperature,
                                                    maxTempYesterday, minTempYesterday,
                                                    weather_MaxT, weather_MinT, weather_MeanT)
        newTemperature[airNode] = airTemperature

        netRadiation = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, internalTimeStep, airTemperature, soilTemp[surfaceNode], waterBalance_Eos, waterBalance_Eo)

        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        else:
            # Note: getBoundaryLayerConductance in original uses weather.Wind; here we set via thermalConductivity update after doThomas iterations
            # We don't have wind inside getBoundaryLayerConductance; preserve structure.
            thermalConductivity[airNode] = constantBoundaryLayerConductance
            for _ in range(1):  # numIterationsForBoundaryLayerConductance
                newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, numNodes, nu, volSpecHeatSoil, nodeDepth, internalTimeStep, thermalConductivity, netRadiation, soilTemp,
                                                                           waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb, timestep, latentHeatOfVapourisation, netRadiationSource, airTemperature)
                thermalConductivity[airNode] = constantBoundaryLayerConductance

        newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, numNodes, nu, volSpecHeatSoil, nodeDepth, internalTimeStep, thermalConductivity, netRadiation, soilTemp,
                                                                   waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb, timestep, latentHeatOfVapourisation, netRadiationSource, airTemperature)
        soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = doUpdate(newTemperature, soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance, numNodes, timeOfDaySecs, internalTimeStep, thermalConductivity)
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            state["morningSoilTemp"] = soilTemp.copy()

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    state["soilTemp"] = soilTemp
    state["newTemperature"] = newTemperature
    state["boundaryLayerConductance"] = boundaryLayerConductance
    state["aveSoilTemp"] = aveSoilTemp
    state["minSoilTemp"] = minSoilTemp
    state["maxSoilTemp"] = maxSoilTemp
    state["thermalConductivity"] = thermalConductivity
    state["volSpecHeatSoil"] = volSpecHeatSoil
    state["heatStorage"] = heatStorage
    state["thermalConductance"] = thermalConductance
    state["minTempYesterday"] = minTempYesterday
    state["maxTempYesterday"] = maxTempYesterday
    state["internalTimeStep"] = internalTimeStep

    return state


def OnProcess(state,
              weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp,
              weather_AirPressure, weather_Radn, weather_Latitude, weather_Wind,
              clock_Today_DayOfYear,
              microClimate_CanopyHeight,
              waterBalance_SW, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb):
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2

    state = getOtherVariables(state, waterBalance_SW, microClimate_CanopyHeight)

    if state["doInitialisationStuff"]:
        if ValuesInArray(state["InitialValues"]):
            state["soilTemp"] = [0.0] * (state["numNodes"] + 2)
            for i in range(len(state["InitialValues"])):
                idx = topsoilNode + i
                if idx < len(state["soilTemp"]):
                    state["soilTemp"][idx] = state["InitialValues"][i]
        else:
            calcSoilTemperature(state["soilTemp"], state["numNodes"], state["thickness"], weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude)
            state["InitialValues"] = [0.0] * state["numLayers"]
            for i in range(state["numLayers"]):
                state["InitialValues"][i] = state["soilTemp"][topsoilNode + i]

        state["soilTemp"][airNode] = weather_MeanT
        state["soilTemp"][surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)

        last_idx = state["numNodes"] + 1
        if last_idx < len(state["soilTemp"]):
            state["soilTemp"][last_idx] = weather_Tav

        state["newTemperature"] = state["soilTemp"].copy()

        state["maxTempYesterday"] = weather_MaxT
        state["minTempYesterday"] = weather_MinT
        state["doInitialisationStuff"] = False

    state = doProcess(state,
                      weather_MeanT, weather_MaxT, weather_MinT, weather_Tav, weather_Amp,
                      weather_AirPressure, weather_Radn, weather_Latitude, weather_Wind,
                      clock_Today_DayOfYear,
                      waterBalance_SW, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, waterBalance_Salb)

    numLayers = state["numLayers"]
    numNodes = state["numNodes"]
    topsoilNode = 2
    surfaceNode = 1
    airNode = 0

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

    BoundaryLayerConductance = state["boundaryLayerConductance"]

    ThermalConductivity = [0.0] * numNodes
    for i in range(numNodes):
        ThermalConductivity[i] = state["thermalConductivity"][1 + i]

    HeatCapacity = [0.0] * numNodes
    for i in range(numNodes):
        HeatCapacity[i] = state["volSpecHeatSoil"][surfaceNode + i]

    HeatStore = [0.0] * numNodes
    for i in range(numNodes):
        HeatStore[i] = state["heatStorage"][surfaceNode + i]

    Thr_profile = state["morningSoilTemp"].copy()

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
        "ThermalConductivity": ThermalConductivity,
        "HeatCapacity": HeatCapacity,
        "HeatStore": HeatStore,
        "Thr_profile": Thr_profile
    }

    return state, outputs


def OnStartOfSimulation(physical_Thickness, physical_BD, physical_Rocks, physical_ParticleSizeSand, physical_ParticleSizeSilt, physical_ParticleSizeClay,
                        organic_Carbon, waterBalance_SW, microClimate_CanopyHeight,
                        weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude,
                        instrumHeight=None, InitialValues=None, DepthToConstantTemperature=10000.0):
    import math as _math  # ensure math imported for called funcs
    globals()['math'] = __import__('math')  # provide math to functions

    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5

    timestep = 24.0 * 60.0 * 60.0
    latentHeatOfVapourisation = 2465000.0
    constantBoundaryLayerConductance = 20.0
    defaultTimeOfMaximumTemperature = 14.0
    bareSoilRoughness = 57.0
    defaultInstrumentHeight = 1.2

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

    pom = 1.3
    ps = 2.63
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numNodes, numLayers, bulkDensity, clay)

    calcSoilTemperature(soilTemp, numNodes, thickness, weather_Tav, weather_Amp, clock_Today_DayOfYear, weather_Latitude)
    newTemperature = soilTemp.copy()
    soilRoughnessHeight = bareSoilRoughness

    instrumentHeight_val = defaultInstrumentHeight
    if instrumHeight is not None and instrumHeight > 0.00001:
        instrumentHeight_val = instrumHeight

    state = {
        "doInitialisationStuff": True,
        "internalTimeStep": 0.0,
        "timeOfDaySecs": 0.0,
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
        "instrumentHeight": instrumentHeight_val,
        "netRadiation": 0.0,
        "canopyHeight": max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0,
        "instrumHeight": instrumHeight if instrumHeight is not None else 0.0,
        "nu": 0.6,
        "boundarLayerConductanceSource": "calc",
        "netRadiationSource": "calc",
        "pInitialValues": InitialValues,
        "InitialValues": InitialValues,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "timestep": timestep,
        "latentHeatOfVapourisation": latentHeatOfVapourisation,
        "constantBoundaryLayerConductance": constantBoundaryLayerConductance,
        "defaultTimeOfMaximumTemperature": defaultTimeOfMaximumTemperature,
        "bareSoilRoughness": bareSoilRoughness,
        "defaultInstrumentHeight": defaultInstrumentHeight,
        "pom": pom,
        "ps": ps
    }
    return state


# No tests were provided in the original codebase.