def initialize_SoilTemperature(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    weather_Tav,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Amp,
    weather_Latitude,
    weather_Radn,
    weather_AirPressure,
    weather_Wind,
    clock_Today_DayOfYear,
    InitialValues=None,
    DepthToConstantTemperature=10000.0,
    instrumHeight=0.0
):
    # constants
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    defaultInstrumentHeight = 1.2
    bareSoilRoughness = 57.0

    def Sum(values, startIndex, endIndex):
        result = 0.0
        index = -1
        for value in values:
            index += 1
            if index >= startIndex:
                result += value
            if index == endIndex:
                break
        return result

    def ToCumThickness(Thickness):
        CumThickness = [0.0] * len(Thickness)
        if len(Thickness) > 0:
            CumThickness[0] = Thickness[0]
            for Layer in range(1, len(Thickness)):
                CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
        return CumThickness

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

    # profile setup
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
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * (thickness[layer - 1] if (layer - 1) < len(thickness) else 0.0), thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]

    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        carbon[i + 1] = organic_Carbon[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]

    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        rocks[i + 1] = physical_Rocks[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]

    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        sand[i + 1] = physical_ParticleSizeSand[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]

    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        silt[i + 1] = physical_ParticleSizeSilt[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]

    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        clay[i + 1] = physical_ParticleSizeClay[i]
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

    soilRoughnessHeight = bareSoilRoughness
    instrumentHeight = instrumHeight if instrumHeight > 1e-5 else defaultInstrumentHeight

    def ToCumThickness_nodes(thickness_arr):
        return ToCumThickness(thickness_arr)

    def calcSoilTemperature(soilTempIO):
        cumulativeDepth = ToCumThickness_nodes(thickness)
        import math
        w = 2 * math.pi / (365.25 * 24 * 3600)
        dh = 0.6
        zd = (2 * dh / w) ** 0.5
        offset = 0.25
        if weather_Latitude > 0.0:
            offset = -0.25
        tmp = [0.0] * (numNodes + 2)
        for nodes in range(1, numNodes + 1):
            tmp[nodes] = weather_Tav + weather_Amp * math.exp(-1.0 * cumulativeDepth[nodes] / zd) * math.sin(((clock_Today_DayOfYear / 365.0) + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)
        # copy into soilTempIO starting at surfaceNode
        for i in range(0, numNodes):
            dst = surfaceNode + i
            if dst < len(soilTempIO) and i < len(tmp):
                soilTempIO[dst] = tmp[i]

    if InitialValues is not None and ValuesInArray(InitialValues):
        for i in range(len(InitialValues)):
            idx = topsoilNode + i
            if idx < len(soilTemp):
                soilTemp[idx] = InitialValues[i]
    else:
        calcSoilTemperature(soilTemp)

    def calcSurfaceTemperature(waterBalance_Salb_, weather_MeanT_, weather_MaxT_, weather_Radn_):
        surfaceT = (1.0 - waterBalance_Salb_) * (weather_MeanT_ + (weather_MaxT_ - weather_MeanT_) * ((max(weather_Radn_, 0.1) * 23.8846 / 800.0) ** 0.5)) + waterBalance_Salb_ * weather_MeanT_
        return surfaceT

    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav
    for i in range(len(newTemperature)):
        newTemperature[i] = soilTemp[i]

    state = {
        "timestep": 24.0 * 60.0 * 60.0,
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
        "canopyHeight": max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0,
        "nu": 0.6,
        "boundarLayerConductanceSource": "calc",
        "netRadiationSource": "calc",
        "defaultTimeOfMaximumTemperature": 14.0
    }
    return state


def process_SoilTemperature(
    state,
    waterBalance_SW,
    waterBalance_Eo,
    waterBalance_Eos,
    waterBalance_Es,
    waterBalance_Salb,
    microClimate_CanopyHeight,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Radn,
    weather_AirPressure,
    weather_Wind,
    weather_Tav,
    weather_Latitude,
    clock_Today_DayOfYear
):
    def Divide(value1, value2, errVal):
        if value2 != 0:
            return value1 / value2
        return errVal

    def Zero(arr):
        if arr is not None:
            for i in range(len(arr)):
                arr[i] = 0.0

    def Sum(values, startIndex, endIndex):
        result = 0.0
        index = -1
        for value in values:
            index += 1
            if index >= startIndex:
                result += value
            if index == endIndex:
                break
        return result

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

    airNode = 0
    surfaceNode = 1
    topsoilNode = 2

    # Update other variables (soil water etc)
    if waterBalance_SW is not None:
        for i in range(len(waterBalance_SW)):
            layer = 1 + i
            if layer < len(state["soilWater"]):
                state["soilWater"][layer] = waterBalance_SW[i]
    state["soilWater"][state["numNodes"]] = state["soilWater"][state["numLayers"]]
    state["canopyHeight"] = max(microClimate_CanopyHeight, state["soilRoughnessHeight"]) / 1000.0
    state["instrumentHeight"] = max(state["instrumentHeight"], state["canopyHeight"] + 0.5)

    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn = [0.0] * (interactionsPerDay + 1)

    def doNetRadiation(solarRadn_, ITERATIONSperDAY):
        import math
        TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
        solarConstant = 1360.0
        solarDeclination = 0.3985 * math.sin(4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + 0.03345 * math.sin(6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))
        cD = (1.0 - solarDeclination * solarDeclination) ** 0.5
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
        cloudFr_ = 2.33 - 3.33 * fr
        cloudFr_ = min(max(cloudFr_, 0.0), 1.0)
        for timestepNumber in range(1, ITERATIONSperDAY + 1):
            solarRadn_[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
        cva_ = math.exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
        return cloudFr_, cva_

    cloudFr, cva = doNetRadiation(solarRadn, interactionsPerDay)

    Zero(state["minSoilTemp"])
    Zero(state["maxSoilTemp"])
    Zero(state["aveSoilTemp"])
    state["boundaryLayerConductance"] = 0.0

    state["internalTimeStep"] = round(state["timestep"] / interactionsPerDay)

    def volumetricSpecificHeat(name, layer):
        specificHeatRocks = 7.7
        specificHeatOM = 0.25
        specificHeatSand = 7.7
        specificHeatSilt = 2.74
        specificHeatClay = 2.92
        specificHeatWater = 0.57
        specificHeatIce = 2.18
        specificHeatAir = 0.025
        if name == "Rocks":
            return specificHeatRocks
        elif name == "OrganicMatter":
            return specificHeatOM
        elif name == "Sand":
            return specificHeatSand
        elif name == "Silt":
            return specificHeatSilt
        elif name == "Clay":
            return specificHeatClay
        elif name == "Water":
            return specificHeatWater
        elif name == "Ice":
            return specificHeatIce
        elif name == "Air":
            return specificHeatAir
        return 0.0

    def volumetricFractionRocks(layer):
        return state["rocks"][layer] / 100.0

    def volumetricFractionOrganicMatter(layer):
        pom = 1.3
        return state["carbon"][layer] / 100.0 * 2.5 * state["bulkDensity"][layer] / pom

    def volumetricFractionSand(layer):
        ps = 2.63
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * state["sand"][layer] / 100.0 * state["bulkDensity"][layer] / ps

    def volumetricFractionSilt(layer):
        ps = 2.63
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * state["silt"][layer] / 100.0 * state["bulkDensity"][layer] / ps

    def volumetricFractionClay(layer):
        ps = 2.63
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * state["clay"][layer] / 100.0 * state["bulkDensity"][layer] / ps

    def volumetricFractionWater(layer):
        return (1 - volumetricFractionOrganicMatter(layer)) * (state["soilWater"][layer] if layer < len(state["soilWater"]) else 0.0)

    def volumetricFractionIce(layer):
        return 0.0

    def volumetricFractionAir(layer):
        return 1.0 - volumetricFractionRocks(layer) - volumetricFractionOrganicMatter(layer) - volumetricFractionSand(layer) - volumetricFractionSilt(layer) - volumetricFractionClay(layer) - volumetricFractionWater(layer) - volumetricFractionIce(layer)

    def ThermalConductance(name, layer):
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
            result = (thermalConductanceRocks ** volumetricFractionRocks(layer)) * (thermalConductanceSand ** volumetricFractionSand(layer)) + (thermalConductanceSilt ** volumetricFractionSilt(layer)) + (thermalConductanceClay ** volumetricFractionClay(layer))
        else:
            result = 0.0
        result = volumetricSpecificHeat(name, layer)
        return result

    def shapeFactor(name, layer):
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
            result = 0.333 - 0.333 * volumetricFractionIce(layer) / (volumetricFractionWater(layer) + volumetricFractionIce(layer) + volumetricFractionAir(layer))
            return result
        elif name == "Air":
            result = 0.333 - 0.333 * volumetricFractionAir(layer) / (volumetricFractionWater(layer) + volumetricFractionIce(layer) + volumetricFractionAir(layer))
            return result
        elif name == "Minerals":
            result = shapeFactorRocks * volumetricFractionRocks(layer) + shapeFactorSand * volumetricFractionSand(layer) + shapeFactorSilt * volumetricFractionSilt(layer) + shapeFactorClay * volumetricFractionClay(layer)
        else:
            result = 0.0
        result = volumetricSpecificHeat(name, layer)
        return result

    def mapLayer2Node(layerArray, nodeArray):
        for node in range(surfaceNode, state["numNodes"] + 1):
            layer = node - 1
            depthLayerAbove = Sum(state["thickness"], 1, layer) if layer >= 1 else 0.0
            d1 = depthLayerAbove - (state["nodeDepth"][node] * 1000.0)
            d2 = state["nodeDepth"][node + 1] * 1000.0 - depthLayerAbove
            dSum = d1 + d2
            nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)

    def doVolumetricSpecificHeat():
        volspecHeatSoil_ = [0.0] * (state["numNodes"] + 1)
        soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
        for node in range(1, state["numNodes"] + 1):
            volspecHeatSoil_[node] = 0.0
            for constituentName in soilConstituentNames:
                volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * (state["soilWater"][node] if node < len(state["soilWater"]) else 0.0)
        mapLayer2Node(volspecHeatSoil_, state["volSpecHeatSoil"])

    def doThermalConductivity():
        thermCondLayers = [0.0] * (state["numNodes"] + 1)
        soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
        for node in range(1, state["numNodes"] + 1):
            numerator = 0.0
            denominator = 0.0
            for constituentName in soilConstituentNames:
                shapeFactorConstituent = shapeFactor(constituentName, node)
                thermalConductanceConstituent = ThermalConductance(constituentName, node)
                thermalConductanceWater = ThermalConductance("Water", node)
                k = (2.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** -1) + (1.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** -1)
                sw = state["soilWater"][node] if node < len(state["soilWater"]) else 0.0
                numerator += thermalConductanceConstituent * sw * k
                denominator += sw * k
            thermCondLayers[node] = Divide(numerator, denominator, 0.0)
        mapLayer2Node(thermCondLayers, state["thermalConductivity"])

    def interpolateTemperature(timeHours):
        import math
        time = timeHours / 24.0
        maxT_time = state["defaultTimeOfMaximumTemperature"] / 24.0
        minT_time = maxT_time - 0.5
        if time < minT_time:
            midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (state["maxTempYesterday"] - state["minTempYesterday"]) / 2.0 + (state["maxTempYesterday"] + state["minTempYesterday"]) / 2.0
            tScale = (minT_time - time) / minT_time
            if tScale > 1.0:
                tScale = 1.0
            elif tScale < 0:
                tScale = 0.0
            currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
            return currentTemperature
        else:
            import math
            currentTemperature = math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
            return currentTemperature

    def interpolateNetRadiation(solarRadn_, cloudFr_, cva_):
        surfaceEmissivity = 0.96
        w2MJ = state["internalTimeStep"] / 1000000.0
        emissivityAtmos = (1 - 0.84 * cloudFr_) * 0.58 * (cva_ ** (1.0 / 7.0)) + 0.84 * cloudFr_
        PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
        lwRinSoil = longWaveRadn(emissivityAtmos, state["airTemperature"]) * PenetrationConstant * w2MJ
        lwRoutSoil = longWaveRadn(surfaceEmissivity, state["soilTemp"][surfaceNode]) * PenetrationConstant * w2MJ
        lwRnetSoil = lwRinSoil - lwRoutSoil
        swRin = solarRadn_
        swRout = waterBalance_Salb * solarRadn_
        swRnetSoil = (swRin - swRout) * PenetrationConstant
        return swRnetSoil + lwRnetSoil

    def getBoundaryLayerConductance(TNew_zb):
        vonKarmanConstant = 0.41
        gravitationalConstant = 9.8
        specificHeatOfAir = 1010.0
        stefanBoltzmannConstant = 0.0000000567
        surfaceEmissivity = 0.98
        SpecificHeatAir = specificHeatOfAir * airDensity(state["airTemperature"], weather_AirPressure)
        roughnessFactorMomentum = 0.13 * state["canopyHeight"]
        roughnessFactorHeat = 0.2 * roughnessFactorMomentum
        d = 0.77 * state["canopyHeight"]
        surfaceTemperature = TNew_zb[surfaceNode]
        diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
        radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(state["airTemperature"]) ** 3)
        frictionVelocity = 0.0
        boundaryLayerCond = 0.0
        stabilityParammeter = 0.0
        stabilityCorrectionMomentum = 0.0
        stabilityCorrectionHeat = 0.0
        heatFluxDensity = 0.0
        for _ in range(1, 4):
            frictionVelocity = Divide(weather_Wind * vonKarmanConstant, (math_log(Divide(state["instrumentHeight"] - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum), 0.0)
            boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, (math_log(Divide(state["instrumentHeight"] - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat), 0.0)
            boundaryLayerCond += radiativeConductance
            heatFluxDensity = boundaryLayerCond * (surfaceTemperature - state["airTemperature"])
            stabilityParammeter = Divide(-vonKarmanConstant * state["instrumentHeight"] * gravitationalConstant * heatFluxDensity, (SpecificHeatAir * kelvinT(state["airTemperature"]) * (frictionVelocity ** 3.0)), 0.0)
            if stabilityParammeter > 0.0:
                stabilityCorrectionHeat = 4.7 * stabilityParammeter
                stabilityCorrectionMomentum = stabilityCorrectionHeat
            else:
                # unstable
                stabilityCorrectionHeat = -2.0 * math_log((1.0 + ((1.0 - 16.0 * stabilityParammeter) ** 0.5)) / 2.0)
                stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
        return boundaryLayerCond

    def math_log(x):
        import math
        return math.log(x) if x > 0 else 0.0

    def doThomas(newTemps):
        a = [0.0] * (state["numNodes"] + 2)
        b = [0.0] * (state["numNodes"] + 1)
        c = [0.0] * (state["numNodes"] + 1)
        d = [0.0] * (state["numNodes"] + 1)
        state["thermalConductance"][airNode] = state["thermalConductivity"][airNode]
        for node in range(surfaceNode, state["numNodes"] + 1):
            volumeOfSoilAtNode = 0.5 * (state["nodeDepth"][node + 1] - state["nodeDepth"][node - 1])
            state["heatStorage"][node] = Divide(state["volSpecHeatSoil"][node] * volumeOfSoilAtNode, state["internalTimeStep"], 0.0)
            elementLength = state["nodeDepth"][node + 1] - state["nodeDepth"][node]
            state["thermalConductance"][node] = Divide(state["thermalConductivity"][node], elementLength, 0.0)
        g = 1.0 - state["nu"]
        for node in range(surfaceNode, state["numNodes"] + 1):
            c[node] = (-state["nu"]) * state["thermalConductance"][node]
            a[node + 1] = c[node]
            b[node] = state["nu"] * (state["thermalConductance"][node] + state["thermalConductance"][node - 1]) + state["heatStorage"][node]
            d[node] = g * state["thermalConductance"][node - 1] * state["soilTemp"][node - 1] + (state["heatStorage"][node] - g * (state["thermalConductance"][node] + state["thermalConductance"][node - 1])) * state["soilTemp"][node] + g * state["thermalConductance"][node] * state["soilTemp"][node + 1]
        a[surfaceNode] = 0.0
        sensibleHeatFlux = state["nu"] * state["thermalConductance"][airNode] * newTemps[airNode]
        if state["netRadiationSource"] == "calc":
            radnNet = Divide(state["netRadiation"] * 1000000.0, state["internalTimeStep"], 0.0)
        else:
            radnNet = Divide(waterBalance_Eos * 2465000.0, state["timestep"], 0.0)
        latentHeatFlux = Divide(waterBalance_Es * 2465000.0, state["timestep"], 0.0)
        soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
        d[surfaceNode] += soilSurfaceHeatFlux
        d[state["numNodes"]] += state["nu"] * state["thermalConductance"][state["numNodes"]] * newTemps[state["numNodes"] + 1]
        for node in range(surfaceNode, state["numNodes"]):
            c[node] = Divide(c[node], b[node], 0.0)
            d[node] = Divide(d[node], b[node], 0.0)
            b[node + 1] -= a[node + 1] * c[node]
            d[node + 1] -= a[node + 1] * d[node]
        newTemps[state["numNodes"]] = Divide(d[state["numNodes"]], b[state["numNodes"]], 0.0)
        for node in range(state["numNodes"] - 1, surfaceNode - 1, -1):
            newTemps[node] = d[node] - c[node] * newTemps[node + 1]

    def doUpdate(numInterationsPerDay):
        for i in range(len(state["newTemperature"])):
            state["soilTemp"][i] = state["newTemperature"][i]
        if state["timeOfDaySecs"] < state["internalTimeStep"] * 1.2:
            for node in range(surfaceNode, state["numNodes"] + 1):
                state["minSoilTemp"][node] = state["soilTemp"][node]
                state["maxSoilTemp"][node] = state["soilTemp"][node]
        for node in range(surfaceNode, state["numNodes"] + 1):
            if state["soilTemp"][node] < state["minSoilTemp"][node]:
                state["minSoilTemp"][node] = state["soilTemp"][node]
            elif state["soilTemp"][node] > state["maxSoilTemp"][node]:
                state["maxSoilTemp"][node] = state["soilTemp"][node]
            state["aveSoilTemp"][node] += Divide(state["soilTemp"][node], numInterationsPerDay, 0.0)
        state["boundaryLayerConductance"] += Divide(state["thermalConductivity"][airNode], numInterationsPerDay, 0.0)

    def calcSurfaceTemperature(waterBalance_Salb_, weather_MeanT_, weather_MaxT_, weather_Radn_):
        surfaceT = (1.0 - waterBalance_Salb_) * (weather_MeanT_ + (weather_MaxT_ - weather_MeanT_) * ((max(weather_Radn_, 0.1) * 23.8846 / 800.0) ** 0.5)) + waterBalance_Salb_ * weather_MeanT_
        return surfaceT

    import math

    doVolumetricSpecificHeat()
    doThermalConductivity()

    for timeStepIteration in range(1, interactionsPerDay + 1):
        state["timeOfDaySecs"] = state["internalTimeStep"] * float(timeStepIteration)
        if state["timestep"] < 24.0 * 60.0 * 60.0:
            state["airTemperature"] = weather_MeanT
        else:
            state["airTemperature"] = interpolateTemperature(state["timeOfDaySecs"] / 3600.0)
        state["newTemperature"][airNode] = state["airTemperature"]
        state["netRadiation"] = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva)

        if state["boundarLayerConductanceSource"] == "constant":
            state["thermalConductivity"][airNode] = 20.0
        else:
            state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"])
            for _ in range(1):
                doThomas(state["newTemperature"])
                state["thermalConductivity"][airNode] = getBoundaryLayerConductance(state["newTemperature"])
        doThomas(state["newTemperature"])
        doUpdate(interactionsPerDay)
        if abs(state["timeOfDaySecs"] - 5.0 * 3600.0) <= min(state["timeOfDaySecs"], 5.0 * 3600.0) * 0.0001:
            for i in range(len(state["soilTemp"])):
                state["morningSoilTemp"][i] = state["soilTemp"][i]

    state["minTempYesterday"] = weather_MinT
    state["maxTempYesterday"] = weather_MaxT

    # Prepare outputs
    numLayers = state["numLayers"]
    final_profile = [0.0] * numLayers
    for i in range(numLayers):
        final_profile[i] = state["soilTemp"][topsoilNode + i]
    average_profile = [0.0] * numLayers
    min_profile = [0.0] * numLayers
    max_profile = [0.0] * numLayers
    for i in range(numLayers):
        average_profile[i] = state["aveSoilTemp"][topsoilNode + i]
        min_profile[i] = state["minSoilTemp"][topsoilNode + i]
        max_profile[i] = state["maxSoilTemp"][topsoilNode + i]

    thermal_conductivity_profile = [0.0] * (state["numNodes"])
    for i in range(state["numNodes"]):
        thermal_conductivity_profile[i] = state["thermalConductivity"][i + 1]
    heat_capacity_profile = [0.0] * (state["numNodes"])
    for i in range(state["numNodes"]):
        heat_capacity_profile[i] = state["volSpecHeatSoil"][i + 1]
    heat_store_profile = [0.0] * (state["numNodes"])
    for i in range(state["numNodes"]):
        heat_store_profile[i] = state["heatStorage"][i + 1]
    thr_profile = [0.0] * (state["numNodes"] + 2)
    for i in range(len(thr_profile)):
        thr_profile[i] = state["morningSoilTemp"][i] if i < len(state["morningSoilTemp"]) else 0.0

    outputs = {
        "FinalSoilTemperature": final_profile,
        "FinalSoilSurfaceTemperature": state["soilTemp"][surfaceNode],
        "AverageSoilTemperature": average_profile,
        "AverageSoilSurfaceTemperature": state["aveSoilTemp"][surfaceNode],
        "MinimumSoilTemperature": min_profile,
        "MinimumSoilSurfaceTemperature": state["minSoilTemp"][surfaceNode],
        "MaximumSoilTemperature": max_profile,
        "MaximumSoilSurfaceTemperature": state["maxSoilTemp"][surfaceNode],
        "BoundaryLayerConductance": state["boundaryLayerConductance"],
        "ThermalConductivity": thermal_conductivity_profile,
        "HeatCapacity": heat_capacity_profile,
        "HeatStore": heat_store_profile,
        "Thr_profile": thr_profile
    }
    return state, outputs