from openai import OpenAI
import os
from pathlib import Path
import json
import xml.dom.minidom
import xml.etree.ElementTree as ET
from time import sleep


# TO-DO
#-----------------------------------------------------------------
#1 - Integration into Crop2ML
#1.1   - Add a specific option to use LLM (-ai + name of model unit)
#1.2   - Get arguments from the command line (input folder + output_folder)
#2 - Composite model support (prompt + response from GPT)
#3 - Add error handling for
#3.1   - API key extraction
#3.2   - Instructions extraction
#3.3   - File reading
#3.4   - Prompt creation
#3.5   - Response from GPT (process not found, etc...)
#4 - Documentation
#-----------------------------------------------------------------


#-----------------------------------------------------------------
# Function to extract the file extension from a file path
#-----------------------------------------------------------------
def extract_extension(file_path):
  fichier = Path(file_path)
  return fichier.suffix


#-----------------------------------------------------------------
# Function to determine the programming language based on file extension
# This function checks the file extension and returns the corresponding language name.
#-----------------------------------------------------------------
def language(extension):
  if extension == ".py":
    return "Python"
  elif extension == ".java":
    return "Java"
  elif extension == ".cs":
    return "C#"
  elif extension == ".cpp":
    return "C++"
  elif extension == ".for":
    return "Fortran"
  elif extension == ".f90":
    return "Fortran90"
  elif extension == ".pyx":
    return "Cython"
  else:
    return "Unknown"


#-----------------------------------------------------------------
# Function to connect to OpenAI's API
# This function reads the API key from a file and initializes the OpenAI client.
#-----------------------------------------------------------------
def extract_api_key(API_KEY_PATH):
  api_key = extract_text(API_KEY_PATH)
  try:
    OpenAI(api_key = api_key)
  except Exception as e:
    print(f"An error occurred while connecting to OpenAI: {e}")
    return None
  return api_key


#-----------------------------------------------------------------
# Function to extract text from a file
# This function reads the content of a file and returns it as a string.
#-----------------------------------------------------------------
def extract_text(file_path):
  with open(file_path, "r", encoding="utf-8", errors="replace") as file:
    return file.read()
  

#-----------------------------------------------------------------
# Function to list all files with a specific extension in a folder recursively except one specified file
# This function walks through the directory tree starting from the specified folder path and collects all files with the given extension.
#-----------------------------------------------------------------
def list_files_with_extension(folder_path, extension, exclude_filename):
  file_list = []
  for root, dirs, files in os.walk(folder_path):
    for file in files:
      if extension == ".cpp" and file.lower().endswith(".h") and file != exclude_filename:
        file_list.append(os.path.join(root, file))
      elif file.lower().endswith(extension) and file != exclude_filename:
        file_list.append(os.path.join(root, file))
  return file_list


#-----------------------------------------------------------------
# Function to create a prompt for the GPT model
# This function constructs a prompt based on the number of files in the specified folder, the main file, and the programming language.
#-----------------------------------------------------------------
def create_prompt_all_files(folder_path, main_file, list_files, language_name):
  nb_files = len(list_files)

  main_file_path = None
  for root, dirs, files in os.walk(folder_path):
    if main_file in files:
      main_file_path = os.path.join(root, main_file)
      break
  if main_file_path is None:
    raise FileNotFoundError(f"File {main_file} is not there.")
  
  prompt = f"Analyze the {nb_files + 1} following crop model source code files, written in {language_name} as a single crop model component. The main file is {main_file}.\n"
  prompt += f"Follow the system instructions. Each file is marked clearly with --- FILE: filename --- at the start and --- END FILE --- at the end.\n\n\n"
  prompt += f"--- FILE: {main_file} ---\n{extract_text(main_file_path)}\n--- END FILE ---"
  
  if nb_files > 0:
    for file in list_files:
      prompt += f"\n\n--- FILE: {os.path.basename(file)} ---\n{extract_text(file)}\n--- END FILE ---"
  
  return prompt


#-----------------------------------------------------------------
# Function to create a prompt for the GPT model
# This function constructs a prompt based on the refactored code.
#-----------------------------------------------------------------
def create_prompt_refactor(code_refactored):
  prompt = f"Analyze the following Python crop model source code as a single crop model component.\n"
  prompt += f"Follow the system instructions and output a JSON file describing the model.\n\n"
  prompt += f"{code_refactored}"
  return prompt


#-----------------------------------------------------------------
# Function to send instructions and prompt to OpenAI's model
# This function takes instructions, a prompt, an API key, and a model name and returns the response from the model.
#-----------------------------------------------------------------
def send_to_gpt(instructions, prompt, api_key, model, reasoning_effort, text_format, verbosity):
  client = OpenAI(api_key = api_key)

  response = client.responses.create(
    model=model,
    reasoning={"effort": reasoning_effort},
    store=True,
    text={
      "format": {"type": text_format},
      "verbosity": verbosity,
      },
    input=[
      {"role": "developer",
        "content": [{"type": "input_text", "text": instructions}],
      },
      {
        "role": "user",
        "content": [{"type": "input_text", "text": prompt}],
      }
    ],
  )

  response = response.output_text
  if response.startswith("```json"):
    response = response[7:].lstrip()
  if response.endswith("```"):
    response = response[:-3].rstrip()
  return(response)


#-----------------------------------------------------------------
# Function to convert JSON data to XML format
# This function takes a file path and JSON data, then converts the data into a Crop2ML-friendly XML format.
#-----------------------------------------------------------------
def json_to_xml(file_path, json_metadata, json_code):
  metadata = json_metadata['metadata']
  init = json_code['init']
  process = json_code['process']
  inputs = json_code.get('inputs',[])
  outputs = json_code.get('outputs',[])
  functions = json_code.get('functions', [])
  tests = json_code['tests']

  # Create XML tree
  root = ET.Element('ModelUnit', {
    "modelid": os.path.basename(file_path) + "." + metadata['Title'],
    "name": metadata['Title'],
    "timestep": "1",
    "version":  metadata['Model version']
  })

  # Description section
  desc = ET.SubElement(root, 'Description')
  ET.SubElement(desc, 'Title').text = metadata.get('Title','')
  ET.SubElement(desc, 'Authors').text = metadata.get('Authors', '')
  ET.SubElement(desc, 'Institution').text = metadata.get('Institution', '')
  ET.SubElement(desc, 'URI').text = metadata.get('URI', '')
  ET.SubElement(desc, 'Reference').text = metadata.get('DOI', '')
  ET.SubElement(desc, 'ExtendedDescription').text = metadata.get('Extended description', '')
  ET.SubElement(desc, 'ShortDescription').text = metadata.get('Short description', '')

  # I/O section
  xml_inputs = ET.SubElement(root, 'Inputs')
  add_inputs(xml_inputs, inputs)

  xml_outputs = ET.SubElement(root, 'Outputs')
  add_outputs(xml_outputs, outputs)

  # Initialization
  if init['name'] != '-':
    ET.SubElement(root, 'Initialization', {
      'name': init['name'],
      'language': 'cyml',
      'filename': f"algo/pyx/{init['name']}.pyx"
    })

  # Functions
  for func in functions:
    ET.SubElement(root, 'Function', {
      'name': func['name'],
      'description': func['description'],
      'language': 'cyml',
      'type': 'external',
      'filename': f"algo/pyx/{func['name']}.pyx"
    })

  # Main Algorithm
  ET.SubElement(root, 'Algorithm', {
    'language': 'cyml',
    'platform': '',
    'filename': 'algo/pyx/' + process['name'] + ".pyx"
  })

  # Parametersets
  add_tests(root, tests, inputs)

  return ET.tostring(root, encoding='utf-8')


# Function: convert JSON 'inputs' to Crop2ML-friendly XML inputs
def add_inputs(xml_inputs, json_inputs):
  for input in json_inputs:
    attrs = {
        'name': str(input['name']),
        'description': str(input.get('description', '')),
        'inputtype': str(input.get('inputtype', ''))
    }
    if input.get('inputtype') == 'parameter':
      attrs['parametercategory'] = str(input.get('category', ''))
    else:
      attrs['variablecategory'] = str(input.get('category', ''))
    attrs['datatype'] = str(input.get('datatype', ''))
    if input.get('datatype') == "DOUBLEARRAY" or input.get('datatype') == "DOUBLELIST":
      attrs['len'] = str(input.get('len', ''))
    attrs['max'] = str(input.get('max', ''))
    attrs['min'] = str(input.get('min', ''))
    attrs['default'] = str(input.get('default', ''))
    attrs['unit'] = str(input.get('unit', ''))
    attrs['uri'] = str(input.get('uri', ''))

    ET.SubElement(xml_inputs, 'Input', attrs)


# Function: convert JSON 'outputs' to Crop2ML-friendly XML outputs
def add_outputs(xml_outputs, json_outputs):
  for output in json_outputs:
    attrs = {
      'name': str(output['name']),
      'description': str(output.get('description', '')),
      'variablecategory': str(output.get('category', '')),
      'datatype': str(output.get('datatype', ''))
    }
    if output.get('datatype') == 'DOUBLEARRAY' or output.get('datatype') == 'DOUBLELIST':
      attrs['len'] = str(output.get('len', ''))
    attrs['max'] = str(output.get('max', ''))
    attrs['min'] = str(output.get('min', ''))
    attrs['unit'] = str(output.get('unit', ''))
    attrs['uri'] = str(output.get('uri', ''))
    ET.SubElement(xml_outputs, 'Output', attrs)


# Function: convert JSON 'tests' to Crop2ML-friendly XML test
def add_tests(root_XML, json_tests, json_inputs):
  parametersSets = ET.SubElement(root_XML, 'Parametersets')
  testSets = ET.SubElement(root_XML, 'Testsets')

  if json_tests == [] or json_tests[0] == "-":
    return
  
  parameter_inputs = []
  variable_inputs = []
  inputtype_by_name = {inp['name']: inp.get('inputtype', '') for inp in json_inputs}

  for test in json_tests:
    test_inputs = test['inputs']
    test_outputs = test['outputs']

    for test_input in test_inputs:
      name = test_input['name']
      if inputtype_by_name.get(name) is not None:
        if inputtype_by_name.get(name) == 'parameter':
          parameter_inputs.append(test_input)
        else:
          variable_inputs.append(test_input)

    if len(parameter_inputs) > 0:
      parameterset = ET.SubElement(parametersSets, 'Parameterset', {
        'name': "p_" + test.get('name'),
        'description': test.get('description', '')
      })
      for parameter_input in parameter_inputs:
        ET.SubElement(parameterset, 'Param', name=parameter_input.get('name')).text = str(parameter_input.get('value'))

    if len(variable_inputs) > 0 or len(test_outputs) > 0:
      testset = ET.SubElement(testSets, 'Testset', {
        'name': "t_" + test.get('name'),
        'description': test.get('description', '')
      })
      if len(parameter_inputs) > 0:
        testset.set('parameterset', "p_" + test.get('name'))
      test = ET.SubElement(testset, 'Test', name=testset.get('name'))
      for variable_input in variable_inputs:
        ET.SubElement(test, 'InputValue', name=variable_input.get('name')).text = str(variable_input.get('value'))
      for test_output in test_outputs:
        ET.SubElement(test, 'OutputValue', name=test_output.get('name')).text = str(test_output.get('value'))


#-----------------------------------------------------------------
# Main execution block
# This block extracts the OpenAI API key, the programming language, the system instructions, 
# list the relevant files in the specified folder, and constructs a prompt for the GPT model.
# It then sends the prompt to the model and saves the response as a JSON file and a XML/Crop2ML file. 
#-----------------------------------------------------------------
def main(api_key_path, agent1_path, agent2_path, agent3_path, model, folder_path, output_path, main_file):
  api_key = extract_api_key(api_key_path)
  extension = extract_extension(folder_path + main_file)
  language_name = language(extension)
  list_files = list_files_with_extension(folder_path, extension, main_file)

  instructions_metadata = extract_text(agent1_path)
  instructions_refactor = extract_text(agent2_path)
  instructions_json = extract_text(agent3_path)

  prompt_all_files = create_prompt_all_files(folder_path, main_file, list_files, language_name)
  response_metadata = send_to_gpt(instructions_metadata, prompt_all_files, api_key, model, "high", "json_object", "low")
  response_refactored = send_to_gpt(instructions_refactor, prompt_all_files, api_key, model, "high", "text", "medium")

  prompt_refactored = create_prompt_refactor(response_refactored)
  response_json = send_to_gpt(instructions_json, prompt_refactored, api_key, model, "high", "json_object", "low")

  os.makedirs(output_path, exist_ok=True)
  base, _ = os.path.splitext(main_file)
  json_metadata_path = output_path + base + "_metadata.json"
  json_code_path = output_path + base + "_code.json"
  python_code_path = output_path + base + "_code.py"
  json_metadata = json.loads(response_metadata)
  json_code = json.loads(response_json)

  with open(json_metadata_path, "w", encoding="utf-8") as f:
    json.dump(json_metadata, f, ensure_ascii=False, indent=4)

  with open(json_code_path, "w", encoding="utf-8") as f:
    json.dump(json_code, f, ensure_ascii=False, indent=4)

  with open(python_code_path, "w", encoding="utf-8") as f:
    f.write(response_refactored)

  xml_path = output_path + "unit." + base + ".xml"
  xml_data = json_to_xml(folder_path, json_metadata, json_code)
  dom = xml.dom.minidom.parseString(xml_data)
  with open(xml_path, 'w', encoding='utf-8') as f:
    f.write(dom.toprettyxml())


#-----------------------------------------------------------------
#CONFIGURATION
#-----------------------------------------------------------------
AGENT1_PATH = "./1-agent_metadata.txt"
AGENT2_PATH = "./2-agent_refactor.txt"
AGENT3_PATH = "./3-agent_JSON.txt"
API_KEY_PATH = "./api_key.txt"
MODEL = "gpt-5"

COMPONENTS_DICT = {
  "../Components/ApsimCampbell/": "SoilTemperature.cs",
  "../Components/BiomaSurfacePartonSoilSWATC/": "SoilTemperatureSWAT.cs",
  "../Components/BiomaSurfaceSWATSoilSWATC/": "SoilTemperatureSWAT.cs",
  "../Components/DSSAT_EPICST_standalone/": "STEMP_EPIC.for",
  "../Components/DSSAT_ST_standalone/": "STEMP.for",
  "../Components/Simplace_Soil_Temperature/": "STMPsimCalculator.java",
  "../Components/SQ_Soil_Temperature/": "CalculateSoilTemperature.cs",
  "../Components/Stics_soil_temperature/": "Tempprofile.f90"
}

ENERGY_BALANCE_DIRECTORY = "../Components/SQ_EnergyBalance/"
ENERGY_BALANCE_FILES =  [
    "Evapotranspiration.cs",
    "Netradiation.cs",
    "Netradiationequivalentevaporation.cs",
    "Penman.cs",
    "Potentialtranspiration.cs",
    "Priestlytaylor.cs",
    "Ptsoil.cs",
    "Soilevaporation.cs",
    "Soilheatflux.cs",
    "Canopytemperature.cs",
    "Cropheatflux.cs",
    "Diffusionlimitedevaporation.cs"
  ]

first_index_component = 0
last_index_component = len(COMPONENTS_DICT) - 1
number_iteration = 10
component_keys = list(COMPONENTS_DICT.keys())
component_values = list(COMPONENTS_DICT.values())

#-----------------------------------------------------------------
# Simulation section
#-----------------------------------------------------------------
# Soil temperature section
'''for i in range(first_index_component, last_index_component + 1):
  print(f"Processing {component_keys[i]} : {component_values[i]}")
  for j in range(5, number_iteration + 1): 
    print(f"Iteration {j}")
    output_path = f"{component_keys[i]}{j}/"
    main(API_KEY_PATH, AGENT1_PATH, AGENT2_PATH, AGENT3_PATH, MODEL, component_keys[i], output_path, component_values[i])
'''
# Energy balance section
for i in range(11, len(ENERGY_BALANCE_FILES)):
  print(f"Processing {ENERGY_BALANCE_DIRECTORY} : {ENERGY_BALANCE_FILES[i]}")
  for j in range(1, number_iteration + 1): 
    print(f"Iteration {j}")
    output_path = f"{ENERGY_BALANCE_DIRECTORY}{j}/"
    main(API_KEY_PATH, AGENT1_PATH, AGENT2_PATH, AGENT3_PATH, MODEL, ENERGY_BALANCE_DIRECTORY, output_path, ENERGY_BALANCE_FILES[i])