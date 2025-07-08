from openai import OpenAI
import os
from pathlib import Path
import json
import xml.dom.minidom
import xml.etree.ElementTree as ET

# TO-DO
# - Documentation
# - Add error handling for
#   - API key extraction
#   - File reading
#   - Instructions extraction
#   - Prompt creation
#   - Response from GPT (process not found, etc...)
# - Integration into Crop2ML
#   - Get arguments from the command line (input folder + output_folder)
#   - Add a specific option to use LLM (-ai + name of model unit)
# - Composite model support (prompt + response from GPT)


#CONFIGURATION
INSTRUCTIONS_PATH = "./Instructions_crop2ML.txt"
API_KEY_PATH = "./api_key.txt"
MODEL = "gpt-4.1"

main_file = ""
folder_path = ""
output_path = ""
 

# Function to extract the file extension from a file path
# This function checks the file extension and returns it.
def extract_extension(file_path):
  fichier = Path(file_path)
  return fichier.suffix


# Function to determine the programming language based on file extension
# This function checks the file extension and returns the corresponding language name.
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
  else:
    return "Unknown"


# Function to connect to OpenAI's API
# This function reads the API key from a file and initializes the OpenAI client.
def extract_api_key(API_KEY_PATH):
  api_key = extract_text(API_KEY_PATH)
  try:
    OpenAI(api_key = api_key)
  except Exception as e:
    print(f"An error occurred while connecting to OpenAI: {e}")
    return None
  return api_key


# Function to extract text from a file
# This function reads the content of a file and returns it as a string.
def extract_text(file_path):
  with open(file_path, "r", encoding="utf-8") as file:
    return file.read()
  

# Function to send instructions and prompt to OpenAI's model
# This function takes instructions, a prompt, an API key, and a model name and returns the response from the model.
def send_to_gpt(instructions, prompt, api_key, model):
  client = OpenAI(api_key = api_key)

  response = client.chat.completions.create(
    model=model,
    messages=[
        {"role": "system", "content": instructions},
        {"role": "user", "content": prompt}
    ]
  )
  response = response.choices[0].message.content
  if response.startswith("```json"):
    response = response[7:].lstrip()
  if response.endswith("```"):
    response = response[:-3].rstrip()
  return(response)


# Function to list all files with a specific extension in a folder recursively except one specified file
# This function walks through the directory tree starting from the specified folder path and collects all files with the given extension.
def list_files_with_extension(folder_path, extension, exclude_filename):
  file_list = []
  for root, dirs, files in os.walk(folder_path):
    for file in files:
      if file.lower().endswith(extension) and file != exclude_filename:
        file_list.append(os.path.join(root, file))
  return file_list


# Function to create a prompt for the GPT model
# This function constructs a prompt based on the number of files in the specified folder, the main file, and the programming language.
def create_prompt(folder_path, main_file, list_files, language_name):
  nb_files = len(list_files)

  main_file_path = None
  for root, dirs, files in os.walk(folder_path):
    if main_file in files:
      main_file_path = os.path.join(root, main_file)
      break
  if main_file_path is None:
    raise FileNotFoundError(f"File {main_file} is not there.")
  
  prompt = f"You will be given {nb_files} crop model source code files, each written in {language_name}. Each file will be provided in a separate code block and will be labeled with its filename."
  prompt += f"\n\nThe main crop model source code file is called {main_file}, and its content is as follows: {extract_text(main_file_path)}"
  
  if nb_files > 0:
    number_file = 2
    prompt += "\n\nAdditional source files are as follows:"
    for file in list_files:
      prompt += f"\n\nSource code file {number_file} is {file}, and its content is as follows:\n"
      extract_text(file)
      number_file += 1
  
  prompt += "\nPlease follow the system instructions from the provided code."
  return prompt


# Function to convert JSON data to XML format
# This function takes a file path and JSON data, then converts the data into a Crop2ML-friendly XML format.
def json_to_xml(file_path, json_data):
  metadata = json_data['metadata']
  process = json_data['process']
  process_inputs = process['inputs']
  process_outputs = process['outputs']
  init = json_data['init']
  init_inputs = init['inputs']
  init_outputs = init['outputs']
  functions = json_data.get('functions', [])
  tests = json_data['tests']

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
  inputs = process_inputs
  outputs = process_outputs
  for init_input in init_inputs:
    if init_input.get('name') not in [inp.get('name') for inp in process_inputs]:
      inputs.append(init_input)
  for init_output in init_outputs:
    if init_output.get('name') not in [outp.get('name') for outp in process_outputs]:
      outputs.append(init_output)

  for function in functions:
    if function.get('inputs') :
      for inp in function.get('inputs'):
        if inp.get('name') not in [inp.get('name') for inp in process_inputs]:
          inputs.append(inp)
    if function.get('outputs'):
      for output in function.get('outputs'):
        if output.get('name') not in [outp.get('name') for outp in process_outputs]:
          outputs.append(output)

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

  if json_tests == [] :
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


# Main execution block
# This block extracts the OpenAI API key, the programming language, the system instructions, list the relevant files in the specified folder, and constructs a prompt for the GPT model.
# It then sends the prompt to the model and saves the response as a JSON file and a XML/Crop2ML file. 
if __name__ == "__main__":
  api_key = extract_api_key(API_KEY_PATH)
  extension = extract_extension(folder_path + main_file)
  language_name = language(extension)
  list_files = list_files_with_extension(folder_path, extension, main_file)

  instructions = extract_text(INSTRUCTIONS_PATH)
  prompt = create_prompt(folder_path, main_file, list_files, language_name)
  response = send_to_gpt(instructions, prompt, api_key, MODEL)

  base, _ = os.path.splitext(main_file)
  json_path = output_path + base + ".json"
  json_obj = json.loads(response)

  with open(json_path, "w", encoding="utf-8") as f:
    json.dump(json_obj, f, ensure_ascii=False, indent=4)

  xml_path = output_path + "unit." + base + ".xml"
  xml_data = json_to_xml(folder_path, json_obj)
  dom = xml.dom.minidom.parseString(xml_data)
  with open(xml_path, 'w', encoding='utf-8') as f:
    f.write(dom.toprettyxml())
