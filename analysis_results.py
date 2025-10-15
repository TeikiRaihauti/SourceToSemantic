import os
import xml.etree.ElementTree as ET
import csv
import pandas as pd
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

# -------------------------------------------------
# Function for extracting file extension
# -------------------------------------------------
def extract_extension(file_path):
  fichier = Path(file_path)
  return fichier.suffix


# -------------------------------------------------
# Function for extracting text from a file
# -------------------------------------------------
def extract_text(file_path):
  with open(file_path, "r", encoding="utf-8", errors="replace") as file:
    return file.read()


# -------------------------------------------------
# Function for listing files with a specific extension in a folder, excluding a specific filename
# -------------------------------------------------
def list_files_with_extension(folder_path, extension, exclude_filename):
  file_list = []
  for root, dirs, files in os.walk(folder_path):
    for file in files:
      if extension == ".cpp" and file.lower().endswith(".h") and file != exclude_filename:
        file_list.append(os.path.join(root, file))
      elif file.lower().endswith(extension) and file != exclude_filename:
        file_list.append(os.path.join(root, file))
  return file_list


# -------------------------------------------------
# Function for counting code and comment lines in a file based on comment rules
# -------------------------------------------------
def count_code_and_comment_lines(file_path, comment_rules):
    code_lines = 0
    comment_lines = 0
    total_lines = 0
    in_block_comment = False
    block_start = block_end = None

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            total_lines += 1
            stripped = line.strip()
            if not stripped:
                continue 

            # Check for block comments
            if not in_block_comment and comment_rules["block"]:
                for start, end in comment_rules["block"]:
                    if stripped.startswith(start):
                        in_block_comment = True
                        block_start, block_end = start, end
                        comment_lines += 1
                        break
            elif in_block_comment:
                comment_lines += 1
                if block_end and block_end in stripped:
                    in_block_comment = False
                continue

            # Check for line comments
            if not in_block_comment and any(stripped.startswith(sym) for sym in comment_rules["line"]):
                comment_lines += 1
            elif not in_block_comment:
                code_lines += 1

    return total_lines, code_lines, comment_lines


# -------------------------------------------------
# Function for finding the main file in a folder recursively
# -------------------------------------------------
def find_main_file(folder_path, main_file):
    for root, dirs, files in os.walk(folder_path):
        if main_file in files:
            return os.path.join(root, main_file)
    raise FileNotFoundError(f"{main_file} not found in {folder_path}")


#-----------------------------------------------------------------
# Parse the XML model file and extract inputs, outputs, and functions
#-----------------------------------------------------------------
def parse_model_file(xml_file):
  tree = ET.parse(xml_file)
  root = tree.getroot()
  
  inputs = {inp.attrib['name']: inp.attrib for inp in root.find('Inputs')}
  outputs = {out.attrib['name']: out.attrib for out in root.find('Outputs')}
  functions = {func.attrib['name']: func.attrib for func in root.findall('Function')}
  
  return inputs, outputs, functions


#-----------------------------------------------------------------
# Compare the reference and test XML elements and count the true positive, false positive, and false negative elements
#-----------------------------------------------------------------
def compare_elements(ref_dict, test_dict):
  tp = {}
  fp = {}
  fn = {}

  for key in test_dict:
    if key in ref_dict:
      tp[key] = (ref_dict[key], test_dict[key])
    else:
      fp[key] = test_dict[key]
  
  for key in ref_dict:
    if key not in test_dict:
      fn[key] = ref_dict[key]
  
  return tp, fp, fn

#-----------------------------------------------------------------
# Calculate precision, recall, and F1 score based on true positives, false positives, and false negatives
#-----------------------------------------------------------------
def calculate_metrics(tp, fp, fn):
  precision = len(tp) / (len(tp) + len(fp)) if (len(tp) + len(fp)) > 0 else 0
  recall = len(tp) / (len(tp) + len(fn)) if (len(tp) + len(fn)) > 0 else 0
  f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
  return precision, recall, f1_score


#-----------------------------------------------------------------
# Compare fields of reference and test elements and returns a list with the comparison results
#-----------------------------------------------------------------
def compare_fields(ref_element, test_element):
  comparison = [ref_element.get('name', '')]
  fields = ["inputtype", "variablecategory", "datatype", "len", "max", "min", "default", "unit", "uri"]

  for field in fields:
    ref_value = ref_element.get(field, '')
    test_value = test_element.get(field, '')
    if ref_value == test_value:
      comparison.append("OK")
    else:
      comparison.append("WRONG")
  
  return comparison


#-----------------------------------------------------------------
# Compare fields of reference and test elements and returns a list with the comparison results
#-----------------------------------------------------------------
def compare_fields_aggregated(ref_inputs, test_inputs, tp_in):
  nb_parameter = nb_parameter_ok = nb_state = nb_state_ok = nb_rate = nb_rate_ok = nb_exogenous = nb_exogenous_ok = \
    nb_auxiliary = nb_auxiliary_ok = nb_int = nb_int_ok = nb_double = nb_double_ok = \
    nb_doublearray = nb_doublearray_ok = nb_doublelist = nb_doublelist_ok = \
    nb_char = nb_char_ok = nb_date = nb_date_ok = nb_intarray = nb_intarray_ok = \
    nb_len = nb_len_ok = 0

  for input in tp_in:
    test_element = test_inputs[input]
    ref_element = ref_inputs[input]

    if ref_element.get("inputtype") == "parameter":
      nb_parameter += 1
      test_value = test_element.get("parametercategory", '')
      ref_value = ref_element.get("parametercategory", '')
      if test_value == ref_value:
        nb_parameter_ok += 1

    elif ref_element.get("inputtype") == "variable":
      test_value = test_element.get("variablecategory", '')
      ref_value = ref_element.get("variablecategory", '')

      if ref_value == "state":
        nb_state += 1
        if test_value == ref_value:
          nb_state_ok += 1
      elif ref_value == "rate":
        nb_rate += 1
        if ref_value == ref_value:
          nb_rate_ok += 1
      elif ref_value == "exogenous":
        nb_exogenous += 1
        if test_value == ref_value:
          nb_exogenous_ok += 1
      elif ref_value == "auxiliary":
        nb_auxiliary += 1
        if test_value == ref_value:
          nb_auxiliary_ok += 1

    test_value = test_element.get("datatype", '')
    ref_value = ref_element.get("datatype", '')

    if ref_value == "INT":
      nb_int += 1
      if test_value == ref_value:
        nb_int_ok += 1
    elif ref_value == "DOUBLE":
      nb_double += 1
      if test_value == ref_value:
        nb_double_ok += 1
    elif ref_value == "DOUBLEARRAY":
      nb_doublearray += 1
      nb_len += 1
      if test_value == ref_value:
        nb_doublearray_ok += 1
      if test_element.get("len", '') == ref_element.get("len", ''):
        nb_len_ok += 1
    elif ref_value == "DOUBLELIST":
      nb_doublelist += 1
      nb_len += 1
      if test_value == ref_value:
        nb_doublelist_ok += 1
      if test_element.get("len", '') == ref_element.get("len", ''):
        nb_len_ok += 1
    elif ref_value == "CHAR":
      nb_char += 1
      if test_value == ref_value:
        nb_char_ok += 1
    elif ref_value == "DATE":
      nb_date += 1
      if test_value == ref_value:
        nb_date_ok += 1
    elif ref_value == "INTARRAY":
      nb_intarray += 1
      nb_len += 1
      if test_value == ref_value:
        nb_intarray_ok += 1
      if test_element.get("len", '') == ref_element.get("len", ''):
        nb_len_ok += 1

  return (nb_parameter, nb_parameter_ok, nb_state, nb_state_ok, nb_rate, nb_rate_ok, nb_exogenous, nb_exogenous_ok,
  nb_auxiliary, nb_auxiliary_ok, nb_int, nb_int_ok, nb_double, nb_double_ok,
  nb_doublearray, nb_doublearray_ok, nb_doublelist, nb_doublelist_ok,
  nb_char, nb_char_ok, nb_date, nb_date_ok, nb_intarray, nb_intarray_ok, nb_len, nb_len_ok)


#-----------------------------------------------------------------
# Write the comparative analysis results between two files into a CSV file
# Calculate precision, recall, and F1 score for inputs, outputs, and functions then write the results for each section
#-----------------------------------------------------------------
def write_comparison_csv(output_file, ref_inputs, test_inputs, ref_outputs, test_outputs, ref_func, test_func, 
                         tp_in, fp_in, fn_in, tp_out, fp_out, fn_out, tp_func, fp_func, fn_func):
  
  ### Store the results in a list
  precision_in, recall_in, f1_in = calculate_metrics(tp_in, fp_in, fn_in)
  precision_out, recall_out, f1_out = calculate_metrics(tp_out, fp_out, fn_out)
  precision_func, recall_func, f1_func = calculate_metrics(tp_func, fp_func, fn_func)

  with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter=";")                                             

    ### Inputs section
    writer.writerow(["Inputs stats"])
    writer.writerow(["Crop2ML-generated","AI-generated","True positive","False positive","False negative", "Precision", "Recall", "F1-Score"])
    writer.writerow([len(ref_inputs), len(test_inputs), len(tp_in), len(fp_in), len(fn_in), precision_in, recall_in, f1_in])
    writer.writerow([])

    writer.writerow(["True positive"])
    writer.writerow(["Name","inputtype","variablecategory","datatype","len","max","min","default","unit","uri"])
    for input in tp_in:
      writer.writerow(compare_fields(ref_inputs[input], test_inputs[input]))
    writer.writerow([])

    writer.writerow(["False positive"])
    writer.writerow(["Name","inputtype","variablecategory","datatype"])
    for input in fp_in:
      writer.writerow([fp_in[input].get('name',''), fp_in[input].get('inputtype',''), fp_in[input].get('variablecategory',''), fp_in[input].get('datatype','')])
    writer.writerow([])

    writer.writerow(["False negative"])
    writer.writerow(["Name","inputtype","variablecategory","datatype"])
    for input in fn_in:
      writer.writerow([fn_in[input].get('name',''), fn_in[input].get('inputtype',''), fn_in[input].get('variablecategory',''), fn_in[input].get('datatype','')])
    writer.writerow([])
    writer.writerow([])

    ### Outputs section
    writer.writerow(["Outputs stats"])
    writer.writerow(["Crop2ML-generated","AI-generated","True positive","False positive","False negative", "Precision", "Recall", "F1-Score"])
    writer.writerow([len(ref_outputs), len(test_outputs), len(tp_out), len(fp_out), len(fn_out), precision_out, recall_out, f1_out])
    writer.writerow([])

    writer.writerow(["True positive"])
    writer.writerow(["Name","inputtype","variablecategory","datatype","len","max","min","default","unit","uri"])
    for output in tp_out:
      writer.writerow(compare_fields(ref_outputs[output], test_outputs[output]))
    writer.writerow([])

    writer.writerow(["False positive"])
    writer.writerow(["Name","inputtype","variablecategory","datatype"])
    for output in fp_out:
      writer.writerow([fp_out[output].get('name',''), fp_out[output].get('inputtype',''), fp_out[output].get('variablecategory',''), fp_out[output].get('datatype','')])
    writer.writerow([])

    writer.writerow(["False negative"])
    writer.writerow(["Name","inputtype","variablecategory","datatype"])
    for output in fn_out:
      writer.writerow([fn_out[output].get('name',''), fn_out[output].get('inputtype',''), fn_out[output].get('variablecategory',''), fn_out[output].get('datatype','')])
    writer.writerow([])
    writer.writerow([])

    ### Function section
    writer.writerow(["Functions stats"])
    writer.writerow(["Crop2ML-generated","AI-generated","True positive","False positive","False negative", "Precision", "Recall", "F1-Score"])
    writer.writerow([len(ref_func), len(test_func), len(tp_func), len(fp_func), len(fn_func), precision_func, recall_func, f1_func])
    writer.writerow([])

    writer.writerow(["True positive"])
    writer.writerow(["Name"])
    for func in tp_func:
      writer.writerow([tp_func[func][1].get('name','')])
    writer.writerow([])

    writer.writerow(["False positive"])
    writer.writerow(["Name"])
    for func in fp_func:
      writer.writerow([fp_func[func].get('name','')])
    writer.writerow([])

    writer.writerow(["False negative"])
    writer.writerow(["Name"])
    for func in fn_func:
      writer.writerow([fn_func[func].get('name','')])
    writer.writerow([])


#-----------------------------------------------------------------
# Return the comparative analysis results between two files
#-----------------------------------------------------------------
def result_concatenated(ref_inputs, test_inputs, ref_outputs, test_outputs, ref_func, test_func, 
                         tp_in, fp_in, fn_in, tp_out, fp_out, fn_out, tp_func, fp_func, fn_func):
  result = []
  precision_in, recall_in, f1_in = calculate_metrics(tp_in, fp_in, fn_in)
  precision_out, recall_out, f1_out = calculate_metrics(tp_out, fp_out, fn_out)
  precision_func, recall_func, f1_func = calculate_metrics(tp_func, fp_func, fn_func)  

  nb_parameter, nb_parameter_ok, nb_state, nb_state_ok, nb_rate, nb_rate_ok, nb_exogenous, nb_exogenous_ok, \
  nb_auxiliary, nb_auxiliary_ok, nb_int, nb_int_ok, nb_double, nb_double_ok, \
  nb_doublearray, nb_doublearray_ok, nb_doublelist, nb_doublelist_ok, \
  nb_char, nb_char_ok, nb_date, nb_date_ok, nb_intarray, nb_intarray_ok, nb_len, nb_len_ok = \
    compare_fields_aggregated(ref_inputs, test_inputs, tp_in)

  result += [len(ref_inputs), len(test_inputs), len(tp_in), len(fp_in), len(fn_in), precision_in, recall_in, f1_in]
  result += [nb_parameter, nb_parameter_ok, nb_state, nb_state_ok, nb_rate, nb_rate_ok, nb_exogenous, nb_exogenous_ok,
            nb_auxiliary, nb_auxiliary_ok, nb_int, nb_int_ok, nb_double, nb_double_ok,
            nb_doublearray, nb_doublearray_ok, nb_doublelist, nb_doublelist_ok,
            nb_char, nb_char_ok, nb_date, nb_date_ok, nb_intarray, nb_intarray_ok, 
            nb_len, nb_len_ok]
            
  nb_parameter, nb_parameter_ok, nb_state, nb_state_ok, nb_rate, nb_rate_ok, nb_exogenous, nb_exogenous_ok, \
  nb_auxiliary, nb_auxiliary_ok, nb_int, nb_int_ok, nb_double, nb_double_ok, \
  nb_doublearray, nb_doublearray_ok, nb_doublelist, nb_doublelist_ok, \
  nb_char, nb_char_ok, nb_date, nb_date_ok, nb_intarray, nb_intarray_ok, nb_len, nb_len_ok = \
    compare_fields_aggregated(ref_outputs, test_outputs, tp_out)

  result += [len(ref_outputs), len(test_outputs), len(tp_out), len(fp_out), len(fn_out), precision_out, recall_out, f1_out]
  result += [nb_int, nb_int_ok, nb_double, nb_double_ok,
            nb_doublearray, nb_doublearray_ok, nb_doublelist, nb_doublelist_ok,
            nb_char, nb_char_ok, nb_date, nb_date_ok, nb_intarray, nb_intarray_ok, 
            nb_len, nb_len_ok]
  result += [len(ref_func), len(test_func), len(tp_func), len(fp_func), len(fn_func), precision_func, recall_func, f1_func]

  return result


#-----------------------------------------------------------------
# Filter keys after the first underscore in a dictionary
# This is used to normalize the keys for comparison
#-----------------------------------------------------------------
def filter_keys_after_underscore(d):
  return { (k.split("_", 1)[1] if "_" in k else k).lower(): v for k, v in d.items() }

#-----------------------------------------------------------------
# Main call to parse the model files, compare elements, and write the results to a CSV file
#-----------------------------------------------------------------
def main(REF_FILE, TEST_FILE, output_csv):
  ref_inputs, ref_outputs, ref_func = parse_model_file(REF_FILE)
  test_inputs, test_outputs, test_func = parse_model_file(TEST_FILE)

  ref_inputs = filter_keys_after_underscore(ref_inputs)
  test_inputs = filter_keys_after_underscore(test_inputs)
  ref_outputs = filter_keys_after_underscore(ref_outputs)
  test_outputs = filter_keys_after_underscore(test_outputs)
  ref_func = filter_keys_after_underscore(ref_func)
  test_func = filter_keys_after_underscore(test_func)

  tp_in, fp_in, fn_in = compare_elements(ref_inputs, test_inputs)
  tp_out, fp_out, fn_out = compare_elements(ref_outputs, test_outputs)
  tp_func, fp_func, fn_func = compare_elements(ref_func, test_func)

  write_comparison_csv(output_csv, ref_inputs, test_inputs, ref_outputs, test_outputs, ref_func, test_func, 
                         tp_in, fp_in, fn_in, tp_out, fp_out, fn_out, tp_func, fp_func, fn_func)
  
  return result_concatenated(ref_inputs, test_inputs, ref_outputs, test_outputs, ref_func, test_func, 
                         tp_in, fp_in, fn_in, tp_out, fp_out, fn_out, tp_func, fp_func, fn_func)

#-----------------------------------------------------------------
# Function to count_line
#-----------------------------------------------------------------
def count_line(folder, main_file, iteration):
  PYTHON_DICT = {
    "../Components/ApsimCampbell/": "SoilTemperature_code.py",
    "../Components/BiomaSurfacePartonSoilSWATC/": "SoilTemperatureSWAT_code.py",
    "../Components/BiomaSurfaceSWATSoilSWATC/": "SoilTemperatureSWAT_code.py",
    "../Components/DSSAT_EPICST_standalone/": "STEMP_EPIC_code.py",
    "../Components/DSSAT_ST_standalone/": "STEMP_code.py",
    "../Components/Simplace_Soil_Temperature/": "STMPsimCalculator_code.py",
    "../Components/SQ_Soil_Temperature/": "CalculateSoilTemperature_code.py",
    "../Components/Stics_soil_temperature/": "Tempprofile_code.py"
  }
  
  result = []
  extension = extract_extension(main_file)
  list_files = list_files_with_extension(folder, extension, main_file)
  try:
    main_file_path = find_main_file(folder, main_file)
  except FileNotFoundError:
    return 0, 0, 0, 0, 0, 0
  comment_rules = RULES_BY_EXT.get(extension[1:], RULES_BY_EXT["py"])
  main_total, main_code, main_comments = count_code_and_comment_lines(main_file_path, comment_rules)
  python_total, python_code, python_comments = count_code_and_comment_lines(f"{folder}/{iteration}/{PYTHON_DICT[folder]}", RULES_BY_EXT["py"])
  other_total = other_code = other_comments = 0
  for file in list_files:
    t, c, com = count_code_and_comment_lines(file, comment_rules)
    other_total += t
    other_code += c
    other_comments += com
  total_comments = main_comments + other_comments
  result += [main_code, other_code, total_comments, python_total]
  return result


#-----------------------------------------------------------------
#CONFIGURATION
#-----------------------------------------------------------------
RULES_BY_EXT = {
    "cs":  {"line": ["//"], "block": [("/*", "*/")]},
    "java": {"line": ["//"], "block": [("/*", "*/")]},
    "py":  {"line": ["#"], "block": [("'''", "'''"), ('"""', '"""')]},
    "for": {"line": ["!"], "block": []},
    "f90": {"line": ["!"], "block": []}
}

BENCHMARK_FILES =  {
  "ApsimCampbell": "unit.SoilTemperature.xml",
  "DSSAT_ST_standalone": "unit.STEMP.xml",
  "BiomaSurfacePartonSoilSWATC": "unit.SoilTemperatureSWAT.xml"
}

BENCHMARK_MODELS = {
  "Llama3.1", "OpenAI", "Mistral", "Claude"
}

COMPONENTS_DICT = {
  "ApsimCampbell": "unit.SoilTemperature.xml",
  "BiomaSurfacePartonSoilSWATC": "unit.SoilTemperatureSWAT.xml",
  "BiomaSurfaceSWATSoilSWATC": "unit.SoilTemperatureSWAT.xml",
  "DSSAT_EPICST_standalone": "unit.STEMP_EPIC.xml",
  "DSSAT_ST_standalone": "unit.STEMP.xml",
  "Simplace_Soil_Temperature": "unit.STMPsimCalculator.xml",
  "SQ_Soil_Temperature": "unit.CalculateSoilTemperature.xml",
  "Stics_soil_temperature": "unit.Tempprofile.xml"
}

CODE_DICT = {
  "ApsimCampbell": "AP",
  "BiomaSurfacePartonSoilSWATC": "PS",
  "BiomaSurfaceSWATSoilSWATC": "SW",
  "DSSAT_EPICST_standalone": "DE",
  "DSSAT_ST_standalone": "DS",
  "Simplace_Soil_Temperature": "SA",
  "SQ_Soil_Temperature": "SQ",
  "Stics_soil_temperature": "ST"
}

ENERGY_BALANCE_FILES =  {
   "unit.Evapotranspiration.xml" : "Evapotranspiration",
    "unit.Netradiation.xml" : "Net radiation",
    "unit.Netradiationequivalentevaporation.xml" : "Net radiation equivalent evaporation",
    "unit.Penman.xml" : "Penman",
    "unit.Potentialtranspiration.xml" : "Potential transpiration",
    "unit.Priestlytaylor.xml" : "Priestly-Taylor",
    "unit.Ptsoil.xml" : "Ptsoil",
    "unit.Soilevaporation.xml" : "Soil evaporation",
    "unit.Soilheatflux.xml ": "Soil heat flux",
    "unit.Canopytemperature.xml" : "Canopy temperature",
    "unit.Cropheatflux.xml" : "Crop heat flux",
    "unit.Diffusionlimitedevaporation.xml" : "Diffusion limited evaporation"
}

OUTPUT_FILE_COMPONENTS = "../SoilT_comparison.csv"
ENERGY_BALANCE_DIRECTORY = "SQ_EnergyBalance"
OUTPUT_FILE_ENERGY_BALANCE = "../Energy_balance_comparison.csv"
OUTPUT_FILE_BENCHMARK = "../Benchmark.csv"

number_iteration = 10
component_keys = list(COMPONENTS_DICT.keys())
component_values = list(COMPONENTS_DICT.values())
code_keys = list(CODE_DICT.values())
energy_balance_keys = list(ENERGY_BALANCE_FILES.keys())
energy_balance_values = list(ENERGY_BALANCE_FILES.values())
benchmark_keys = list(BENCHMARK_FILES.keys())
benchmark_values = list(BENCHMARK_FILES.values())

#-----------------------------------------------------------------
# Simulation section
#-----------------------------------------------------------------
# Soil temperature section
with open(OUTPUT_FILE_COMPONENTS, 'w', newline='') as f:
  writer = csv.writer(f, delimiter=";")
  writer.writerow([
    "Model", "Code", "Iteration", "Main line code", "Auxiliary line code", "Comments line", "Refactor line",
    "Crop2ML-generated inputs", "AI-generated inputs", "True positive inputs", "False positive inputs", 
    "False negative inputs", "Precision inputs", "Recall inputs", "F1-Score inputs",
    "nb_parameter", "nb_parameter_ok", "nb_state", "nb_state_ok", "nb_rate", "nb_rate_ok", "nb_exogenous", "nb_exogenous_ok",
    "nb_auxiliary", "nb_auxiliary_ok", "nb_int", "nb_int_ok", "nb_double", "nb_double_ok",
    "nb_doublearray", "nb_doublearray_ok", "nb_doublelist", "nb_doublelist_ok", "nb_char", "nb_char_ok",
    "nb_date", "nb_date_ok", "nb_intarray", "nb_intarray_ok", "nb_len", "nb_len_ok",
    "Crop2ML-generated outputs", "AI-generated outputs", "True positive outputs", 
    "False positive outputs", "False negative outputs", 
    "Precision outputs", "Recall outputs", "F1-Score outputs",
    "nb_int_out", "nb_int_ok_out", "nb_double_out", "nb_double_ok_out", "nb_doublearray_out", "nb_doublearray_ok_out", 
    "nb_doublelist_out", "nb_doublelist_ok_out", "nb_char_out", "nb_char_ok_out", 
    "nb_date_out", "nb_date_ok_out", "nb_intarray_out", "nb_intarray_ok_out", "nb_len_out", "nb_len_ok_out",
    "Crop2ML-generated functions", "AI-generated functions", "True positive functions", 
    "False positive functions", "False negative functions", "Precision functions", 
    "Recall functions", "F1-Score functions"
  ])
  for i in range(0, len(component_keys)):
    ref_file = f"../Components_original/{component_keys[i]}/crop2ml/{component_values[i]}"
    for j in range(1, number_iteration + 1):
      test_file = f"../Components/{component_keys[i]}/{j}/{component_values[i]}"
      output_file = f"../Components/{component_keys[i]}/{j}/{os.path.splitext(os.path.basename(test_file))[0]}.csv"
      result = main(ref_file, test_file, output_file)
      result_line = count_line(f"../Components/{component_keys[i]}/", component_values[i], j)                                         
      writer.writerow([component_keys[i]] + [code_keys[i]] + [j] + result_line + result)


# Energy balance section
with open(OUTPUT_FILE_ENERGY_BALANCE, 'w', newline='') as f:
  writer = csv.writer(f, delimiter=";")
  writer.writerow([
    "Model", "Iteration", 
    "Crop2ML-generated inputs", "AI-generated inputs", "True positive inputs", "False positive inputs", 
    "False negative inputs", "Precision inputs", "Recall inputs", "F1-Score inputs",
    "nb_parameter", "nb_parameter_ok", "nb_state", "nb_state_ok", "nb_rate", "nb_rate_ok", "nb_exogenous", "nb_exogenous_ok",
    "nb_auxiliary", "nb_auxiliary_ok", "nb_int", "nb_int_ok", "nb_double", "nb_double_ok",
    "nb_doublearray", "nb_doublearray_ok", "nb_doublelist", "nb_doublelist_ok", "nb_char", "nb_char_ok",
    "nb_date", "nb_date_ok", "nb_intarray", "nb_intarray_ok", "nb_len", "nb_len_ok",
    "Crop2ML-generated outputs", "AI-generated outputs", "True positive outputs", 
    "False positive outputs", "False negative outputs", 
    "Precision outputs", "Recall outputs", "F1-Score outputs",
    "nb_int_out", "nb_int_ok_out", "nb_double_out", "nb_double_ok_out", "nb_doublearray_out", "nb_doublearray_ok_out", 
    "nb_doublelist_out", "nb_doublelist_ok_out", "nb_char_out", "nb_char_ok_out", 
    "nb_date_out", "nb_date_ok_out", "nb_intarray_out", "nb_intarray_ok_out", "nb_len_out", "nb_len_ok_out",
    "Crop2ML-generated functions", "AI-generated functions", "True positive functions", 
    "False positive functions", "False negative functions", "Precision functions", 
    "Recall functions", "F1-Score functions"
  ])
  for i in range(0, len(energy_balance_keys)):
    ref_file = f"../Components_original/{ENERGY_BALANCE_DIRECTORY}/crop2ml/{energy_balance_keys[i]}"
    for j in range(1, number_iteration + 1):
      test_file = f"../Components/{ENERGY_BALANCE_DIRECTORY}/{j}/{energy_balance_keys[i]}"
      output_file = f"../Components/{ENERGY_BALANCE_DIRECTORY}/{j}/{os.path.splitext(os.path.basename(test_file))[0]}.csv"
      result = main(ref_file, test_file, output_file)                                           
      writer.writerow([energy_balance_values[i]] + [j] + result)


#Benchmark section
with open(OUTPUT_FILE_BENCHMARK, 'w', newline='') as f:
  writer = csv.writer(f, delimiter=";")
  writer.writerow([
    "LLM", "Model", 
    "Crop2ML-generated inputs", "AI-generated inputs", "True positive inputs", "False positive inputs", 
    "False negative inputs", "Precision inputs", "Recall inputs", "F1-Score inputs",
    "nb_parameter", "nb_parameter_ok", "nb_state", "nb_state_ok", "nb_rate", "nb_rate_ok", "nb_exogenous", "nb_exogenous_ok",
    "nb_auxiliary", "nb_auxiliary_ok", "nb_int", "nb_int_ok", "nb_double", "nb_double_ok",
    "nb_doublearray", "nb_doublearray_ok", "nb_doublelist", "nb_doublelist_ok", "nb_char", "nb_char_ok",
    "nb_date", "nb_date_ok", "nb_intarray", "nb_intarray_ok", "nb_len", "nb_len_ok",
    "Crop2ML-generated outputs", "AI-generated outputs", "True positive outputs", 
    "False positive outputs", "False negative outputs", 
    "Precision outputs", "Recall outputs", "F1-Score outputs",
    "nb_int_out", "nb_int_ok_out", "nb_double_out", "nb_double_ok_out", "nb_doublearray_out", "nb_doublearray_ok_out", 
    "nb_doublelist_out", "nb_doublelist_ok_out", "nb_char_out", "nb_char_ok_out", 
    "nb_date_out", "nb_date_ok_out", "nb_intarray_out", "nb_intarray_ok_out", "nb_len_out", "nb_len_ok_out",
    "Crop2ML-generated functions", "AI-generated functions", "True positive functions", 
    "False positive functions", "False negative functions", "Precision functions", 
    "Recall functions", "F1-Score functions"
  ])
  for i in range(0, len(benchmark_keys)):
    ref_file = f"../Components_original/{benchmark_keys[i]}/crop2ml/{benchmark_values[i]}"
    for llm in BENCHMARK_MODELS:
      test_file = f"../Benchmark/{llm}/{benchmark_values[i]}"
      output_file = f"../Benchmark/{llm}/{os.path.splitext(os.path.basename(test_file))[0]}.csv"
      result = main(ref_file, test_file, output_file)                                           
      writer.writerow([llm] + [CODE_DICT.get(benchmark_keys[i])] + result)