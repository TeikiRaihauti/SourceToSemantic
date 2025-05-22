from openai import OpenAI
import os
from pathlib import Path
import urllib.request, urllib.error, urllib.parse
import json
from pprint import pprint

#CONFIGURATION
system_instruction_path = "C:/Users/raihauti/Documents/These/References/Instructions.txt"
file_path = "C:/Users/raihauti/Documents/These/References/SourceCode/xsamara/SOLAR.for"
folder_path = "C:/Users/raihauti/Documents/These/References/SourceCode/dssat-csm-os"
output_path = "C:/Users/raihauti/Documents/These/References/TestSource/"
api_key_OpenAI = ""
model = "gpt-4o-mini"
REST_URL = "https://services.agroportal.lirmm.fr"
api_key_AgroPortal = ""
ontologies = "AGROVOC,CMP"
  
 

def language(file_path):
  fichier = Path(file_path)
  if fichier.suffix == ".py":
    return "Python"
  elif fichier.suffix == ".java":
    return "Java"
  elif fichier.suffix == ".cs":
    return "C#"
  elif fichier.suffix == ".cpp":
    return "C++"
  elif fichier.suffix == ".for":
    return "Fortran"
  elif fichier.suffix == ".f90":
    return "Fortran90"
  else:
    return "Unknown"
  

def afficher_arborescence_avec_fichier_cible(folder_path, file_path):
  folder_path = os.path.abspath(folder_path)
  file_path = os.path.abspath(file_path)

  arbo = []
  for folder, sub_folder, files in os.walk(folder_path):
    niveau = folder.replace(folder_path, '').count(os.sep)
    indent = '│   ' * niveau + '├── '
    arbo.append(f"{indent}{os.path.basename(folder)}/")

    for f in files:
      f_path = os.path.join(folder, f)
      if os.path.abspath(f_path) == file_path:
        indent_fichier = '│   ' * (niveau + 1) + '├── '
        arbo.append(f"{indent_fichier}{f} <== TARGET FILE")

  return "\n".join(arbo)


def extract_text(file_path):
  with open(file_path, "r", encoding="utf-8") as file:
    return file.read()
  

def send_to_gpt(instructions, prompt, api_key, model):
  client = OpenAI(api_key = api_key)

  response = client.chat.completions.create(
    model=model,
    messages=[
        {"role": "system", "content": instructions},
        {"role": "user", "content": prompt}
    ]
  )
  return(response.choices[0].message.content)


def get_json(url):
  opener = urllib.request.build_opener()
  opener.addheaders = [('Authorization', 'apikey token=' + api_key_AgroPortal)]
  return json.loads(opener.open(url).read())


def return_annotations(annotations):
  annotations_result = []
  for result in annotations:
    try:
      get_json(result["annotatedClass"]["links"]["self"])
    except urllib.error.HTTPError:
      print(f"Error retrieving {result['annotatedClass']['@id']}")
      continue

    ontology = result["annotatedClass"]["links"]["ontology"].split("/")[-1]
    score = result["score"]
    for annotation in result["annotations"]:
      annotations_result.append((ontology, annotation["text"], score))
  return annotations_result


def get_description(text):
  if text.startswith("```json"):
    text = text[7:].lstrip()
  if text.endswith("```"):
    text = text[:-3].rstrip()
  json_obj = json.loads(text)
  for function in json_obj.get('functions', []):
    for input in function.get('inputs', []):
      description = input.get('description', '')
      agroportal_result = get_json(REST_URL + "/annotator?ontologies=" + urllib.parse.quote(ontologies) + "&text=" + urllib.parse.quote(description) + "&score=old&lemmatize=t")
      annotations = return_annotations(agroportal_result)
      annotations_json = [
        {"ontology": ontologie, "text": texte, "score": score}
        for ontologie, texte, score in annotations
        ]
      input["annotations"] = annotations_json

    for output in function.get('outputs', []):
      description = output.get('description', '')
      agroportal_result = get_json(REST_URL + "/annotator?ontologies=" + urllib.parse.quote(ontologies) + "&text=" + urllib.parse.quote(description) + "&score=old&lemmatize=t")
      annotations = return_annotations(agroportal_result)
      annotations_json = [
        {"ontology": ontologie, "text": texte, "score": score}
        for ontologie, texte, score in annotations
        ]      
      output["annotations"] = annotations_json
      
    description = function.get('description', '')
    agroportal_result = get_json(REST_URL + "/annotator?ontologies=" + urllib.parse.quote(ontologies) + "&text=" + urllib.parse.quote(description) + "&score=old&lemmatize=t")
    annotations = return_annotations(agroportal_result)
    annotations_json = [
      {"ontology": ontologie, "text": texte, "score": score}
      for ontologie, texte, score in annotations
      ]
    function["annotations"] = annotations_json
  return json_obj


if __name__ == "__main__":
  instructions = extract_text(system_instruction_path)
  language_name = language(file_path)
  text = extract_text(file_path)
  arbo = afficher_arborescence_avec_fichier_cible(folder_path, file_path)
  prompt = f"Follow the instructions based on the contents of this {language_name} file. Here is the directory structure for context and the target file is marked accordingly: \n\n{arbo} \n\n Here is the content:\n\n{text}\n\n"
  response = send_to_gpt(instructions, prompt, api_key_OpenAI, model)
  json_obj = get_description(response)

  base, _ = os.path.splitext(file_path)
  json_path = base + ".json"
  output_file = os.path.join(output_path, os.path.basename(json_path))

  with open(output_file, "w", encoding="utf-8") as f:
    json.dump(json_obj, f, ensure_ascii=False, indent=4)



