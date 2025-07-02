from openai import OpenAI
import os
from pathlib import Path
import urllib.request, urllib.error, urllib.parse
import json
from pprint import pprint

#CONFIGURATION
system_instruction_path = "./Instructions_crop2ML.txt"
main_file = ""
folder_path = ""
output_path = ""
api_key_OpenAI = ""
model = "gpt-4.1"
  
 
# Function to determine the programming language based on file extension
# This function checks the file extension and returns the corresponding language name.
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
  

# Function to extract text from a file
# This function reads the content of a file and returns it as a string.
def extract_text(file_path):
  with open(file_path, "r", encoding="utf-8") as file:
    return file.read()
  

# Function to send instructions and prompt to OpenAI's GPT model
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
  return(response.choices[0].message.content)


# Main execution block
# This block extracts the system instructions, determines the programming language, counts the number of files in the specified folder, and constructs a prompt for the GPT model.
# It then sends the prompt to the model and saves the response as a JSON file. 
if __name__ == "__main__":
  instructions = extract_text(system_instruction_path)
  language_name = language(folder_path + main_file)
  nb_files = len([f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))])

  if nb_files < 1:
    raise ValueError("There should be at least 1 file in the folder to run this script.")
  prompt = f"You will be given {nb_files} crop model source code files, each written in {language_name}. Each file will be provided in a separate code block and will be labeled with its filename."
  prompt += f"\n\nThe main crop model source code file is called {main_file}, and its content is as follows: {extract_text(folder_path + main_file)}"
  
  if nb_files > 1:
    number_file = 2
    prompt += "\n\nAdditional source files are as follows:"
    for file in os.listdir(folder_path):
      file_path = os.path.join(folder_path, file)
      if os.path.isfile(file_path):
        prompt += f"\n\nSource code file {number_file} is {file}, and its content is as follows:\n"
        extract_text(file_path)
        number_file += 1

  response = send_to_gpt(instructions, prompt, api_key_OpenAI, model)

if response.startswith("```json"):
  response = response[7:].lstrip()
if response.endswith("```"):
  response = response[:-3].rstrip()

  base, _ = os.path.splitext(main_file)
  json_path = base + ".json"
  output_file = os.path.join(output_path, os.path.basename(json_path))
  json_obj = json.loads(response)

  with open(output_file, "w", encoding="utf-8") as f:
    json.dump(json_obj, f, ensure_ascii=False, indent=4)
