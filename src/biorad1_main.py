
from biorad1_computing import ComputingModel
import os


def system_environment_setting():


   print("Biorad1 version 0.9 beta")
   print('Working directory: ' + os.getcwd() + '"')


def execute_Biorad1():
   input_file = "../inputs/UZ_conf_file.yaml"
   system_environment_setting()
   print(60 * "-")
   model = ComputingModel()
   if not model.load_yaml(input_file): 
      print("Error of yaml loading.")
      return False
   print(60 * "-")
   if not model.load_data_model():
      print("Error of data model loading.")
      return False
   print(60 * "-")
   if not model.compute():
      print("Error of computing.")
      return False
   return True
