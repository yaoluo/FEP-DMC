#!/usr/bin/env python3
# This Python script generates PERTURBO input files. The user should specify the calculation type (-c option)

import datetime,os,sys
import numpy as np
import datetime
import argparse
import re
from ast import literal_eval

from yaml import load, dump
try:
   from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
   from yaml import Loader, Dumper


##########################################################################
def set_default(param_name,param_data): 
   if 'default' in param_data[executable][param_name].keys():
      return param_data[executable][param_name]['default']
      
   elif 'typical' in param_data[executable][param_name].keys():
      return param_data[executable][param_name]['typical']

   else:
      sys.exit('Each variable must contain \'default\' or \'typical field\' in yaml file')

##########################################################################
def get_description(description):
      desc = re.sub(r'<a.*?>', '', description, flags=re.DOTALL)
      desc = re.sub(r'<.*?code>', '', desc, flags=re.DOTALL)
      return desc

##########################################################################
def create_arg_namespace(param_data,input_data):
   
   list_calc_mode = input_data.keys()
   
   help_description = 'Generate an exemplified input file for the PERTURBO code'
   
   parser = argparse.ArgumentParser(description=help_description)
   
   parser.add_argument('-c','--calc_mode', help='Calculation mode.', choices=list_calc_mode, required=True)
   parser.add_argument('-i','--input_name', help='Name of the input file.')

   # add all perturbo parameters as arguments, expcept calc_mode
   for key,value in param_data['perturbo'].items():
      if key == 'calc_mode':
         continue
      if value['type'] == 'family':
         continue
      
      # description for the help. remove <code>, </code>, <a>*</a>
      description = get_description(value['description'])
      
      parser.add_argument('--'+key, help=description)

   # add all qe2pert parameters as arguments, expcept repeated_parameters
   for key,value in param_data['qe2pert'].items():
      repeated_parameters = ['prefix','debug']
      if key in repeated_parameters:
         continue
      if value['type'] == 'family':
         continue
      
      # description for the help. remove <code>, </code>, <a>*</a>
      description = get_description(value['description'])
      
      parser.add_argument('--'+key, help=description)

   return parser.parse_args()

##########################################################################
def write_parameter_to_input(finput,param,executable,args,param_data,optional=False):


   if optional:
      if args.__dict__[param] is None:
         str_key = '! ' + str(param)
         param_value = set_default(param,param_data)
      else:
         str_key = ' ' + str(param)
         param_value = args.__dict__[param]
   else:
      str_key = ' ' + str(param)
      if args.__dict__[param] is None:
         param_value = set_default(param,param_data)
      else:
         param_value = args.__dict__[param]
      
   param_type = param_data[executable][param]['type']
   
   #
   # Dimensions
   #
   if 'dimensions' in param_data[executable][param].keys():
      dimensions = int(param_data[executable][param]['dimensions'].replace('(','').replace(')',''))
   else:
      dimensions = 1

   #
   # Units
   #
   if 'units' in param_data[executable][param].keys():
      units = param_data[executable][param]['units']
   else:
      units = None

   if param_type == 'string':
      
      if param_value == "''":
         finput.write('{:20} = \'\'\n'.format(str_key))
      else:
         # Replace prefix if it was specified
         if args.__dict__['prefix'] is not None:
            param_value = re.sub('prefix',args.__dict__['prefix'],param_value)
         finput.write('{:20} = \'{}\'\n'.format(str_key,param_value))
   
   elif param_type == 'logical':
      true_false = param_value
      if type(true_false) is not bool:
         if 'true' in true_false.lower():
            true_false = True
         else:
            true_false = False   

      if true_false: 
         finput.write('{:20} = .true.\n'.format(str_key))
      else:
         finput.write('{:20} = .false.\n'.format(str_key))

   elif dimensions != 1:
      default_tuple = literal_eval(param_value)
      for dim in range(dimensions):
         finput.write('{:20} = {}\n'.format(str_key+'('+str(dim+1)+')',default_tuple[dim]) )
   else: 
      if units:
         finput.write('{:20} = {:<15} ! {}\n'.format(str_key,param_value,units) )
      else:
         finput.write('{:20} = {:}\n'.format(str_key,param_value) )


##########################################################################
##########################################################################

script_path = os.path.dirname(os.path.abspath(__file__))
docs_path   = os.path.abspath(os.path.join(script_path, '..', 'docs'))

param_qe2pert  = os.path.join(docs_path,'input_parameters_qe2pert.yml')
param_perturbo = os.path.join(docs_path,'input_parameters_perturbo.yml')
input_template = os.path.join(docs_path,'input_template.yml')

param_data = {}

# Read the yaml files
with open(param_qe2pert,'r') as stream:
   param_data['qe2pert'] = load(stream,Loader=Loader)

with open(param_perturbo,'r') as stream:
   param_data['perturbo'] = load(stream,Loader=Loader)

with open(input_template,'r') as stream:
   input_data = load(stream,Loader=Loader)

# Parse the command line
args = create_arg_namespace(param_data,input_data)

# Check if the user did not specify a variabe which is not used in a given calc_mode
list_used_param = []

# Parameters that should not be checked
general_arg_list = ['calc_mode', 'input_name', 'prefix']

if 'mandatory' in input_data[args.calc_mode].keys():
   list_used_param = list_used_param + list(input_data[args.calc_mode]['mandatory'].keys())

if 'optional' in input_data[args.calc_mode].keys():
   list_used_param = list_used_param + list(input_data[args.calc_mode]['optional'].keys())

for arg in vars(args):
   if arg in general_arg_list:
      continue
   if args.__dict__[arg] is not None:
      if arg not in list_used_param:
         print('WARNING: '+arg+' input parameter is not used in '+args.calc_mode+' calculation mode.')
##########################################################################

# Write the input
if not args.input_name:
   if args.calc_mode == 'qe2pert':
      input_name = 'qe2pert.in'
   else:
      input_name = 'pert.in'
else:
   input_name = args.input_name

finput = open(input_name,'w')

finput.write('! This input file for PERTURBO was generated by {} script'.format(os.path.basename(__file__))+'\n')
finput.write('! Date: {} '.format(datetime.datetime.now().strftime("%B %d, %Y  %H:%M"))+'\n'*2)

if args.calc_mode == 'qe2pert':
   finput.write('&qe2pert'+'\n')
else:
   finput.write('&perturbo'+'\n')
   
# ===========MANDATORY=================  
finput.write('! ***Mandatory parameters***'+'\n')

if args.calc_mode == 'qe2pert':
   executable = 'qe2pert'
else:
   executable = 'perturbo'
   finput.write('{:20} = \'{}\'\n'.format(' calc_mode',args.calc_mode))

write_parameter_to_input(finput,'prefix',executable,args,param_data)

for key,value in input_data[args.calc_mode]['mandatory'].items():
   if value:
      finput.write('! '+value+'\n')

   write_parameter_to_input(finput,key,executable,args,param_data)
   
# ===========OPTIONAL=================  
if 'optional' in input_data[args.calc_mode].keys():
   finput.write('\n'*2+'! ***Optional parameters***'+'\n')
   for param,value in input_data[args.calc_mode]['optional'].items():
      write_parameter_to_input(finput,param,executable,args,param_data,optional=True)



finput.write('/\n')

finput.close()

print('File '+input_name+' is generated.')

