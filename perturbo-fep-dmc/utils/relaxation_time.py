#!/usr/bin/env python3
# This Python script calculates the relaxation times and scattering rates from the imsigma calculation

import os,sys
import numpy as np
import argparse

description = 'Calculation of the relaxation times from the imsigma PERTURBO calculation'

parser = argparse.ArgumentParser(description=description)

parser.add_argument('-i','--imsigma_file',help='The output file from the imsigma calculation')

#
# parse the command line
#
args = parser.parse_args()

#
# if the imsigma file was not specified, try to find it
#
if args.imsigma_file is None:
   imsigma_files = []
   for filename in os.listdir("."):
      if filename.endswith(".imsigma"):
         imsigma_files.append(filename)
   if len(imsigma_files) == 0:
      sys.exit('Error: No .imsigma files found.')

   elif len(imsigma_files) > 1:
      sys.exit('Error: More then one .imsigma files found in the directory. Specify the file with -i option.')
   else:
      imsigma_file = imsigma_files[0]

else:
   imsigma_file = args.imsigma_file


print('Processing file: '+imsigma_file)

relax_file = 'relaxation_time.dat'
frel = open(relax_file,'w')
frel.write(' # it   ik  ibnd      E(ibnd)(eV)   Relaxation time(in fs)   Scattering rate (in THz)\n')

# Constants
ryd2mev  = 13.605698066 * 1.0E3 # Rydberg to meV
timeunit = 4.8377687E-17        # Rydberg atomic unit t0 in s

header_block = True
with open(imsigma_file,'r') as fims:
   # Skip the header block
   while True:
      line = fims.readline()
      if not line.strip().startswith('#'):
         break
   
   # Process all the lines after the header block
   while line != "": #for line in fims:
      # comments in the "body" of file
      #line = fims.readline()
      if line.strip().startswith('#'):
         frel.write(line)
      else:
         parsing = line.split()
         if len(parsing) != 5:
            sys.exit('Error: Wrong format of the '+imsigma_file+' file.')
         else:
            imsigma = float(parsing[4])
            # convert to Ryd
            imsigma = imsigma / ryd2mev 
            
            # calculate the relaxation time
            rel_time = ( timeunit * 1.0E15 ) / ( imsigma * 2.0)

            # keep all the fields in the line untouched. Only replace the last field (imsigma) by the relaxation time
            new_line = line.rsplit(' ',1)[0] + (' {:23.16e}'*2).format(rel_time,1.0 / rel_time * 1.0E3) + '\n'
            frel.write(new_line)

      line = fims.readline()

frel.close()

print('The '+relax_file+' file was created.')
