#!/usr/bin/env python2

#-------------------------------------------------------------------------------
# File: merge_imsigma.py
# Author: jjchou <jjchou.comphy@gmail.com>
# Create Date: 2017-06-10 13:10:47
#
# Comment:
#-------------------------------------------------------------------------------
import numpy as np

class imsigma(object):
   """
   parse prefix.imsigma or prefix.imsigma_mode files.
   """
   def __init__(self, imsigma_file):
      with open(imsigma_file, 'r') as inpf:
         #read head
         self.head = []
         for i in range(4):
            self.head.append(inpf.readline())
         data = self.head[-1].strip().split()
         self.nkpt = np.int( data[ 2] )
         self.nbnd = np.int( data[ 4] )
         self.ntmp = np.int( data[-3] )
         if 'NO.modes' in self.head[-1]:
            self.nmod = np.int( data[-1] )
         else:
            self.nmod = 1 

         #read scattering rate
         self.t_head = []
         self.rates = np.zeros( (self.ntmp, self.nkpt, self.nbnd, self.nmod), np.float)
         self.eig   = np.zeros( (self.nkpt, self.nbnd), np.float )

         for it in range(self.ntmp):
            head = []
            for i in range(5):
               head.append( inpf.readline() )
            self.t_head.append( head )
            
            for ik in range(self.nkpt):
               for ib in range(self.nbnd):
                  for imod in range(self.nmod):
                     data = inpf.readline().strip().split()
                     self.rates[it,ik,ib,imod] = np.float( data[-1] )
                     if it == 0 and imod == 0 : self.eig[ik, ib] = np.float(data[3])
               #skip separator #----
               if self.nmod > 1: inpf.readline()

   def output(self, output):
      outf =  open(output,'w')
      #write head
      outf.writelines(self.head)
      #write scattering rate
      for it in range(self.ntmp):
         outf.writelines(self.t_head[it])
         for ik in range(self.nkpt):
            for ib in range(self.nbnd):
               if self.nmod > 1:
                  for imod in range(self.nmod):
                     ol = "{0:3d}  {1:6d}  {2:4d}  {3:12.6f}  {4:4d}  {5:23.16E}"\
                        .format(it+1, ik+1, ib+1, self.eig[ik,ib], \
                                 imod+1, self.rates[it, ik, ib, imod])
                     outf.writelines(ol+'\n')
               else:
                  ol = "{0:3d}  {1:6d}  {2:4d}  {3:12.6f}  {4:23.16E}"\
                     .format(it+1, ik+1, ib+1, self.eig[ik,ib], \
                     self.rates[it, ik, ib, 0])
                  outf.writelines(ol+'\n')
            if(self.nmod > 1): outf.writelines('#-----------------------\n')
      outf.close()


"""
Merge multiply imsigma files
"""
if __name__ == '__main__':
   from optparse import OptionParser
   usage = "usage: %prog [options] prefix.imsigma-1 prefix.imsigma-2"
   parser = OptionParser(usage=usage)
   parser.add_option('-w', '--weight', type="float", nargs=3, dest='wgt', default=(1.0, 1.0, 2.0),\
         help="(w1, w2, tot): weight factors for each imsigma file w1/tot, w2/tot.")
   parser.add_option('-o', '--output', type="string", dest='opf', default='merge.imsigma',\
      help="output filename")

   options, args = parser.parse_args()
   # number of files to merge
   numf = len(args)
   if numf != 2: raise NameError("Can't merge one file. Needs two files.")

   wfactor = np.array( [ options.wgt[0]/options.wgt[2],  options.wgt[1]/options.wgt[2] ] )

   #parse the first imsigma file
   imsga = imsigma( args[0] )
   # multiply the weight factor
   imsga.rates = imsga.rates*wfactor[0]
      
   tmp_imsga = imsigma( args[1] )
      
   # check consistence
   if tmp_imsga.head[3].strip() != imsga.head[3].strip():
      raise NameError( 'Mismatch:\n' + tmp_imsga.head[3] + imsga.head[3] )
   for it in range(imsga.ntmp):
      if tmp_imsga.t_head[it][1].strip() != imsga.t_head[it][1].strip():
         raise NameError( 'Mismatch:\n' + tmp_imsga.t_head[it] + imsga.t_head[it] )
   
   imsga.rates += tmp_imsga.rates*wfactor[1]

   #output merged files
   imsga.output(options.opf)
