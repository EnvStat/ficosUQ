from __future__ import print_function
import h5py
import numpy
import os
import sys

# Janne Ropponen / SYKE 10.6.2020
#
# Changelog:
# 2023-09-14: Addded 'r+' argument to h5py.File as the code could not write
#             results otherwise on my laptop - Karel Kaurila
# 2020-06-10: First version

### Get arguments
if len(sys.argv)!=5:
   print("Error: Invalid arguments")
   print("Changes given nutrients (N, P) in all blocks by multiplier.")
   print("Usage: python changeloads.py <loadfile> <loadname> <N-multiplier> <P-multiplier>")
   print("  Eg.: python changeloads.py loading.hdf5 atmdep 1.5 1.0")
   print("       Multiplies all loads beginning with atmdep by a factor of 1.5 for N, 1.0 for P")
   print("  Eg.: python changeloads.py loading.hdf5 -atmdep 2.0 0.0")
   print("       Multiplies all loads EXCEPT those beginning with atmdep by a factor of 2 for N, zero for P")
   print("Note: Negative multipliers are not allowed.")
   sys.exit(1)

loadtochange = sys.argv[2]
multiplierN = float(sys.argv[3])
multiplierP = float(sys.argv[4])

if multiplierN<0.0 or multiplierP<0.0:
   print("Error: negative multiplier")
   sys.exit(1)

# 
# root = h5py.File(sys.argv[1])
root = h5py.File(sys.argv[1],'r+')

reverse = False
if loadtochange[0] == '-':
   reverse = True
   loadtochange = loadtochange[1:]

changedloads = 0


for blockname in root['Blocks'].keys():
   for load in root['Blocks/'+blockname+'/Loads'].keys():
      if ((not reverse and load[0:len(loadtochange)] == loadtochange)
        or (reverse and load[0:len(loadtochange)] != loadtochange)):
           loaddata = root['Blocks/'+blockname+'/Loads/'+load]
           loaddata[:] *= numpy.array([multiplierN,multiplierP,multiplierN,multiplierP])
           print("Modified load",load,"in block",blockname)
           changedloads += 1

root.close()

if changedloads == 0:
   print("Load not found in dataset.")
else:
   print("Modified a total of", changedloads, "loads.")
