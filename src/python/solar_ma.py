# Function for smoothing solar radiation by moving average
import h5py
import numpy as np
import os
import sys

def solar_ma(n=1):
    fname = "hd_files_ma_%d.hdf5" % n
    cmd = "cp hd_files.hdf5 hd_files_ma_%d.hdf5" % n
    os.system(cmd)
    with h5py.File('hd_files.hdf5','r') as f:
        with h5py.File(fname,'r+') as f_ma:
            sol_set = f['Solar_radiation']
            sol_set_ma = f_ma['Solar_radiation']
            n_sol_set = len(sol_set)
            for i in range(n_sol_set):
                ind = range(i-n,min(n_sol_set-1,i+n))
                sol_set_ma[i] = np.mean(sol_set[ind])

if __name__ == "__main__":
    n = int(sys.argv[1])
    solar_ma(n)

    
