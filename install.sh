#!/bin/bash
# main installation script


# Update submodules
# -----------------------------
# Updates submodules to match what the main (superproject) expects.
# Essentially ensures that external libraries from github are fixed to the commits that
# are known to work in this project.
# see 'man gitsubmodules(7)' and 'man git-submodule'
git submodule init
git submodule update

# Install matlab dependencies --------

matlab -nodisplay -nodesktop -r "run('matlabSetup.m'); exit;"

# HDF5 libraries ------

if [ ! -f "bin/hdf5/hdf5-1.8.19/hdf5/bin/h5fc" ]; then
  bin/install_hdf5_libraries.sh
fi

# FICOS simulator ------
if [ ! -f "bin/ficos/wqficos" ]; then
  cd bin/ficos
  make
  cd ../..
fi

# download hdf5 data -----
# zenodo_get ...
# ---


# Install python virtual environment
src/python/python_install.sh
