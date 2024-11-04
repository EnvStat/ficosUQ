#!/bin/bash
# main installation script

# Setting up git submodules for gpstuff and BrewerMap ------
gpstuffRef=114937ec0a201306489a66cbba38283e722fb998
BrewerMapRef=1773eb111abe7db7a4197e09a53ca42e223f28eb
if [ -d ".git" ]
then
  git submodule init
  git submodule update
else
  # Not installed through git (e.g. through an archive in Zenodo)
  # cloning the submodules instead and checking out specific versions
  cd src/submodules
  # gpstuff -------------------------------------
  git clone https://github.com/gpstuff-dev/gpstuff.git gpstuff
  cd gpstuff
  git checkout $gpstuffRef
  cd ..

  # BrewerMap
  git clone https://github.com/DrosteEffect/BrewerMap.git BrewerMap
  cd BrewerMap
  git checkout $BrewerMapRef
  cd ../../..
fi

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


# Install python virtual environment
src/python/python_install.sh

# Download data from Zenodo
src/python/download_data.sh
