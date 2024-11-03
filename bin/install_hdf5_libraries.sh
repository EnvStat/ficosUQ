#!/bin/bash

# Install hdf5 libraries needed by FICOS simulator

# run this from the project root folder ('ficosUQ')

ficosUQroot=$(pwd)
binDir=$ficosUQroot/bin
hdf5Dir=$binDir/hdf5

if [ ! -d "$hdf5Dir" ]; then
  mkdir $hdf5Dir
fi

cd $hdf5Dir

if [ ! -d "hdf5-1.8.19" ]; then
  if [ ! -f "hdf5-1.8.19.tar.gz" ]; then
    # see also: https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/obtain51819.html
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/hdf5-1.8.19.tar.gz
  fi


  tar zxvf hdf5-1.8.19.tar.gz
  rm hdf5-1.8.19.tar.gz
fi

cd hdf5-1.8.19

./configure CC=gcc-9 FC=gfortran-9 -enable-fortran -enable-fortran2003 FCFLAGS=-fdefault-integer-8  
make
make install
