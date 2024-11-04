#!/bin/bash
# launch this from ficosUQ project folder
# local python installation set to version 3.8.10
# following guide in:
# https://help.dreamhost.com/hc/en-us/articles/115000702772-Installing-a-custom-version-of-Python-3

# custom python 3.8.10 installation ---------------------------

localPyDir='pyEnv'

cd src/python

if [ ! -d "$localPyDir" ]; then
  mkdir $localPyDir
fi
cd $localPyDir
localPyFullPath=$(pwd)

pySrcDir="Python-3.8.10_source"
pyInstallDir="Python-3.8.10"
pyInstallFullDir=$localPyFullPath/$pyInstallDir
if [ ! -d "$pyInstallDir" ]; then
  
  if [ ! -f "Python-3.8.10.tgz" ]; then
    wget https://www.python.org/ftp/python/3.8.10/Python-3.8.10.tgz
  fi
  
  if [ ! -d "$pySrcDir" ]; then
    mkdir $pySrcDir
    tar zxvf Python-3.8.10.tgz --directory $pySrcDir
    rm Python-3.8.10.tgz
  fi
  
  mkdir $pyInstallDir
  
  cd $pySrcDir/Python-3.8.10
  ./configure --prefix=$pyInstallFullDir
  make
  make install
  
  cd $localPyFullPath
fi

customPyRoot=$pyInstallFullDir/bin
customPy3=$customPyRoot/python3
customPip3=$customPyRoot/pip3


# virtual environment -------------------------------


customVirtualenv=$customPyRoot/virtualenv

if [ ! -f "$customVirtualenv" ]; then
  # install virtual environment
  $customPip3 install virtualenv

  # create virtual environment
  $customVirtualenv -p $customPy3 $localPyFullPath/venv
fi

# paths to python3 and pip3 within virtual environment
venvDir=$localPyFullPath/venv/bin
py3venv=$venvDir/python3
pip3venv=$venvDir/pip3
pipvenv=$venvDir/pip

# install python dependencies within virtual environment ---------

installCmd="$pip3venv install --require-virtualenv --force-reinstall"
$installCmd h5py==2.10.0
# install numpy last, as installing h5py afterwards may lead to installing newer numpy
$installCmd numpy==1.17.4

