#!/bin/bash
# Script for automatically downloading data
# launch this from ficosUQ project folder
# The download may fail if Zenodo is busy,
# in this case the download is retried 3 times with 10 second wait in between

dataRec=14028484
hdf5Rec=14028433
dataDir=data
zRetries=3
zRetryWait=10

# https://pypi.org/project/zenodo-get/
pyBins=src/python/pyEnv/venv/bin
zGet=$pyBins/zenodo_get

if [ ! -f "$zGet" ]; then
  # install zenodo-get
  if [ ! -d "$pyBins" ]; then
    echo "Python virtual environment not found, please install it first with src/python/python_install.sh" 1>&2
    exit 64
  fi
  pipPath=$pyBins/pip
  $pipPath install zenodo-get --require-virtualenv
fi

# download files
zGetCmd="$zGet -o $dataDir -R $zRetries -p $zRetryWait"
$zGetCmd -r $dataRec
$zGetCmd -r $hdf5Rec
