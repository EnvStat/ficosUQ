#!/bin/bash
# download example results from GitHub
# run this from the project root folder

# Define paths -----------
prjRootDir=$(pwd)

exResDir=results/calibration_example
exFigDir=figures
# url for the release with example results
exResUrl=https://github.com/kkaurila/ficosUQ/releases/download/v1.0.2RC
exPostPredFile=postPredSummary.mat
exPostPredUrl=$exResUrl/$exPostPredFile
# chla scenario summary
exChlaSumFile=example_scenarioPredictionTable.mat
exChlaSumUrl=https://github.com/kkaurila/ficosUQ/releases/download/v1.0.3/$exChlaSumFile

exFig3File=example_jointPosterior.png
exFig4File=example_postPred_Uto.png
exFig3Url=$exResUrl/$exFig3File
exFig4Url=$exResUrl/$exFig4File

# Download files -----------

# example results
cd $exResDir
if [ ! -f "$exPostPredFile" ];then
  wget $exPostPredUrl
fi
if [ ! -f "$exChlaSumFile" ];then
  wget $exChlaSumUrl
fi
cd $prjRootDir

# example figures
if [ ! -d "$exFigDir" ];then
  mkdir $exFigDir
fi
cd $exFigDir

if [ ! -f "$exFig3File" ];then
  wget $exFig3Url
fi

if [ ! -f "$exFig4File" ];then
  wget $exFig4Url
fi

cd $prjRootDir
