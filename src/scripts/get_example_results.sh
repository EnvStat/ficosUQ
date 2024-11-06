#!/bin/bash
# download examle results automatically from GitHub
# run this from the project root folder

# Define paths -----------
prjRootDir=$(pwd)

exResDir=results/calibration_example
exFigDir=figures
# url for the release with example results
exResUrl=https://github.com/kkaurila/ficosUQ/releases/download/v1.0.2RC
exResFile=postPredSummary.mat
exPostPred=$exResUrl/$exResFile

exFig3File=example_jointPosterior.png
exFig4File=example_postPred_Uto.png
exFig3=$exResUrl/$exFig3File
exFig4=$exResUlr/$exFig4File

# Download files -----------

# example results
cd $exResDir
if [ ! -f "$exPostPred" ];then
  wget $exPostPred
fi
cd $prjRootDir

# example figures
if [ ! -d "$exFigDir" ];then
  mkdir $exFigDir
fi
cd $exFigDir

if [ ! -f "$exFig3" ];then
  wget $exFig3
fi

if [ ! -f "$exFig4" ];then
  wget $exFig4
fi

cd $prjRootDir
