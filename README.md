# ficosUQ

## Usage

### Setup 

When first using this project, check that you have installed the required compilers and Matlab toolboxes (see the Dependencies section below) and then run the installation script 'install.sh' on a terminal from the repository root folder ('ficosUQ'). This script will install all the remaining external libraries.

### FICOS Calibration

After installing the required libraries, the ficos calibration code can be run by starting a Matlab session in the 'ficosUQ' folder and running 'ficos_calibration();', which will calibrate the FICOS parameters with all the default options. This may take a long time depending on the amount of cores available. The calibration script ends with posterior sampling for the calibrated parameters.

Note: the version of FICOS being calibrated uses the disk quite heavily as it reads and writes frequently from/to the .hdf5 input and result files. You may want to reduce the amount of simulations run in parallel by setting the 'cores' option to 'ficos_calibration', e.g. run 'ficos_calibration('cores', [2, 5]);' to run between 2 to 5 simulations in parallel (default is [2, 20]).

### Posterior inference

#### Parameter posterior

The FICOS calibration script above will automatically sample from the posterior distribution for the simulator parameters. The samples will be saved in the 'calibration_result.mat' file in the variables 'ws_rt'. The posterior can be visualized with the function 'plotPostJointMargDens'. When called without arguments, this will plot the example results in 'results/calibration_example'. To plot different results, call this function with the path to the 'calibration_results.mat' produced by the calibration script. 

#### Posterior predictive simulations

You can also use the For posterior sampling, start a Matlab session in the 'ficosUQ' folder and run 'postPredSimulations'. By default this will start posterior sampling from an example result file. You may sample instead from your own results by calling this function with the 'sampleFile'-argument, e.g. postPredSimulations('sampleFile','<path_to_your_results_file.mat>');. 

#### Posterior predictive summaries

After all posterior predictive simulations have finished, their results can be collected and summarised with 
 >> postPredTimeseries(postPredPath);

where postPredPath is the root folder for the posterior predictive simulation results. 
 

## Project structure

This project folder is structured as follows:

- ficosUQ/
 - bin/
  - ficos/         : FICOS simulator source code
   - README.md     : separate README for the FICOS simulator
  - hdf5/          : external hdf5 binaries needed to run ficos
 - data/           : calibration data and .hdf5 input files for the simulator
 - src/            : Source code for this project organized into sub folders by topic
 - README.md       : this file
 - startup.m       : startup script run by Matlab
 - install.sh      : script for installing all the dependencies needed by the calibration code
 - matlabSetup.m   : matlab script used as part of the install.sh to install Matlab specific dependencies

## Dependencies

### Compilers

- gcc-9 and gfortran-9
  - simulator compilation has been tested with versions 9.4.0 and 9.5.0 for both gcc and gfortran
  - on Ubuntu, install both with 'sudo apt-get install gfortran-9', which should also install gcc-9 of the same version.
  - version can be checked with 'gfortran-9 --version' and 'gcc-9 --version' respectively
  
### Matlab

FICOS Calibration and posterior predictions have been run on a server with Matlab R2022b, while the posterior summaries and figures were done on a laptop with Matlab R2024a.

The main calibration script requires the following Matlab toolboxes:
 - Control System Toolbox
 - Optimization Toolbox
 - Statistics and Machine Learning Toolbox
 - Curve Fitting Toolbox
 - Parallel Computing Toolbox

Toolbox dependencies were checked with matlab.codetools.requiredFilesAndProducts(). See also the script in src/util/listLocalDependencies.m .

### External libraries

#### Installation script

All of the external libraries below can be installed by running the install.sh script in the project root directory. The script assumes you have installed the compilers gcc-9 and gfortran-9 above. 

#### hdf5 1.8.19

The hdf5 library version 1.8.19 and its associated fortran compiler script is installed via the installation script in bin/install_hdf5_libraries.sh.
This requires the gcc-9 and gfortran-9 compilers mentioned above.

#### Python 3.8.10

The python functions in src/python have been tested with Python 3.8.10 with using the following libraries:
- numpy 1.17.4
- h5py  2.10.0

Python 3.8.10 as well as the above libraries are installed into a virtual environment using the script in src/python/python_install.sh.


## Data

The data used for calibrating FICOS is stored in two Zenodo repositories, 10.5281/zenodo.14028484 for the calibration data and 10.5281/zenodo.14028433 for the .hdf5 input files. All of the following files need to be placed in the ficosUQ/data folder:

- intensiveStationsData.xlsx
- hd_files.hdf5
- hd_files_ma_5.hdf5
- loading.hdf5
- loading_Aurajoki.hdf5


