# ficosUQ
[![DOI](https://zenodo.org/badge/874713046.svg)](https://zenodo.org/badge/latestdoi/874713046)

This repository contains all the code needed to reproduce the calibration and uncertainty quantification of the Finnish Coastal Nutrient Load Model (FICOS) as described in [https://arxiv.org/abs/2410.02448](https://arxiv.org/abs/2410.02448) 

## Usage

The preferred way to use this repository is by cloning it from GitHub, since git and GitHub have been used to develop this project. This repository is also available as an archive on Zenodo, see the link in the DOI badge above. 

### Setup 

When first using this project, check that you have installed the required compilers and Matlab toolboxes (see the Dependencies section below) then open a terminal on the repository root folder (referred to as *ficosUQ* in this README) and run the installation script: 
>install.sh

This script installs all the remaining external libraries and downloads the required data from Zenodo.

### FICOS Calibration

After installing the required libraries, the ficos calibration code can be run by starting a Matlab session in this repository's root folder and running 
>ficos_calibration();

which will calibrate the FICOS parameters with all the default options. This may take a long time depending on the amount of cores available. The calibration script ends with posterior sampling for the calibrated parameters.

#### Note on hard disk usage

The version of FICOS being calibrated uses the disk quite heavily as it reads and writes frequently from/to the .hdf5 input and result files. You may want to reduce the amount of simulations run in parallel with the 'cores' option. For example, you can set the amount of parallel simulations to be between 2 to 5 (default is between 2 to 20) with
>ficos_calibration('cores', [2, 5]);

### Posterior inference

#### Parameter posterior

The FICOS calibration script above will automatically sample from the posterior distribution for the simulator parameters. The samples will be saved in the 'calibration_result.mat' file in the variables 'ws_rt'. The posterior can be visualized with the function 'plotPostJointMargDens'. When called without arguments, this will plot the example results in 'results/calibration_example'. To plot different results, call this function with the path to the 'calibration_results.mat' produced by the calibration script:
>plotPostJointMargDens('<path_to_results_folder>/calibration_results.mat')

#### Posterior predictive simulations

You can also use the For posterior sampling, start a Matlab session in the 'ficosUQ' folder and run 'postPredSimulations'. By default this will start posterior sampling from an example result file. You may sample instead from your own results by calling this function with the named 'sampleFile' argument, e.g. 
>postPredSimulations('sampleFile','<path_to_results_folder>/calibration_results.mat');. 

#### Posterior predictive summaries

After all posterior predictive simulations have finished, their results can be collected and summarised with 
 > postPredTimeseries(postPredPath);

where postPredPath is the root folder for the posterior predictive simulation results. 
 

## Project structure

This project folder is structured as follows:

* ficosUQ/ (project root)
    * README.md (this file)
    * **install.sh** (script for installing all dependencies) 
    * bin/
        * ficos/
            * (FICOS simulator source code)
            * README.md
                * (separate README for the FICOS simulator, **includes license information for FICOS**)
        * hdf5/
            * (external hdf5 binaries needed to run ficos)
    * data/
        * (calibration data and .hdf5 input files for the simulator, see *Data* section below)
    * src/
        * **ficos_calibration.m**
            * (Matlab function for launching the main calibration script)
        * python
            * python_install.sh
                * (script for installing python dependencies in a virtual environment)
            * download_data.sh
                * (script for downloading required data using Python library *zenodo_get*)
            * (python source code used in this project)
        * submodules
            * (git sub modules for external Matlab packages *gpstuff* and *BrewerMap*)  
        * (Other folders containing Matlab source code for this project, organized by topic)
    

## Dependencies

### Compilers

#### gcc and gfrotran

The FICOS simulator and the hdf5 libraries have been compiled and tested using versions **9.4.0** and **9.5.0** for both *gcc* and *gfortran*. The script for compiling FICOS assumes compatible versions of these compilers are installed as *gcc-9* and *gfortran-9*. On Ubuntu, both compilers these can be installed with
>sudo apt-get install gfortran-9

For other Linux distributions, you may need to replace 'apt-get' above with another package manager (e.g. 'yum', 'dnf').

You can check the version of these compilers with 
>gfortran-9 --version
>gcc-9 --version

If you have installed the above compilers in some other manner, you may need to change the 'CC' and 'FC' options in the ./configure in 'bin/install_hdf5_libraries.sh'. The compiler configuration may be updated to be more flexible in a future release.

  
### Matlab

FICOS Calibration and posterior predictions have been run on a server with Matlab R2022b, while the posterior summaries and figures were done on a laptop with Matlab R2024a.

The main calibration script requires the following Matlab toolboxes:
 - Control System Toolbox
 - Optimization Toolbox
 - Statistics and Machine Learning Toolbox
 - Curve Fitting Toolbox
 - Parallel Computing Toolbox

Toolbox dependencies were checked with 
>matlab.codetools.requiredFilesAndProducts('src/ficos_calibration.mat'). Most of the Matlab toolboxes above are required by the gpstuff package used for the Gaussian Process emulator. Parallel Computing Toolbox is needed for running multiple simulations in parallel.

See also 'src/util/findDependencies.m'.

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

Python 3.8.10 as well as the above libraries are installed into a virtual environment using the following script 
>src/python/python_install.sh.


## Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14028484.svg)](https://doi.org/10.5281/zenodo.14028484) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14028433.svg)](https://doi.org/10.5281/zenodo.14028433)

The data used for calibrating FICOS is stored in two Zenodo repositories: 10.5281/zenodo.14028484 for the calibration data and 10.5281/zenodo.14028433 for the .hdf5 input files. All of the following files need to be placed in the ficosUQ/data folder:

- intensiveStationsData.xlsx
- hd_files.hdf5
- hd_files_ma_5.hdf5
- loading.hdf5
- loading_Aurajoki.hdf5

### Automatic download

The main installation script attempts to download all of the above files automatically. This may fail is the Zenodo servers are busy. In this case you may try again later with the following command (executed on the terminal on the project root folder):
>src/python/download_data.sh

### Manual download

You may also download the files manually from the following DOIs:
- https://doi.org/10.5281/zenodo.14028484
- https://doi.org/10.5281/zenodo.14028433
