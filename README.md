# ficosUQ
[![DOI](https://zenodo.org/badge/874713046.svg)](https://zenodo.org/badge/latestdoi/874713046)

This repository contains all the code for implementing the methods described in _Bayesian calibration and uncertainty quantification for a large nutrient load impact model_, published online in Ecological Informatics [https://doi.org/10.1016/j.ecoinf.2024.102976](https://doi.org/10.1016/j.ecoinf.2024.102976) (open access). 

## Usage

The preferred way to use this repository is by cloning it from GitHub, since git and GitHub have been used to develop this project. This repository is also available as an archive on [Zenodo](https://doi.org/10.5281/zenodo.14031270). Note that this link always points to the **most recent release**, which may not match *this* version.

### Setup 

When first using this project, check that you have installed the required compilers and Matlab toolboxes (see the Dependencies section below) then open a terminal on the repository root folder (referred to as *ficosUQ* in this README) and run the installation script: 
`./install.sh`

This script installs all the remaining external libraries and downloads the required data from Zenodo.

### FICOS Calibration

After installing the required libraries, the ficos calibration code can be run by starting a Matlab session in this repository's root folder and running `ficos_calibration();`

which will calibrate the FICOS parameters with all the default options. This may take a long time depending on the amount of cores available. The calibration script ends with posterior sampling for the calibrated parameters.

#### Note on hard disk usage

The version of FICOS being calibrated uses the disk quite heavily as it reads and writes frequently from/to the .hdf5 input and result files. You may want to reduce the amount of simulations run in parallel with the named `'cores'` option. For example, you can set the amount of parallel simulations to be between 2 to 5 (default is between 2 to 20) with `ficos_calibration('cores', [2, 5]);`.

### Posterior inference

#### Parameter posterior

The FICOS calibration script above will automatically sample from the posterior distribution for the simulator parameters. The samples will be saved in the `calibration_result.mat` file in the variables `ws_rt`. The posterior can be visualized with the function `plotPostJointMargDens`. When called without arguments, this will plot the example results in `results/calibration_example`. To plot different results, call this function with the path to the `calibration_results.mat` produced by the calibration script:

`plotPostJointMargDens('<path_to_results_folder>/calibration_results.mat')`.

The figure plotted with the example results corresponds to Figure 3 in the [manuscript](https://doi.org/10.1016/j.ecoinf.2024.102976).

#### Posterior predictive simulations

You can also use the For posterior sampling, start a Matlab session in the project rood folder and run `postPredSimulations`. By default this will start posterior sampling from an example result file. You may sample instead from your own results by calling this function with the named `sampleFile` argument, e.g. 

`postPredSimulations('sampleFile','<path_to_results_folder>/calibration_results.mat');`. 

#### Posterior predictive summaries

After all posterior predictive simulations have finished, their results can be collected and summarised with 

`postPredTimeseries(postPredPath);`

where `postPredPath` is the root folder for the posterior predictive simulation results.

This function will also plot a figure for the summary. Once you have summarised results with `postPredTimeseries` and saved the summary to a file, such as `results/calibration_example/postPredSummary.mat`, you can produce the plots again with

`plotPostPredTimeseries('results/calibration_example/postPredSummary.mat');`.

This figure plotted with the above function corresponds to Figure 4 in the [manuscript](https://doi.org/10.1016/j.ecoinf.2024.102976).

The example result file `results/calibration_example/postPredSummary.mat` is included as a release asset for [v1.0.2RC](https://github.com/kkaurila/ficosUQ/releases/tag/v1.0.2RC) and can be downloaded with the script `src/scripts/get_example_results.sh`.

#### Load reduction scenarios

Load reduction scenarios, where catchment area nutrients loads are reduced by 20, 40, 60 and 80 percent can be launched with 

`runCatchmentScenarios`.

With default arguments this will use posterior samples from the example results, to use other calibration results or sampes, see `help runCatchmentScenarios`.

Once all the simulations for the reduction scenarios have finished, the results can be collected and summarised with 

`tbChlaSum = summariseChlaScenarios(resultDir)`,

where `resultDir` is the folder containing predictions for all of the load reductions scenarios from `runCatchmentScenarios` . This folder should contain sub folders `scenario_1/`, `scenario_2/`, etc. which themselves contain predictions for each sample.

Finally, the resulting summary table `tbChlaSum` can be plotted with `plotScenarioChlaDist(tbChlaSum)`.

An example summary table is provided as a [release asset](https://github.com/kkaurila/ficosUQ/releases/download/v1.0.3/example_scenarioPredictionTable.mat).

## Project structure

### Outline

This project folder is structured as follows:

* `ficosUQ/` project root
    * `README.md` *this file*
    * `install.sh` **script for installing all dependencies** 
    * `bin/` Compiled binaries
        * `ficos/` FICOS simulator
            * `README.md` separate README for the FICOS simulator, **includes license information for FICOS**
            * the source code and compiled binaries for FICOS
        * `hdf5/` external hdf5 binaries needed to run FICOS
    * `data/`
        * calibration data and .hdf5 input files for the simulator, see *Data* section below
    * `src/` Source code for this project
        * `ficos_calibration.m` **Matlab function for launching the main calibration script**
        * `python/` python source code used in this project
            * `python_install.sh` script for installing python dependencies in a virtual environment
            * `download_data.sh` script for downloading required data using Python library [zenodo_get](https://github.com/dvolgyes/zenodo_get)
            * `changeloadsnutrients.py` Python script for changing nutrient loads.
            * `solar_ma.py` Python script for smoothing boundary condition solar radiation.
            * `pyEnv/` Python virtual environment will be installed here 
        * `submodules/` git sub modules for external Matlab packages
            * `gpstuff/` [GitHub repository](https://github.com/gpstuff-dev/gpstuff/tree/114937ec0a201306489a66cbba38283e722fb998)
            * `BrewerMap/`  [GitHub repository](https://github.com/DrosteEffect/BrewerMap/tree/1773eb111abe7db7a4197e09a53ca42e223f28eb)
        * Other folders containing Matlab source code for this project, organized by topic (`calibration`, `simulator`, etc.).
    

## Dependencies

### Operating System

The FICOS calibration code has been developed and tested under `Ubuntu 20.04` and `Ubuntu 22.04` (both 64 bit). We have not tested other operating systems, but assume the code will also work on similar versions of other Linux distributions. 

Mac OS and Windows are not supported. The majority of the MATLAB code in this project (excluding *external* code in sub modules) may in principle work with these operating systems as well, apart from a few functions using common Linux programs, such as `sed`.

### Compilers

#### gcc and gfortran

The FICOS simulator and the hdf5 libraries have been compiled and tested using versions **9.4.0** and **9.5.0** for both `gcc` and `gfortran`. The script for compiling FICOS assumes compatible versions of these compilers to be executable with commands `gcc-9` and `gfortran-9`. On Ubuntu, both compilers these can be installed with

`sudo apt-get install gfortran-9`

For other Linux distributions, you may need to replace `apt-get` above with another package manager (e.g. `yum`, `dnf`).

You can check which version of these compilers you have installed with 

`gfortran-9 --version` and

`gcc-9 --version`.

If you have installed the above compilers in some other manner, you may need to change the `CC` and `FC` options in the line starting with `./configure` in `bin/install_hdf5_libraries.sh`. The compiler configuration may be updated to be more flexible in a future release.

### Matlab

FICOS calibration and posterior predictions have been run on a server with `MATLAB R2022b`, while the posterior summaries and figures were done on a laptop with `MATLAB R2024a`.

The main calibration script requires the following Matlab toolboxes:
 - Control System Toolbox
 - Optimization Toolbox
 - Statistics and Machine Learning Toolbox
 - Curve Fitting Toolbox
 - Parallel Computing Toolbox

Toolbox dependencies were checked with 
`matlab.codetools.requiredFilesAndProducts('src/ficos_calibration.m')`

Most of the Matlab toolboxes above are required by the `gpstuff` package used for the Gaussian Process emulator. Parallel Computing Toolbox is needed for running multiple simulations in parallel.

See also `src/util/findDependencies.m`.

### External libraries

#### Installation script

All of the external libraries below can be installed by running the install.sh script in the project root directory. The script assumes you have installed the compilers `gcc-9` and `gfortran-9` and Matlab Toolboxes above.

#### hdf5 1.8.19

The hdf5 library version `1.8.19` and its associated Fortran compiler script is installed via the installation script `bin/install_hdf5_libraries.sh`.
This requires the `gcc-9` and `gfortran-9` compilers mentioned above.

#### Python 3.8.10

The python functions in src/python have been tested with Python 3.8.10 with using the following libraries:
- `numpy 1.17.4`
- `h5py  2.10.0`

Python 3.8.10 as well as the above libraries are installed into a virtual environment using the script `src/python/python_install.sh`.

## Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14028484.svg)](https://doi.org/10.5281/zenodo.14028484) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14028433.svg)](https://doi.org/10.5281/zenodo.14028433)

The data used for calibrating FICOS is stored in two Zenodo repositories: [10.5281/zenodo.14028484](https://doi.org/10.5281/zenodo.14028484) for the calibration data and [10.5281/zenodo.14028433](https://doi.org/10.5281/zenodo.14028433) for the .hdf5 input files. All of the following files need to be placed in the `data/` folder:

- `intensiveStationsData.xlsx`
- `hd_files.hdf5`
- `hd_files_ma_5.hdf5`
- `loading.hdf5`
- `loading_Aurajoki.hdf5`

### Automatic download

The main installation script attempts to download all of the above files automatically. This may fail if the Zenodo servers are busy, in which case you may try again later by executing `src/python/download_data.sh` on a terminal in the project root folder.

### Manual download

You may also download the files manually from the following addresses:
- https://doi.org/10.5281/zenodo.14028484
- https://doi.org/10.5281/zenodo.14028433
