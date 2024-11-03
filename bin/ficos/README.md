# About

This folder contains a pre-release of the FICOS water quality simulator developed by the Finnish Environment Institute (Syke), see the License section.

# Requirements

## hdf5-1.8.19

Before installing the FICOS simulator, make sure you have installed the hdf5-libraries in ficosUQ/bin/hdf5/hdf5-1.8.19/ . These libraries need to be compiled with a specific configuration and may require using older gcc and gfortran compilers.
The required hdf5 libraries and FICOS can both be installed with the script ficosUQ/install_hdf5_libraries.sh or with the general installation script ficosUQ/install.sh. Note that all installation scripts need to be executed from the 'ficosUQ' folder.

## hdf5 input files

When using FICOS, you also need two hdf5 input files named hd_files.hdf5 and loading.hdf5. The former contains the boundary conditions of the simulation, while the latter contains nutrient loadings. The calibration package assumes these files to be located in the folder 'ficosUQ/data'.

## System

This version of FICOS has been tested on Linux, Ubuntu 22.04 specifically. All instructions and installation scripts assume they are used on Linux. Although it may be possible to use this on other operating systems, such as Windows or Mac, with some additional steps, these are not supported. 

# Installing FICOS

The preferrable method for installing FICOS and all other components of the calibration package is to execute the installation script ficosUQ/install.sh. FICOS may also be installed independently, as long as the requirements above are met with the following steps:

1. In a terminal, navigate to the 'ficosUQ' folder
2. execute the following command: 'bin/ficos/make'

# Using FICOS

This version of FICOS is meant to be used with the Matlab package with functions 'run_sim' and 'runFICOSparallel'.
Assuming you are in the 'ficosUQ' folder, FICOS can also be used independently of the calibration code with the following command:

bin/ficos/wqficos bin/ficos/ficossettings.ini [resultdir] data/hd_files.hdf5 data/loading.hdf5 bin/ficos/wqficos.ini test 0 0

This will launch the simulation using the default parametrisation and the outputs will be saved in [resultdir]/results.hdf5.

**Do not alter any of the files in the command above.** The ficossettings.ini file is used as the template for other parametrisations. To try a different parametrisation, use a copy of this file instead. The Matlab package provides functions for constructing input files for alternate parametrisations and scenarios. 

# License

This pre-release of the FICOS water quality code is provided solely for the purposes of the publication of Kaurila et al. with the permission of FICOS developers and does not in any way limit or bind the authors of the FICOS code in regard to future official publication or licensing of the code. The FICOS code provided here can be freely used with the proper attribution (Syke/FICOS team).

The FICOS code will be officially released later (tentatively in 2025) by the FICOS team with proper documentation, licensing and attributions after some clean-up of the code. The core equations will not be affected by this process.




