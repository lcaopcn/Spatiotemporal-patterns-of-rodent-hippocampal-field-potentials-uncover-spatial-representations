# Spatiotemporal-patterns-of-rodent-hippocampal-field-potentials-uncover-spatial-representations
This repository contains the source data and code for the manuscript "Spatiotemporal patterns of rodent hippocampal field potentials uncover spatial representations".

# Dependencies
1. For the MATLAB code, all required functions are included in the “Utilities” folder. Simply add this folder to your MATLAB path. The MATLAB code has been tested on MATLAB 2020a on a Windows desktop.
   
2. The real-time online decoding C++ code requires a C++ compiler that supports the C++11 standard, an OpenMP library, Intel® Math Kernel Library, and Intel® Integrated Performance Primitives. The C++ code has been tested on a Windows desktop using the Visual Studio 2019 and Intel® Parallel Studio XE environment.z
   
3. The Python code requires the pyhsmm_spiketrains package (the packages and installation instructions can be found at https://github.com/lcaopcn/pyhsmm_packages). These packages are compatible with Linux/Mac OS only. Windows machines may encounter some issues. You can use the Windows Subsystem for Linux (WSL) if you only have a Windows computer.

# How to run
In each folder, run “Figure*_code.m” to visualize the results. For instance, by running Figure1bc_code.m in the “Figure1_source_data” folder, you will reproduce Figures 1b and 1c from our manuscript.

# How to cite
Our manuscript has not yet been published. However, a preprint version of this manuscript can be found on BioRxiv: https://doi.org/10.1101/828467.
