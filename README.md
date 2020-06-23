# Spatiotemporal-patterns-of-rodent-hippocampal-field-potentials-uncover-spatial-representations
Source data and code for manuscript "Spatiotemporal patterns of rodent hippocampal field potentials uncover spatial representations".

# Dependencies
1. For MATLAB code, all required functions are included in "Utilities" folder. Just include this folder into you MALTAB path. All MATLAB code were tested on MATLAB 2020a on a Windows desktop.

2. For real-time online decoding c++ code, it requires an c++ complier that support c++11 standard, an openmp openmp library, Intel® Math Kernel Library and Intel® Integrated Performance Primitives. C++ code were tested on a Windows desktop, with Visual Studio 2019 and Intel® Parallel Studio XE environment.

3. For python code, it requires pyhsmm_spiketrains package (packages and install instructions are at https://github.com/lcaopcn/pyhsmm_packages). These packages are worked on linux/Mac OS only, Windows machine may experience some problems. You can use wsl if you only have a Windows compuater. 

# How to run
In each floder, run Figrue*_code.m to visualize the results. For example, run Figure1bc_code.m in "Figure1_source_data", you will get the same figure 1b and 1c as in our manuscript.

# How to cite
Our manuscript has not been published yet. But preprint version of this manuscript is avaliable in BioRxiv:https://doi.org/10.1101/828467.



