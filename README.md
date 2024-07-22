# CPFEM-VUMAT
Crystal plasticity finite element code, VUMAT files for Abaqus

## Introduction
The repository contains Abaqus crystal plasticity VUMAT code for the following papers, 
1. AM-316steel-NC/		Chen, Wen, et al. ["Microscale residual stresses in additively manufactured stainless steel."](https://www.nature.com/articles/s41467-019-12265-8) Nature Communications, 10, 1-13, (2019)

2. AM-HEA/ Ren, Jie, et al. ["Strong yet ductile nanolamellar high-entropy alloys by additive manufacturing."](https://www.nature.com/articles/s41586-022-04914-8) Nature 608, 62-68, (2022)

3. AM-dislocation-cell/Zhang, Yin, et al. ["Modeling of microscale internal stresses in additively manufactured stainless steel."](https://iopscience.iop.org/article/10.1088/1361-651X/ac8698) Modelling and Simulation in Materials Science and Engineering, 30, 074001 (2022).

4. Crack/ Baolin, Wang, Baolin, et al. ["Modeling of crack tip fields and fatigue crack growth in fcc crystals."](https://www.sciencedirect.com/science/article/pii/S0022509624001571) Journal of the Mechanics and Physics of Solids, 188, 105691 (2024).

## Code Structure
The VUMAT codes are written in Fortran, other postprocessing scripts are using Python3 and MATLAB. The main crystal plasticity code is vumat_*.f, this VUMAT code requires external input files such as slip system definition and crystal orientation. Please remember to change the path of these files (variables of FILE1, FILE2 and etc.) when using this code. The sample running script can be found in the 'run.sh'

After the calculation, the postprocessing script could be used to get the lattice strain response, such as the commands in 'fit.sh'. The outputs will be written into 'fitting.dat', 'lattice_strain_LD.dat' and 'lattice_strain_TD.dat'.

## Contact
Please contact Dr. Yin Zhang yinzhang@pku.edu.cn for more technical details.
