## CPFEM-VUMAT
Crystal plasticity finite element code, VUMAT files for Abaqus

# Introduction
The repository contains Abaqus crystal plasticity VUMAT code for the following papers, 
1. 316steel_NC/		Chen, Wen, et al. "Microscale residual stresses in additively manufactured stainless steel." Nature communications 10.1 (2019): 1-12.

2. AM-HEA/ Ren, Jie, et al. "Strong yet ductile nanolamellar high-entropy alloys by additive manufacturing." under review.

3. Crack/ Baolin, Wang, et al. in preparation.

# Code Structure
The VUMAT codes are written in Fortran, other postprocessing scripts are using Python3 and MATLAB. The main crystal plasticity code is vumat_*.f, this VUMAT code requires external input file such as slip system defination and crystal orientation. Please remember to change the path of these files (variables of FILE1, FILE2 and etc.) when using this code. The sample runing script can be found in the 'run.sh'

After the calculation, the postprocessing script could be used to get the lattice strain response, such as the commands in 'fit.sh'. The outputs will be written into 'fitting.dat', 'lattice_strain_LD.dat' and 'lattice_strain_TD.dat'.

# Contact
Contact Yin Zhang, yzhang951@gatech.edu for more technical details.