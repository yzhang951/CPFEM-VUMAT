# Crystal Plasticity Finite Elmenet (CPFE) model for AM-HEA
Ubuntu operating system 18.04
Abaqus 6.13

## Input code
aeuler, slipsys, J_8000.inp, vumat_dual_phase.f, run_comp_8000.sh

## File usage
vumat_dual_phase.f		: main VUMAT file

aeuler					: orientation file, euler angles in deg.
						  find more in Appendix A of https://doi.org/10.1016/0022-5096(92)80003-9

slipsys					: slip system infomation

J_8000.inp				: example input file for Abaqus

run_comp_8000			: shell script to run the whole simulation

loading_curve.py		: post-processing python script to extract stress-strain curve


## Important things in VUMAT file
NOEL in line 191		: number of elements
FILE 1~5 in line		: path of aeuler file, slipsys file and other output


## aeuler file
The euler angles information start in line 5. 
Column 1-3				: euler angles in degree
Column 4-6				: eigenstrain to simulate type-II internal stress


## J_8000.inp
We use mass scaling to accerlate the simulations.
Lower the density if your have a convergence issue.

## Post-processing script
loading_curve.py, fit.sh, plot.plg