#!/bin/bash
rm	lattice_strain_LD.dat
rm	lattice_strain_TD.dat
/opt/Commands/abq6131 python loading_curve.py
#./plot.plg
#/opt/Commands/abq6131 python lattice_strain.py
