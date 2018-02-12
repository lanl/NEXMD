#!/bin/bash

dipole_type=$1
nstates=$2

if [ $dipole_type = 0 ]
   then
       while read line cstep; do
	   grep -A$line -A $nstates -m 1 'Ground State Molecular Dipole Moment (A.U.)' md.out >> gsdipole.out
       done < dipline.out
fi

if [ $dipole_type = 1 ]
   then
       while read line cstep; do
	   grep -A$line -A $nstates -m 1 'Frequencies (eV) and Transition Dipole Moments (AU)' md.out >> transdipole.out
       done < dipline.out
fi

if [ $dipole_type = 2 ]
   then
       while read line cstep; do
	   grep -A$line -A $nstates -m 1 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' md.out >> excdipole.out
       done < dipline.out
fi
