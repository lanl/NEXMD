# NEXMD (Nonadiabatic EXcited-state Molecular Dynamics
NEXMD is a program for excited-state molecular dynamics. It includes efficient algorithms for nonadiabatic dynamics of molecules in dielectric environments. It is written in Fortran 90, with scripts in Python 2.7 for preparing input files and running the program on parallel systems.

<hr/>

## Useage

The program is run by ```./nexmd``` in a working directory which includes the input file ```input.ceon```. Input files for multiple trajectories can be prepared using the getexcited.py script.

### Prerequisites

The following packages must be installed and configured locally:
* BLAS [http://www.netlib.org/blas/]
* LAPACK [http://www.netlib.org/lapack]

### Set-up

Run ```Make Install``` to compile the program. 

### Run

Run ```nexmd.exe > [output file]``` in the directory with input.ceon. The code is usually run in trivially parallel form with multiple trajectories prepared with the getexcited.py script. See the included manual for more information. (add more)

<hr/>

## Introduction

NEXMD simulates the photoinduced adiabatic and non-adiabatic ground- and excited-state molecular dynamics of organic chromophores. It uses the CEO (collective electronic oscillator) package with a variety of semiempirical methods from the SQM package. Tully’s fewest-switches surface hopping approach to quantum transitions is employed, with instantaneous decoherence and a Min-Cost algorithm for the detection of trivial unavoided crossings. Several TDSCF (time-dependent self-consistent field) QM/continuum models are available for including the effects of a solvent.

## Architecture
The program has several modules/sections which are based on the SQM, Amber, and CEO packages.
* qmmm
* davidson
* others

## Authors

* Walter Malone, Benjamin Nebgen, Alexander White, Yu Zhang, Huajing Song, Josiah Bjorgaard, Andrew Sifain, Beatriz Rodriguez-Hernandez,   Sebastian Fernandez-Alberti, Adrian E. Roitberg, Tammie Nelson, Sergei Tretiak

## Acknowledgments

* Los Alamos National Lab (LANL), Center for Nonlinear Studies (CNLS), Center for Integrated Nanotechnologies (CINT)
* CONICET
* UNQ
* ANPCyT

## License
