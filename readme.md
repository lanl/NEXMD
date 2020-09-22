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

## Copyright Notice

© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
Triad National Security, LLC as management and operations contractor for Los Alamos National Laboratory,
plans to release the NexMD Software under an open source license at https://github.com/lanl/NEXMD. 
This software was co-authored with several individuals that include Sebastian Fernandez Alberti as professor
at Universidad Nacional de Quilmes (UNQ) and Adrian Roitberg as professor at the University of Florida (UF). 
Triad acknowledges UNQ and UF’s role in co-authorship of the software.

## License

This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:  
1.  Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
 
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
 
3.  Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
