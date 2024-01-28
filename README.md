
# Non-adiabatic EXcited-state Molecular Dynamics (NEXMD)

`nexmd.exe` is a command-line Fortran90 program designed to understand the non-radiative relaxation of
a main-block molecule with 100s of atoms and 10s of relevant excited states
on the picosecond time scale after a vertical photoexcitation.
With this design goal, `nexmd.exe` has implementations of various
non-adiabatic molecular dynamics(NAMD) algorithms,
with efficiency and stability improvements to reach these time scales. 
To provide the electronic structure parameters for these NAMD algorithms,
`nexmd.exe` also performs the self-consistent electronic structure, linear response,
quadratic response and 'Z-vector' calculations 'on-the-fly' ;
using the CEO (collective electronic oscillator) package with semiempirical Hamiltonians
from the SQM package. 
For the classical nuclei and Tully density matrix propagation,
only the classical nuclear geometry, nuclear velocity and the Tully
density matrix at $t=0$ (when the vertical excitation occurs)  need be provided.
NEXMD's modular design allows for individual calculations to be carried out
with different control words in the input file, thus NEXMD is capable of calculating,
among many possibilities, energy minimized geometries and absorption spectra,
at the same level of theory as the NAMD.

------

NEXMD is a scientific program that assumes that the input it is given
is worth running through its algorithms to give reproducible output.

------

There is a helper Python3 script `getexcited.py` to aid in the preparation of the input file, submission
of a swarm of trajectories (parallelization) and processing the ouput (statistics and visualization),
with a copy in this repository and a more frequently updated pypi package<create hyperlink>.  

------

## Highlighted Features

- Solvation correction to the total energy for a single point nuclear geometry.
- Sampling nuclear geometries of the ground-state potential energy surface (PES) under NVE or NVT (through Langevin dynamics).
- Vertical excitation stick and broadened spectrum for a single point nuclear geometry or an ensemble of nuclear geometries.
- Gaussian cube files for the transition density matrix and natural transition orbitals for a vertical excitation. (Relies on external scripts as of now. Integrate it?)
- Minimum energy nuclear geometry and vibrational force constants on any calculated ground or excited state PES.
- Born-Oppenheimer molecular dynamics trajectories on any calculated ground or excited state PES.
- Non-adiabatic molecular dynamics trajectories following the algorithms of trajectory surface hopping (TSH), Ehrenfest, or ab-initio multiple cloning
  (AIMC).
- Geometry constraints for all molecular dynamics simulations. Currently supports freezing bond distances or normal modes of the molecules.

Please refer to the NEXMD manual for a complete list of supported calculations.

## How to cite

When using NEXMD results in academic work, please cite the package:

> Walter Malone, Benjamin Nebgen, Alexander White, Yu Zhang, Huajing Song, 
Josiah A. Bjorgaard, Andrew E. Sifain, Beatriz Rodriguez-Hernandez, 
Victor M. Freixas, Sebastian Fernandez-Alberti, Adrian E. Roitberg, 
Tammie R. Nelson, and Sergei Tretiak,
NEXMD Software Package for Nonadiabatic Excited State 
Molecular Dynamics Simulations,
*Journal of Chemical Theory and Computation* **2020** *16* (9), 5771-5783,
DOI: 10.1021/acs.jctc.0c00248 

> Victor M. Freixas, Walter Malone, Xinyang Li, Huajing Song, Hassiel Negrin-Yuvero, Royle Pérez-Castillo, Alexander White,
Tammie R. Gibson, Dmitry V. Makhov, Dmitrii V. Shalashilin, Yu Zhang, Nikita Fedik, Maksim Kulichenko, Richard Messerly, Luke Nambi Mohanam,
Sahar Sharifzadeh, Adolfo Bastida, Shaul Mukamel, Sebastian Fernandez-Alberti, and Sergei Tretiak,
NEXMD v2.0 Software Package for Nonadiabatic Excited State Molecular Dynamics Simulations,
*Journal of Chemical Theory and Computation* **2023** *19* (16), 5356-5368,
https://doi.org/10.1021/acs.jctc.3c00583

as well as the 2020 review of the methods implemented in NEXMD:

> Tammie R. Nelson, Alexander J. White, Josiah A. Bjorgaard, 
Andrew E. Sifain, Yu Zhang, Benjamin Nebgen, Sebastian Fernandez-Alberti, 
Dmitry Mozyrsky, Adrian E. Roitberg, and Sergei Tretiak,
Non-adiabatic Excited-State Molecular Dynamics: Theory and Applications 
for Modeling Photophysics in Extended Molecular Materials,
*Chemical Reviews* **2020** *120* (4), 2215-2287,
DOI: 10.1021/acs.chemrev.9b00447 

we encourage users to cite the references and libraries used by
NEXMD, which are detailed in the NEXMD manual.

## Usage

### Prerequisites

The following dependencies must be installed and configured before
compiling `nexmd.exe`: (1) a build automation tool such as `make`, (2) a 
FORTRAN90 compiler, (3) and a BLAS/LAPack library compatible with the
compiler and associated libraries. 

1. Build automation tool

`make` is a standard tool

2. FORTRAN90 compiler

- It is recommended to compile the package with `ifort`. If the Intel compiler is
not installed already, check [this
link](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.**8srmr7**)
out.

3. BLAS/LAPACK libraries 

- [Intel
  MKL](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html)

Or equivalently

- [BLAS](http://www.netlib.org/blas/)
- [LAPACK](http://www.netlib.org/lapack)


`getexcited.py` dependencies are listed in its README.

### Set up

- Download the [stable release](https://github.com/lanl/nexmd/releases/latest) or clone the
  developmental version via `git clone https://github.com/lanl/nexmd.git`.
- Go to the root directory of the repo. For the stable version, unzip the
  compressed file first. If cloning, ensure you have checked out the desired branch
  with `git checkout <branch name>`
- Run `make`, which compiles the code with the compiler and libraries specified by options
  in the command line. These options can be checked and modified by
  editing the [Makefile](./Makefile). Custom paths to libraries can be added to line 22
  of the [Makefile](./Makefile) for all options; 
  for further customization of other compilation flags 
  and compilers, it is advisable to fill out the option labeled custom at line xxx. 
  The default target is `ic_mkl`; running `make` is 
  equivalent to `make ic_mkl`, which will compile NEXMD with the
  `ifort` compiler and `mkl` BLAS/LAPack libraries in their default location.
  
- **Optional:** Add the resulting executable (default name `nexmd.exe`) from
  to your PATH. A quick way to do so in the root directory of the repo
  would be `export PATH="$PWD:$PATH"`

### Tests

Example input files and their expected output are in `./tests`, more information on
the tests can be found [here](./tests/README.md). Users are encouraged to
verify the executable runs correctly by comparing the output of their
executable with the expected outputs provided.

### Run

It is helpful to allow `nexmd.exe` to access all of the available stack:

```shell
ulimit -s unlimited
```

If the default input filename `input.ceon` is used, go to the directory
contains the input file and run

```shell
nexmd.exe > [output file]
```

Alternatively, you can run the program with an arbitrary input filename with

```shell
nexmd.exe <[input file] > [output file]
```

`getexcited.py` is typically used to call multiple `nexmd.exe` executions in parallel.
More information on nexmd.exe can be found in the
[manual](./manual/documentation.pdf), more information on getexcited.py can be found
in its repo.

## Bug reporting

If you find a bug in the code, feel free to [open a new
issue](https://github.com/lanl/NEXMD/issues/new) or send an email to <nexmd-users@lanl.gov>.

## Architecture

The program has several modules/sections which are based on the SQM, Amber, and
CEO packages.

- qmmm
- davidson
- others

## Authors

Walter Malone, Victor M. Freixas, Xinyang Li, Hassiel Negrin-Yuvero,
Royle Pérez-Castillo, Dmitry V. Makhov, Dmitrii V. Shalashilin,
Nikita Fedik, Maksim Kulichenko, Richard Messerly, Luke Nambi Mohanam,
Sahar Sharifzadeh, Adolfo Bastida, Shaul Mukamel, Benjamin Nebgen,
Alexander White, Yu Zhang, Huajing Song, Josiah Bjorgaard,
Andrew Sifain, Beatriz Rodriguez-Hernandez, Sebastian Fernandez-Alberti,
Adrian E. Roitberg, Tammie Nelson, Sergei Tretiak

## Acknowledgments

- Los Alamos National Lab (LANL), Center for Nonlinear Studies (CNLS), Center
  for Integrated Nanotechnologies (CINT)
- CONICET
- UNQ
- ANPCyT
  
## Copyright Notice

© 2020. Triad National Security, LLC. All rights reserved. This program was
produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC
for the U.S. Department of Energy/National Nuclear Security Administration. All
rights in the program are reserved by Triad National Security, LLC, and the
U.S. Department of Energy/National Nuclear Security Administration. The
Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to
reproduce, prepare derivative works, distribute copies to the public, perform
publicly and display publicly, and to permit others to do so. Triad National
Security, LLC as management and operations contractor for Los Alamos National
Laboratory, plans to release the NexMD Software under an open source license at
<https://github.com/lanl/NEXMD>. This software was co-authored with several
individuals that include Sebastian Fernandez Alberti as professor at
Universidad Nacional de Quilmes (UNQ) and Adrian Roitberg as professor at the
University of Florida (UF). Triad acknowledges UNQ and UF’s role in
co-authorship of the software.

## License

This program is open source under the BSD-3 License. Redistribution and use in
source and binary forms, with or without modification, are permitted provided
that the following conditions are met:  

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
1. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
1. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
