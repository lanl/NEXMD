#/usr/bin/python

'''
 _________________________________________________________
|                                                         |
| The path below must be in string format, meaning        |
| pathtopath = '<path>', where <path> is in single        |
| quotations.                                             |
|_________________________________________________________|

'''

pathtopack = 'Insert path to getexcited_package here, do not include "getexcited_package" in path'

'''
         _______________________________________
        |#######################################|
        |##|                                 |##|
        |##|        USE WITH CAUTION!        |##|
        |##|                                 |##|
        |##| Any questions or suggestions    |##|
        |##| regarding this script, feel     |##|
        |##| free to contact sifain@usc.edu. |##|
        |##|_________________________________|##|
        |#######################################|
 _________________________________________________________
|                                                         |
| Currently, there are 13 main functions:                 |
|                                                         |
| (1) prepare inputs for single-point calculations        |
| (2) generate a combined optical spectrum from           |
| single-point calculations                               |
| (3) prepare inputs for non-adiabatic excited-state      |
| molecular dynamics (NEXMD)                              |
| (4) prepare input files for an adiabatic simulation,    |
| with geometries taken from NEXMD                        |
| (5) collect populations as a function of time from hops |
| and  quantum coefficients                               |
| (6) collect PESs and NACTs from NEXMD                   |
| (7) prepare restart input files for NEXMD               |
| (8) clean the directories of unfinished trajectories    |
| (9) access options for geometry analysis                |
| (10) access options for dipole analysis                 |
| (11) access options for transition density analysis     |
| (12) access options for pump-push-probe spectroscopy    |
| (13) access code testing tools                          |
|                                                         |
| NOTE: Ground-state trajectory must be completed first.  |
| Coordinates and velocities are selected from a single   |
| ground-state trajectory. The files that contain         |
| coordinates and velocities are coords.xyz and           |
| velocity.out.                                           |
|                                                         |
| NOTE: Main output file from NEXMD program must be       |
| called 'md.out'.  Currently, this is how it is defined  |
| in the main submission script.  If this output file is  |
| renamed, 'md.out' in the 'collectceo.sh' script must    |
| be changed accordingly in order to generate an optical  |
| spectrum from single-point calculations.  Also,         |
| 'md.out' in 'optspec.py' must be changed.               |
|                                                         |
| NOTE: This script can generate a new set of random      |
| seeds or use a list provided by the user.  The latter   |
| option may be important for code testing or             |
| benchmarking purposes.  The user-defined list must      |
| strictly be a list of random seeds, with no header or   |
| footer, and the number of random seeds must be equal to |
| or greater than the number of trajectories requested.   |
|_________________________________________________________|

'''

import sys
import os
if not os.path.exists(pathtopack):
    print 'You must provide the path to getexcited_package in getexcited.py (pathtopack).'
    sys.exit()
sys.dont_write_bytecode = True
sys.path.append('%s' % (pathtopack))
from getexcited_package.spcalc import spcalc
from getexcited_package.optspec import optspec
from getexcited_package.nexmd import nexmd
from getexcited_package.population import population
from getexcited_package.pesnact import pesnact
from getexcited_package.restart import restart
from getexcited_package.newsim import newsim
from getexcited_package.cleandir import cleandir
from getexcited_package.angle import angle
from getexcited_package.dihedral import dihedral
from getexcited_package.bondlength import bondlength
from getexcited_package.bla import bla
from getexcited_package.timing import timing
from getexcited_package.permdipole import permdipole
from getexcited_package.tdiagonal import tdiagonal
from getexcited_package.header import header

funq = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations\n[2] Generate an optical spectrum from single-point calculations\n[3] Prepare input files for NEXMD\n[4] Prepare input files for adiabatic dynamics with geometries from NEXMD\n[5] Collect populations from NEXMD\n[6] Collect pess and nacts from NEXMD\n[7] Prepare restart input files for NEXMD\n[8] Clean out the directories of NEXMD trajectories that are incomplete\n[9] Access options for geometry analysis\n[10] Access options for dipole analysis\n[11] Access options for transition density analysis\n[12] Access options for pump-push-probe spectroscopy (*** UNDER DEVELOPMENT, DO NOT USE ***)\n[13] Access code testing tools\n\nEnter the number corresponding to the desired task: ')
if funq not in [1,2,3,4,5,6,7,8,9,10,11,12,13]:
    print 'Answer must be 1 through 13.'
    sys.exit()
if funq == 1:
    spcalc(header)
if funq == 2:
    optspec(pathtopack,header)
if funq == 3:
    nexmd(header)
if funq == 4:
    newsim(header)
if funq == 5:
    population(header)
if funq == 6:
    pesnact(header)
if funq == 7:
    restart(pathtopack,header)
if funq == 8:
    cleandir(header)
if funq == 9:
    advq = input('\nSelect a task from the following list:\n\n[1] Calculate dihedral angle\n[2] Calculate bond lengths\n[3] Calculate bond length alternation\n[4] Calculate angle between two bonds\n\nEnter the number corresponding to the desired task: ')
    if advq not in [1,2,3,4]:
        print 'Answer must be 1 through 4.'
        sys.exit()
    if advq == 1:
        dihedral(header)
    if advq == 2:
        bondlength(header)
    if advq == 3:
        bla(header)
    if advq == 4:
        angle(header)
if funq == 10:
    advq = input('\nSelect a task from the following list:\n\n[1] Collect excited-state permanent dipole moment\n\nEnter the number corresponding to the desired task: ')
    if advq != 1:
        print 'Answer must be 1.'
        sys.exit()
    if advq == 1:
        permdipole(pathtopack,header)
if funq == 11:
    advq = input('\nSelect a task from the following list:\n\n[1] Analyze induced charge from diagonal elements of the transition density matrix\n\nEnter the number corresponding to the desired task: ')
    if advq not in [1]:
        print 'Answer must be 1.'
        sys.exit()
    if advq == 1:
        tdiagonal(header)
if funq == 12:
    sys.exit()
    advq = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations after pump-push delay time\n[2] Generate optical spectrum from single-point calculations after pump-push delay time\n[3] Prepare input files for NEXMD after push pulse\n\nEnter the number corresponding to the desired task: ')
    if advq not in [1,2,3]:
        print 'Answer must be 1 through 3.'
        sys.exit()
    if advq == 1:
        spcalc_push()
    if advq == 2:
        optspec_push()
    if advq == 3:
        nexmd_push()
if funq == 13:
    advq = input('Select a task from the following list:\n\n[1] Collect timing data from trajectories\n\nEnter the number corresponding to the desired task: ')
    if advq != 1:
        print 'Answer must be 1.'
        sys.exit()
    if advq == 1:
        timing(pathtopack)
