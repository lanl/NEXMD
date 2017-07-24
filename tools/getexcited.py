#/usr/bin/python

'''
 _________________________________________________________
|                                                         |
| The path below must be in string format, meaning        |
| PATHTOPATH = '<PATH>', where <PATH> is in single        |
| quotations.                                             |
|_________________________________________________________|

'''

PATHTOPACK = 'INSERT PATH TO getexcited_package HERE, DO NOT INCLUDE "getexcited_package" IN PATH'

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
| There are nine main functions: (1) prepare inputs for   |
| single-point calculations, (2) generate a combined      |
| optical spectrum from single-point calculations, (3)    |
| prepare inputs for non-adiabatic excited-state          |
| molecular dynamics (NEXMD), (4) collect populations as  |
| a function of time from hops and  quantum coefficients, |
| (5) collect PESs and NACTs from NEXMD, (6) prepare      |
| restart input files for NEXMD, (7) clean the            |
| directories of unfinished trajectories, (8) access      |
| options for geometry analysis, and (9) access options   |
| for pump-push-probe spectroscopy.                       |
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
if not os.path.exists(PATHTOPACK):
    print 'You must provide the path to getexcited_package in getexcited.py (Line 3).'
    sys.exit()
sys.dont_write_bytecode = True
sys.path.append('%s' % (PATHTOPACK))
from getexcited_package.spcalc import SPCALC
from getexcited_package.optspec import OPTSPEC
from getexcited_package.nexmd import NEXMD
from getexcited_package.population import POPULATION
from getexcited_package.pesnact import PESNACT
from getexcited_package.restart import RESTART
from getexcited_package.cleandir import CLEANDIR
from getexcited_package.dihedral import DIHEDRAL
from getexcited_package.bondlength import BONDLENGTH

FUNQ = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations\n[2] Generate an optical spectrum from single-point calculations\n[3] Prepare input files for NEXMD\n[4] Collect populations from NEXMD\n[5] Collect PESs and NACTs from NEXMD\n[6] Prepare restart input files for NEXMD\n[7] Clean out the directories of NEXMD trajectories that are incomplete\n[8] Access options for geometry analysis\n[9] Access options for pump-push-probe spectroscopy (*** UNDER DEVELOPMENT, DO NOT USE ***)\n\nEnter the number corresponding to the desired task: ')
if FUNQ not in [1,2,3,4,5,6,7,8,9]:
    print 'Answer must be 1 through 9.'
    sys.exit()
if FUNQ == 1:
    SPCALC()
if FUNQ == 2:
    OPTSPEC(PATHTOPACK)
if FUNQ == 3:
    NEXMD()
if FUNQ == 4:
    POPULATION()
if FUNQ == 5:
    PESNACT()
if FUNQ == 6:
    RESTART(PATHTOPACK)
if FUNQ == 7:
    CLEANDIR()
if FUNQ == 8:
    ADVQ = input('Select a task from the following list:\n\n[1] Calculate a dihedral angle\n[2] Calculate a bond length\n\nEnter the number corresponding to the desired task: ')
    if ADVQ not in [1,2]:
        print 'Answer must be 1 or 2.'
        sys.exit()
    if ADVQ == 1:
        DIHEDRAL()
    if ADVQ == 2:
        BONDLENGTH()
if FUNQ == 9:
    sys.exit()
    ADVQ = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations after pump-push delay time\n[2] Generate optical spectrum from single-point calculations after pump-push delay time\n[3] Prepare input files for NEXMD after push pulse\n\nEnter the number corresponding to the desired task: ')
    if ADVQ not in [1,2,3]:
        print 'Answer must be 1 through 3.'
        sys.exit()
    if ADVQ == 1:
        SPCALC_PUSH()
    if ADVQ == 2:
        OPTSPEC_PUSH()
    if ADVQ == 3:
        NEXMD_PUSH()
