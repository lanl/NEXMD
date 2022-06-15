#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares input files for non-adiabatic              |
| excited-state molecular dynamics (NEXMD) after a push pulse.      |
| This function is used when simulating pump-push-probe             |
| spectroscopy.                                                     |
|                                                                   |
| A general header called 'header' must be in the NEXMD directory   |
| (e.g. NEXMD) and must have all inputs set except for:             |
|                                                                   |
| 1) Random seed (rnd_seed)                                         |
| 2) Initial excited state (exc_state_init_flag)                    |
| 3) Initial nuclear coordinates and velocities (nucl_coord_veloc)  |
| 4) Initial quantum amplitudes and phase (quant_amp_phase)         |
|                                                                   |
| In parentheses shown above, are flags in 'header' that label      |
| these inputs.  This function finds these labels and fills them    |
| in accordingly.                                                   |
|                                                                   |
| NOTE: All NEXMD folders and rseedslists, inside the NEXMD         |
| directory, will be deleted if this function is completely         |
| executed!                                                         |
|___________________________________________________________________|

'''

import numpy as np
import random
import os
import sys
import shutil
import glob
import fileinput

CWD = os.getcwd()

def NEXMD():
    
    print('Preparing input files for NEXMD after push pulse.')
    
## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory after pump pulse: ')
    if not os.path.exists(NEXMDIR):
        print('Path %s does not exist.' % (NEXMDIR))
        sys.exit()
    SPDIR = raw_input('Single-point calculations directory before push pulse: ')
    if not os.path.exists(SPDIR):
        print('Path %s does not exist.' % (SPDIR))
        sys.exit()
    OUTDIR = raw_input('Output directory [e.g. NEXMD_push]: ')
    if not os.path.exists(OUTDIR):
        print('Path %s does not exist.' % (OUTDIR))
        sys.exit()

## DELETE PREVIOUS NEXMD FOLDERS AND RSEEDSLISTS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (OUTDIR))
    NEXMDS.sort()
    if len(NEXMDS) != 0:
        CONTQ = input('** WARNING ** All NEXMD folders and rseedslists inside %s will be deleted!\nDo you want to continue? Answer YES [1] or NO [0]: ' % (OUTDIR))
        if CONTQ not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if CONTQ == 0:
            sys.exit()
    for NEXMD in NEXMDS:
        print('Deleting', '%s' % (NEXMD))
        shutil.rmtree(NEXMD)
    RSEEDSLIST = glob.glob('%s/rseedslist*' % (OUTDIR))
    for RSEEDS in RSEEDSLIST:
        os.remove(RSEEDS)
        print('Deleting', '%s' % (RSEEDS))

## CHECK RUNNING DYNAMICS AND GET NSTATES ##
    if not os.path.exists('%s/header' % (OUTDIR)):
        print('Path %s/header does not exist.' % (OUTDIR))
        sys.exit()
    HEADER = open('%s/header' % (OUTDIR),'r')
    HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1])
    if TSMAX <= 0:
        print('Must change n_class_steps in %s/header to greater than 0 for dynamics.' % (OUTDIR))
        sys.exit()

## CHOOSE GEOMETRIES ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print('There are no NEXMD folders in %s.' % (SPDIR))
        sys.exit()
    with open('%s/totdirlist' % (OUTDIR),'w') as DATA:
        for NEXMD in NEXMDS:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print('Path %sdirlist1 does not exist.' % (NEXMD))
                sys.exit()
            INPUT = fileinput.input('%s/dirlist1' % (NEXMD))
            DATA.writelines(INPUT)
    DIRLIST1 = np.int_(np.genfromtxt('%s/totdirlist' % (OUTDIR)))
    if isinstance(DIRLIST1,int) == True:
        DIRLIST1 = np.array([DIRLIST1])
    os.remove('%s/totdirlist' % (OUTDIR))
    NTRAJ = len(DIRLIST1)
    NTRAJQ = input('How many trajectories for NEXMD? Enter a number no greater than %d: ' % (NTRAJ))
    if isinstance(NTRAJQ, int) == False:
        print('Number of trajectories must be integer.')
        sys.exit()
    if NTRAJQ == 0:
        print('Number of trajectories must be positive integer.')
        sys.exit()
    if np.abs(NTRAJQ) > NTRAJ:
        print('Number of trajectories must be less than or equal to %d.' % (NTRAJ))
        sys.exit()
    if NTRAJQ < 0:
        NTRAJ = np.abs(NTRAJQ)
    else:
        INTERVAL = np.int(np.ceil(NTRAJ/np.float(NTRAJQ)))
        if INTERVAL*NTRAJQ > NTRAJ:
            INTERVAL = INTERVAL - 1
            NTRAJ = len(DIRLIST1[0::INTERVAL])
        else:
            NTRAJ = NTRAJQ
        COORDSQ = input('You have requested %d evenly-spaced coordinate files in the range %d to %d for NEXMD.\nContinue? Answer YES [1] or NO [0]: ' % (NTRAJ,DIRLIST1[0],DIRLIST1[0::INTERVAL][-1]))
        if COORDSQ not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if COORDSQ == 0:
            sys.exit()

## SPLIT GEOMETRIES ##
    SPLIT = input('Number of trajectories per NEXMD folder: ')
    if isinstance(SPLIT, int) == False:
        print('Number of trajectories per NEXMD folder must be integer.')
        sys.exit()
    if SPLIT <= 0:
        print('Number of trajectories per NEXMD folder must be integer greater than zero.')
        sys.exit()
    DIRSPLIT = SPLIT*np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT)+1)
    DIRSPLIT[-1] = DIRLIST1[-1]
    if NTRAJ < 0:
        DIRSPLIT = np.split(np.sort(random.sample(DIRLIST1,NTRAJ)),DIRSPLIT)
    else:
        DIRSPLIT = np.split(DIRLIST1[0::INTERVAL],DIRSPLIT)

## CHOOSE RANDOM SEEDS ##
    RANDQ = input('New random seeds? Answer YES [1] or NO [0]: ')
    if RANDQ not in [1,0]:
        print('Answer must be 1 or 0.')
        sys.exit()
    if RANDQ == 1:
        RSEEDS = random.sample(np.arange(1,1000001), NTRAJ)
    else:
        RSEEDSLIST = raw_input('Path to random-seeds list: ')
        if not os.path.exists(RSEEDSLIST):
            print('Path %s does not exist.' % (RSEEDSLIST))
            sys.exit()
        RSEEDSLIST = open('%s' % (RSEEDSLIST),'r')
        RSEEDSLIST = RSEEDSLIST.readlines()
        LEN = len(RSEEDSLIST)
        if LEN < NTRAJ:
            print('Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories requested is %d.' % (LEN,NTRAJ))
            sys.exit()
        RSEEDS = np.zeros(LEN)
        INDEX = 0
        for LINE in RSEEDSLIST:
            VAL = LINE.split()
            RSEEDS[INDEX] = np.float(VAL[0])
            INDEX += 1
    RSEEDS = np.int_(RSEEDS)

## INDICES FOR EXCITED-TO-EXCITED OSCILLATOR STRENGTHS ##
    INDICES = np.cumsum(np.insert(np.arange(1,NSTATES + 2)[::-1],[0],0,axis=0))

## PREPARE NEXMD INPUT FILES ##
    NDNEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NDNEXMDS.sort()
    SPNEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    SPNEXMDS.sort()
    ERROR = open('%s/NEXMD_push.err' % (CWD),'w')
    EXCEN = input('Laser excitation energy in eV: ')
    EXCSD = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
    TRAJ = 0
    INDEX = 0
    PSHFLAG = 0
    for NEXMD in np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT)+1):
        os.makedirs('%s/NEXMD%d' % (OUTDIR,NEXMD))
        DIRLIST = open('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD),'w')
        for DIR in DIRSPLIT[INDEX]:
            os.makedirs('%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR))
## GET COORDINATES, VELOCITIES, AND EXCITED-STATE BEFORE PUSH PULSE ##
            for NDNEXMD in NDNEXMDS:
                if os.path.exists('%s/%04d/restart.out' % (NDNEXMD,DIR))
                    DATA = open('%s/%04d/restart.out' % (NDNEXMD,DIR),'r')
                    IPSHFLAG = 0
                    break
                else:
                    IPSHFLAG == 1
            if IPSHFLAG == 1:
                ERROR.wirte( '%s%04d/restart.out' % (NDNEXMD,DIR), 'does not exist')
                PSHFLAG = 1
                TRAJ += 1
                continue
            DATA = DATA.readlines()
            CINDEX = 0
            ARRAY = np.array([])
            for LINE in DATA:
                if 'State' in LINE:
                    ISTATE = np.int(LINE.split()[-1])
                if '$COORD' in LINE:
                    ARRAY = np.append(ARRAY,CINDEX)
                if '$ENDCOORD' in LINE:
                    ARRAY = np.append(ARRAY,CINDEX)
                if '$VELOC' in LINE:
                    ARRAY = np.append(ARRAY,CINDEX)
                if '$ENDVELOC' in LINE:
                    ARRAY = np.append(ARRAY,CINDEX)
                CINDEX += 1
            ARRAY = np.int_(ARRAY)
            if len(ARRAY) != 4:
                ERROR.wirte( '%s%04d/restart.out' % (NEXMD,DIR), 'is incomplete')
                PSHFLAG = 1
                TRAJ += 1
                continue
            COORDS = DATA[ARRAY[0]:ARRAY[1]+1:1]
            VELOCS = DATA[ARRAY[2]:ARRAY[3]+1:1]
## GET EXCITED-TO-EXCITED DIPOLE MOMENTS ##
            for SPNEXMD in SPNEXMDS:
                if os.path.exists('%s/%04d/muab.out' % (SPNEXMD,DIR)):
                    DATA = np.genfromtxt('%s/%04d/muab.out' % (SPNEXMD,DIR))[INDICES[ISTATE-1]:INDICES[ISTATE]-1:1]
                    IPSHFLAG = 0
                    break
                else:
                    IPSHFLAG = 1
            if IPSHFLAG == 1:
                ERROR.wirte( '%s%04d/muab.out' % (SPNEXMD,DIR), 'does not exist')
                PSHFLAG = 1
                TRAJ += 1
                continue
## DETERMINE NEW EXCITED-STATE AND QUANTUM COEFFICIENTS ##
            LEN = len(DATA)
            QPOP = np.zeros(LEN)
            OINDEX = 0
            for LOW, HIGH, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                QPOP[OINDEX] = OSCSTREN*np.exp(-(ENERGY - EXCEN)**(2.0)/(2.0*EXCSD**(2.0)))/np.sqrt(2.0*np.pi*EXCSD**(2.0))
                OINDEX += 1
            QPOP = QPOP/np.sum(QPOP)
            STATE = np.searchsorted(np.cumsum(QPOP),np.random.uniform()) + ISTATE
            QPOP = np.zeros(LEN) + ISTATE - 1
            QPOP[STATE-1] = 1.0
## FILL IN INPUT FILE ##
            INPUT = open('%s/NEXMD%d/%04d/input.ceon' % (OUTDIR,NEXMD,DIR),'w')
            for LINE in HEADER:
                if 'rnd_seed' in LINE:
                    INPUT.write('   rnd_seed=%d, ! Seed for the random number generator\n' % (RSEEDS[TRAJ]))
                else:
                    if 'exc_state_init_flag' in LINE:
                        INPUT.write('   exc_state_init=%d, ! Initial excited state (0 - ground state) [0]\n' % (STATE))
                    else:
                        if 'nucl_coord_veloc' in LINE:
                            for LINE in COORDS:
                                INPUT.write(LINE)
                            INPUT.write('\n')
                            for LINE in VELOCS:
                                INPUT.write(LINE)
                        else:
                            if 'quant_amp_phase' in LINE:
                                INPUT.write('&coeff\n')
                                for LINE in QPOP:
                                    INPUT.write('  %.3f  %.3f\n' % (LINE,0.0))
                                INPUT.write('&endcoeff\n')
                            else:
                                INPUT.write(LINE)
            DIRLIST.wirte( '%04d' % (DIR))
            print('%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR))
            TRAJ += 1
        DIRLIST.close()
        shutil.copyfile('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD), '%s/NEXMD%d/dirlist1' % (OUTDIR,NEXMD))
        INDEX += 1
    if PSHFLAG == 1:
        print('One or more NEXMD trajectories cannot be prepared, check NEXMD_push.err.')
    np.savetxt('%s/rseedslist' % (OUTDIR), np.transpose(RSEEDS[0:TRAJ:1]))
