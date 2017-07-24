#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares input files for non-adiabatic              |
| excited-state molecular dynamics (NEXMD).                         |
|                                                                   |
| A general header called 'header' must be in the NEXMD directory   |
| (e.g. nexmd) and must have all inputs set except for:             |
|                                                                   |
| 1) Random seed (rnd_seed)                                         |
| 2) Initial excited state (exc_state_init_flag)                    |
| 3) Initial nuclear coordinates and velocities (nucl_coord_veloc)  |
| 4) Initial quantum amplitudes and phase (quant_amp_phase)         |
|                                                                   |
| In parentheses shown above, are flags in 'header' that label      |
| these inputs. This function finds these labels and fills them in  |
| accordingly. A 'ceo.err' file will be generated in case there     |
| exists incomplete single-point calculations.                      |
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
    
    print 'Preparing input files for NEXMD.'
    
## DIRECTORY NAMES ##
    GSDIR = raw_input('Ground-state dynamics directory: ')
    if not os.path.exists(GSDIR):
        print 'Path %s does not exist.' % (GSDIR)
        sys.exit()
    SPDIR = raw_input('Single-point calculations directory: ')
    if not os.path.exists(SPDIR):
        print 'Path %s does not exist.' % (SPDIR)
        sys.exit()
    OUTDIR = raw_input('Output directory [e.g. nexmd]: ')
    if not os.path.exists(OUTDIR):
        print 'Path %s does not exist.' % (OUTDIR)
        sys.exit()

## DELETE PREVIOUS NEXMD FOLDERS AND RSEEDSLISTS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (OUTDIR))
    NEXMDS.sort()
    if len(NEXMDS) != 0:
        CONTQ = input('** WARNING ** All NEXMD folders and rseedslists inside %s will be deleted!\nDo you want to continue? Answer YES [1] or NO [0]: ' % (OUTDIR))
        if CONTQ not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if CONTQ == 0:
            sys.exit()
    for NEXMD in NEXMDS:
        print 'Deleting', '%s' % (NEXMD)
        shutil.rmtree(NEXMD)
    RSEEDSLIST = glob.glob('%s/rseedslist*' % (OUTDIR))
    for RSEEDS in RSEEDSLIST:
        os.remove(RSEEDS)
        print 'Deleting', '%s' % (RSEEDS)

## CHECK RUNNING DYNAMICS ##
    if not os.path.exists('%s/header' % (OUTDIR)):
        print 'Path %s/header does not exist.' % (OUTDIR)
        sys.exit()
    HEADER = open('%s/header' % (OUTDIR),'r')
    HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1])
            break
    if TSMAX <= 0:
        print 'Must change n_class_steps in %s/header to greater than 0 for dynamics.' % (OUTDIR)
        sys.exit()

## FIND GEOMETRIES ##
    if not os.path.exists('%s/coords.xyz' % (GSDIR)):
        print 'Path %s/coords.xyz does not exist.' % (GSDIR)
        sys.exit()
    DATAC = open('%s/coords.xyz' % (GSDIR),'r')
    DATAC = DATAC.readlines()
    LENC = len(DATAC)
    if not os.path.exists('%s/velocity.out' % (GSDIR)):
        print 'Path %s/velocity.out does not exist.' % (GSDIR)
        sys.exit()
    DATAV = open('%s/velocity.out' % (GSDIR),'r')
    DATAV = DATAV.readlines()
    LENV = len(DATAV)
    NCOORDS = 0
    INDEX = 0
    ARRAYC = np.array([])
    for LINE in DATAC:
        if 'time' in LINE:
            if NCOORDS == 0:
                TINIT = np.float(LINE.split()[-1])
            else:
                TIME = np.float(LINE.split()[-1])
            NCOORDS += 1
            ARRAYC = np.append(ARRAYC,INDEX)
        INDEX += 1
    ARRAYC = np.append(ARRAYC,LENC + 1)
    ARRAYC = np.int_(ARRAYC)
    if NCOORDS == 0:
        print 'No coordinates were found.'
        sys.exit()
    if NCOORDS == 1:
        print 'Only initial coordinates, at %.2f fs, were found.' % (TINIT)
        sys.exit()
    if NCOORDS > 1:
        INDEX = 0
        ARRAYV = np.array([0])
        for LINE in DATAV:
            if 'time' in LINE:
                ARRAYV = np.append(ARRAYV,INDEX)
            INDEX += 1
        ARRAYV = np.append(ARRAYV,LENV)
        ARRAYV = np.int_(ARRAYV)
    TINC = TIME/(NCOORDS - 1)
    print 'A total of %d coordinate files, ranging from %.2f to %.2f fs in increments of %.2f fs, were found.\nNote: Only coordinate files used for single-point calculations can be used for NEXMD.' % (NCOORDS,TINIT,TIME,TINC)

## CHOOSE GEOMETRIES ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (SPDIR)
        sys.exit()
    with open('%s/totdirlist' % (OUTDIR),'w') as DATA:
        for NEXMD in NEXMDS:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
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
        print 'Number of trajectories must be integer.'
        sys.exit()
    if NTRAJQ == 0:
        print 'Number of trajectories must be positive integer.'
        sys.exit()
    if np.abs(NTRAJQ) > NTRAJ:
        print 'Number of trajectories must be less than or equal to %d.' % (NTRAJ)
        sys.exit()
    if NTRAJQ < 0:
        NTRAJ = np.abs(NTRAJQ)
        COORDSQ = input('You have requested %d randomly-selected coordinate files in the range %d to %d.\nContinue? Answer YES [1] or NO [0]: ' % (NTRAJ,DIRLIST[0],DIRLIST[-1]))
    else:
        INTERVAL = np.int(np.ceil(NTRAJ/np.float(NTRAJQ)))
        if INTERVAL*NTRAJQ > NTRAJ:
            INTERVAL = INTERVAL - 1
            NTRAJ = len(DIRLIST1[0::INTERVAL])
        else:
            NTRAJ = NTRAJQ
        COORDSQ = input('You have requested %d evenly-spaced coordinate files in the range %d to %d for NEXMD.\nContinue? Answer YES [1] or NO [0]: ' % (NTRAJ,DIRLIST1[0],DIRLIST1[0::INTERVAL][-1]))
    if COORDSQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if COORDSQ == 0:
        sys.exit()

## SPLIT GEOMETRIES ##
    SPLIT = input('Number of trajectories per NEXMD folder: ')
    if isinstance(SPLIT, int) == False:
        print 'Number of trajectories per NEXMD folder must be integer.'
        sys.exit()
    if SPLIT <= 0:
        print 'Number of trajectories per NEXMD folder must be integer greater than zero.'
        sys.exit()
    DIRSPLIT = SPLIT*np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT)+1)
    DIRSPLIT[-1] = DIRLIST1[-1]
    if NTRAJQ < 0:
        DIRSPLIT = np.split(np.sort(random.sample(DIRLIST1,NTRAJ)),DIRSPLIT)
    else:
        DIRSPLIT = np.split(DIRLIST1[0::INTERVAL],DIRSPLIT)

## EXTRACT ATOMIC NUMBERS ##
    if not os.path.exists('%s/restart.out' % (GSDIR)):
        print 'Path %s/restart.out does not exist.' % (GSDIR)
        sys.exit()
    ANUM = open('%s/restart.out' % (GSDIR),'r')
    ANUM = ANUM.readlines()
    TOP = None
    BOTTOM = None
    INDEX = 0
    for LINE in ANUM:
        if '$COORD' in LINE:
            TOP = INDEX
        if '$ENDCOORD' in LINE:
            BOTTOM = INDEX
            break
        INDEX += 1
    if isinstance(TOP, int) == True and isinstance(BOTTOM, int) == True:
        ANUM = [ LINE.split()[0] for LINE in ANUM[TOP+1:BOTTOM:1] ]
    else:
        print 'There is a problem with %s/restart.out.' % (GSDIR)
        sys.exit()

## CHOOSE RANDOM SEEDS ##
    RANDQ = input('New random seeds? Answer YES [1] or NO [0]: ')
    if RANDQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if RANDQ == 1:
        RSEEDS = random.sample(np.arange(1,1000001), NTRAJ)
    else:
        RSEEDS = raw_input('Path to random-seeds list: ')
        if not os.path.exists(RSEEDS):
            print 'Path %s does not exist.' % (RSEEDS)
            sys.exit()
        RSEEDS = np.int_(np.genfromtxt('%s' % (RSEEDS)))
        if isinstance(RSEEDS,int) == True:
            RSEEDS = np.array([RSEEDS])
        LEN = len(RSEEDS)
        if LEN < NTRAJ:
            print 'Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories requested is %d.' % (LEN,NTRAJ)
            sys.exit()
    for RSEED in RSEEDS:
        if RSEED < 0:
            print 'A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (RSEED)
            sys.exit()
    RSEEDS = np.int_(RSEEDS)

## PREPARE NEXMD INPUT FILES ##
    SPNEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    SPNEXMDS.sort()
    ERROR = open('%s/ceo.err' % (CWD),'w')
    STYPE = input('Spectral lineshape? Answer GAUSSIAN [0] or LORENTZIAN [1]: ')
    if STYPE not in [0,1]:
        print 'Answer must be 0 or 1.'
        sys.exit()
    EXCEN = input('Laser excitation energy in eV: ')
    if isinstance(EXCEN, int) == False and isinstance(EXCEN, float) == False:
        print 'Excitation energy must be integer or float.'
        sys.exit()
    if STYPE == 0:
        SPECB = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
    else:
        SPECB = input('Spectral broadening (i.e. Lorentzian FWHM) in eV [e.g. 0.36]: ')
    if isinstance(SPECB, int) == False and isinstance(SPECB, float) == False:
        print 'Spectral broadening must be integer or float.'
        sys.exit()
    TRAJ = 0
    INDEX = 0
    CEOFLAG = 0
    for NEXMD in np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT)+1):
        os.makedirs('%s/NEXMD%d' % (OUTDIR,NEXMD))
        DIRLIST = open('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD),'w')
        for DIR in DIRSPLIT[INDEX]:
            os.makedirs('%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR))
            for SPNEXMD in SPNEXMDS:
                if os.path.exists('%s/%04d/ceo.out' % (SPNEXMD,DIR)):
                    DATA = np.genfromtxt('%s/%04d/ceo.out' % (SPNEXMD,DIR))
                    ICEOFLAG = 0
                    break
                else:
                    ICEOFLAG = 1
            if ICEOFLAG == 1:
                print >> ERROR, '%s%04d/ceo.out' % (SPNEXMD,DIR), 'does not exist'
                RSEEDS[TRAJ] = -123456789
                CEOFLAG = 1
                ICEOFLAG = 0
                TTRAJ += 1
                continue
            LEN = len(DATA)
            QPOP = np.zeros(LEN)
            OINDEX = 0
            if STYPE == 0:
                for STATE, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    QPOP[OINDEX] = OSCSTREN*np.exp(-(ENERGY - EXCEN)**(2.0)/(2.0*SPECB**(2.0)))/np.sqrt(2.0*np.pi*SPECB**(2.0))
                    OINDEX += 1
            else:
                for STATE, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    QPOP[OINDEX] = OSCSTREN/(1.0 + ((ENERGY - EXCEN)**2.0)/(SPECB/2.0)**(2.0))/(SPECB*np.pi/2.0)
                    OINDEX += 1
            QPOP = QPOP/np.sum(QPOP)
            STATE = np.searchsorted(np.cumsum(QPOP),np.random.uniform()) + 1
            QPOP = np.zeros(LEN)
            QPOP[STATE-1] = 1.0
            COORDS = DATAC[ARRAYC[DIR]+1:ARRAYC[DIR+1]-1:1]
            VELOCS = DATAV[ARRAYV[DIR]+2:ARRAYV[DIR+1]-1:1]
            INPUT = open('%s/NEXMD%d/%04d/input.ceon' % (OUTDIR,NEXMD,DIR),'w')
            for LINE in HEADER:
                if 'rnd_seed' in LINE:
                    INPUT.write('   rnd_seed=%d, ! Seed for the random number generator\n' % (RSEEDS[TRAJ]))
                else:
                    if 'exc_state_init_flag' in LINE:
                        INPUT.write('   exc_state_init=%d, ! Initial excited state (0 - ground state) [0]\n' % (STATE))
                    else:
                        if 'nucl_coord_veloc' in LINE:
                            INPUT.write('&coord\n')
                            AINDEX = 0
                            for LINE in COORDS:
                                VAL = LINE.split()
                                INPUT.write('{:>6}  {:>12}  {:>12}  {:>12}'.format(ANUM[AINDEX],VAL[1],VAL[2],VAL[3]))
                                INPUT.write('\n')
                                AINDEX += 1
                            INPUT.write('&endcoord\n\n&veloc\n')
                            for LINE in VELOCS:
                                INPUT.write(LINE)
                            INPUT.write('&endveloc\n')
                        else:
                            if 'quant_amp_phase' in LINE:
                                INPUT.write('&coeff\n')
                                for LINE in QPOP:
                                    INPUT.write('  %.3f  %.3f\n' % (LINE,0.0))
                                INPUT.write('&endcoeff\n')
                            else:
                                INPUT.write(LINE)
            print >> DIRLIST, '%04d' % (DIR)
            print '%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR)
            TRAJ += 1
        DIRLIST.close()
        shutil.copyfile('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD), '%s/NEXMD%d/dirlist1' % (OUTDIR,NEXMD))
        INDEX += 1
    np.savetxt('%s/rseedslist' % (OUTDIR), np.transpose(RSEEDS[0:TRAJ:1]))
    if CEOFLAG == 1:
        print 'One or more NEXMD trajectories cannot be prepared, check ceo.err.'
        sys.exit()
    else:
        os.remove('%s/ceo.err' % (CWD))
