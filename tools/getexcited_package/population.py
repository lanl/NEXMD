#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function collects NEXMD populations from surface hopping     |
| and from quantum populations.                                     |
|                                                                   |
| Two files are generated in the current working directory if       |
| collecting NEXMD populations are requested, 'pop.out' and         |
| 'pop.err'.  In 'pop.out', first column is time in fs, followed by |
| the average population on each PES, followed by the sum of all    |
| PES populations (should be = 1.0), followed by the average        |
| populations from quantum coefficients, followed by the sum of     |
| quantum populations (should be = 1.0).  In 'pop.err', first       |
| column is directory of trajectory, third column is the time at    |
| which the trajectory has ended in fs, or 'does not exist'.  The   |
| 'pop.err' file is not generated if all trajectories are complete. |
|                                                                   |
| 'Completed' trajectories are trajectories that have completed     |
| within the time defined by user while executing this function.    |
| 'Excellent' trajectories are trajectories that have completed     |
| within the time defined in 'header' located in the NEXMD          |
| directory.                                                        |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import glob

CWD = os.getcwd()

def POPULATION():

    print 'Collecting populations.'

## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()

## USER-DEFINED LENGTH OF ANALYSIS AND INITIALIZE ARRAYS ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header' % (NEXMDIR),'r')
    for LINE in HEADER:
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1]) + 1
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
    TCOLL = input('Collect populations up to what time in femtoseconds: ')
    if isinstance(TCOLL, int) == False and isinstance(TCOLL, float) == False:
        print 'Time must be integer or float.'
        sys.exit()
    if TCOLL < 0:
        print 'Time must be integer or float greater than zero.'
        sys.exit()
    TCOLL = np.float(TCOLL)
    if TCOLL > (TSMAX - 1)*DT:
        TCOLL = (TSMAX - 1)*DT
    TSCOL = 0
    while TSCOL*DT*ODATA <= TCOLL:
        TSCOL += 1
    FPOPH = np.zeros((TSCOL,NSTATES))
    FPOPC = np.zeros((TSCOL,NSTATES))

## COLLECT POPULATIONS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()
    OUTPUT = open('%s/pop.out' % (CWD),'w')
    ERROR = open('%s/pop.err' % (CWD),'w')
    TTRAJ = 0
    CTRAJ = 0
    ETRAJ = 0
    POPFLAG = 0
    for NEXMD in NEXMDS:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'Path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        for DIR in DIRLIST1:
            if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,DIR)):
                print >> ERROR, '%s%04d/coeff-n.out' % (NEXMD,DIR), 'does not exist'
                POPFLAG = 1
                TTRAJ += 1
                continue
            DATA = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
            DATA = DATA.readlines()
            TSTEPS = len(DATA)
            if TSTEPS >= TSCOL:
                POPH = np.zeros((TSCOL,NSTATES))
                POPC = np.zeros((TSCOL,NSTATES))
                INDEX = 0
                for LINE in DATA[0:TSCOL]:
                    VAL = LINE.split()
                    PES = np.int(VAL[0])
                    POPH[INDEX][PES-1] = 1.0
                    POPC[INDEX] = np.float_(VAL[2:2+NSTATES])
                    INDEX += 1
                FPOPH += POPH
                FPOPC += POPC
                print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                CTRAJ += 1
                if TSTEPS == TSMAX:
                    ETRAJ += 1
            else:
                print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                POPFLAG = 1
            TTRAJ += 1
    if CTRAJ == 0:
        print 'No trajectories were completed within %0*.2f.' % (len(str(TSMAX)),TCOLL)
    else:
        FPOPH = FPOPH/CTRAJ
        FPOPC = FPOPC/CTRAJ
        print 'Total Trajectories:', '%04d' % (TTRAJ)
        print 'Completed Trajectories:', '%04d' % (CTRAJ)
        print 'Excellent Trajectories:', '%04d' % (ETRAJ)
        print >> OUTPUT, 'Total Trajectories: ', '%04d' % (TTRAJ)
        print >> OUTPUT, 'Completed Trajectories: ', '%04d' % (CTRAJ)
        print >> OUTPUT, 'Excellent Trajectories: ', '%04d' % (ETRAJ)
        for TSTEP in np.arange(TSCOL):
            print >> OUTPUT, '%0*.2f' % (len(str((TSMAX))) + 2,DT*TSTEP), ' '.join(str('%.3f' % (x)) for x in FPOPH[TSTEP]), '%.3f' % (np.sum(FPOPH[TSTEP])), ' '.join(str('%.3f' % (x)) for x in FPOPC[TSTEP]), '%.3f' % (np.sum(FPOPC[TSTEP]))
    if POPFLAG == 1:
        print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pop.err.' % (len(str(TSMAX)),TCOLL)
    else:
        os.remove('%s/pop.err' % (CWD))

