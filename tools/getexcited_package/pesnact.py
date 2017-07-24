#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function collects PESs and NACTs from NEXMD.                 |
|                                                                   |
| A total of four output files will be generated in the current     |
| working directory if collecting PESs and NACTs is requested.  In  |
| 'pesall.out' and 'nactall.out', the PESs and NACTs at all         |
| time-steps are shown, respectively.  In 'pesall.out', columns     |
| from left to right are: directory of trajectory, current state,   |
| new state, followed by all PESs.  Likewise, 'nactall.out' is in   |
| similar format.  The columns showing NACTs are consecutive rows   |
| of the NACT matrix, same as that shown in 'nact.out', located in  |
| the directory of each trajectory.  In 'peshop.out' and            |
| 'nacthop.out' are same data as '...all.out', but only at          |
| time-steps where hops occur.  An error file will be generated if  |
| certain files do not exist or if trajectories did not finish up   |
| to the user-defined length of analysis.                           |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import glob

CWD = os.getcwd()

def PESNACT():

    print 'Collecting PESs and NACTs.'

## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()

## USER-DEFINED LENGTH AND TYPE OF ANALYSIS ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header'% (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    NUM = 0
    TLINE = len(HEADER)
    VERB = None
    for LINE in HEADER:
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1]) + 1
        if 'n_quant_steps' in LINE:
            NQSTEP = np.int(LINE.split()[0][len('n_quant_steps='):-1])
            if NQSTEP == 0:
                NQSTEP = 1
        if '&moldyn' in LINE:
            TLINE = NUM
        if 'verbosity' in LINE and NUM > TLINE and VERB is None:
            VERB = np.int(LINE.split()[0][len('verbosity='):-1])
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
            if ODATA == 0:
                print 'No data has been printed to files because out_data_steps = 0 in header.'
                sys.exit()
        NUM += 1
    TCOLL = input('Collect PESs and NACTs up to what time in femtoseconds: ')
    if isinstance(TCOLL, int) == False and isinstance(TCOLL, float) == False:
        print 'Time must be integer or float.'
        sys.exit()
    if TCOLL < 0:
        print 'Time must be integer or float greater than zero.'
        sys.exit()
    TCOLL = np.float(TCOLL)
    NSTEPS = 0
    while NSTEPS*DT < TCOLL:
        NSTEPS += 1
    TCOLL = (NSTEPS - 1)*DT
    if TCOLL > (TSMAX - 1)*DT:
        TCOLL = (TSMAX - 1)*DT
    TSCOL = 0
    while TSCOL*DT*ODATA <= TCOLL:
        TSCOL += 1
    NACTQ = input('Collect PESs and NACTs from ALL CLASSICAL TIME-STEPS [1] or AT HOPS ONLY [2]: ')
    if NACTQ not in [1,2]:
        print 'Answer must be 1 or 2.'
        sys.exit()
    if VERB == 3:
        if NACTQ == 1 and NQSTEP != 1:
            LINENUMS = np.arange(NQSTEP, TSCOL*NQSTEP - (NQSTEP - 1), NQSTEP)
        else:
            LINENUMS = np.arange(1, TSCOL)
    else:
        LINENUMS = np.arange(1, TSCOL)

## INDICES TO CUT EXTRANEOUS NACT DATA ##
    INDICES = np.array([])
    INDEX = NSTATES
    for TERM in np.split(np.arange(NSTATES*NSTATES),NSTATES):
        INDICES = np.append(INDICES,TERM[-INDEX::])
        INDEX -= 1
    INDICES = np.int_(np.insert(INDICES + 1, 0, 0, 0))

## COLLECT PESs AND NACTs ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()
    if NACTQ == 1:
        PESALL = open('%s/pesall.out' % (CWD),'w')
        NACTALL = open('%s/nactall.out' % (CWD),'w')
    PESHOP = open('%s/peshop.out' % (CWD),'w')
    NACTHOP = open('%s/nacthop.out' % (CWD),'w')
    ERROR = open('%s/pesnact.err' % (CWD),'w')
    TTRAJ = 0
    CTRAJ = 0
    ETRAJ = 0
    ERRFLAG = 0
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
                ERRFLAG = 1
                TTRAJ += 1
                continue
            HOPS = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
            HOPS = HOPS.readlines()
            TSTEPS = len(HOPS)
            if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                ERRFLAG = 1
                TTRAJ += 1
                continue
            PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
            PES = PES.readlines()
            if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                ERRFLAG = 1
                TTRAJ += 1
                continue
            NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
            NACT = NACT.readlines()
            if TSTEPS >= TSCOL:
                CSTATE = np.int(HOPS[0].split()[0])
                for LINE in LINENUMS:
                    NSTATE = np.int(HOPS[LINE].split()[0])
                    PESS = np.float_(PES[LINE].split())
                    NACTS = np.float_(NACT[LINE-1].split())[INDICES]
                    if NACTQ == 1:
                        print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                        print >> NACTALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                    if NSTATE != CSTATE:
                        print >> PESHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                        print >> NACTHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                    CSTATE = NSTATE
                print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                CTRAJ += 1
                if TSTEPS == TSMAX:
                    ETRAJ += 1
            else:
                print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                ERRFLAG = 1
            TTRAJ += 1

## PRINT SUMMARY ##
    print 'Total Trajectories:', '%04d' % (TTRAJ)
    print 'Completed Trajectories:', '%04d' % (CTRAJ)
    print 'Excellent Trajectories:', '%04d' % (ETRAJ)
    if ERRFLAG == 1:
        print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(TSMAX)),TCOLL)
    else:
        os.remove('%s/pesnact.err' % (CWD))
