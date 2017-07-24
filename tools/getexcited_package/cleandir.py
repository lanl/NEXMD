#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function cleans out the directories of unfinished            |
| trajectories.  At times, you may experience strange behavior on   |
| you HPC cluster.  A complete restart of unfinished trajectories   |
| may be desired when this occurs.                                  |
|                                                                   |
| If cleaning out the directories of unfinished trajectories is     |
| requested, this function determines the number of classical steps |
| from the 'header' in the NEXMD directory and searches all         |
| trajectories.  Completed trajectories are not affected.  In the   |
| directories of unfinished trajectories, all files are deleted     |
| except for the original 'input.ceon' file.  This option should    |
| only be used if the restart option is not more suitable!          |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import shutil
import glob

def CLEANDIR():
    
    print 'Cleaning directories of unfinished trajectories.'
    
## CHECK QUESTION ##
    CHECKQ = input('Are you sure you want to delete all unfinished trajectories?\nAnswer YES [1] or NO [0]: ')
    if CHECKQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if CHECKQ == 0:
        sys.exit()

## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()

## DETERMINE CLASSICAL TIME-STEPS ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header' % (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1]) + 1
            break
    CONTQ = input('Trajectories less than %d classical time-steps will be deleted.\nContinue? Answer YES [1] or NO [0]: ' % (TSMAX - 1))
    if CONTQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if CONTQ == 0:
        sys.exit()

## CLEAN DIRECTORIES ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()
    CLNFLAG = 0
    TRAJ = 0
    for NEXMD in NEXMDS:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'Path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        DIRLIST = open('%s/dirlist' % (NEXMD),'w')
        for DIR in DIRLIST1:
            if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,DIR)):
                CLNFLAG = 1
            else:
                DATA = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                TSTEPS = len(DATA)
                if TSTEPS != TSMAX:
                    CLNFLAG = 1
            if CLNFLAG == 1:
                FILES = next(os.walk('%s/%04d' % (NEXMD,DIR)))[1]
                for FILE in FILES:
                    shutil.rmtree('%s/%04d/%s' % (NEXMD,DIR,FILE))
                FILES = next(os.walk('%s/%04d' % (NEXMD,DIR)))[2]
                for FILE in FILES:
                    if FILE != 'input.ceon':
                        os.remove('%s/%04d/%s' % (NEXMD,DIR,FILE))
                if not os.path.exists('%s/%04d/input.ceon' % (NEXMD,DIR)):
                    print 'Path %s%04d/input.ceon does not exist.' % (NEXMD,DIR)
                    sys.exit()
                INPUT = open('%s/%04d/input.ceon' % (NEXMD,DIR),'r')
                NINPUT = open('%s/%04d/ninput.ceon' % (NEXMD,DIR),'w')
                for LINE in INPUT:
                    if 'n_class_steps' in LINE:
                        NINPUT.write('   n_class_steps=%d, ! Number of classical steps [1]\n' % (TSMAX - 1))
                    else:
                        NINPUT.write(LINE)
                os.rename('%s/%04d/ninput.ceon' % (NEXMD,DIR), '%s/%04d/input.ceon' % (NEXMD,DIR))
                CLNFLAG = 0
                TRAJ += 1
                print >> DIRLIST, '%04d' % (DIR)
                print '%s%04d' % (NEXMD,DIR)
    print 'The contents of %d trajectories have been deleted.' % (TRAJ)
