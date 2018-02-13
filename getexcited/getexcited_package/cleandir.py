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

def cleandir(header):
    
    print 'Cleaning directories of unfinished trajectories.'
    
    ## Check to delete question ##
    checkq = input('Are you sure you want to delete all unfinished trajectories?\nAnswer yes [1] or no [0]: ')
    if checkq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if checkq == 0:
        sys.exit()

    ## Directory names ##
    NEXMDir = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print 'Path %s does not exist.' % (NEXMDir)
        sys.exit()
    ## Check if NEXMD folders exist ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDir)
        sys.exit()

    ## Information from header ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print 'Path %s/header does not exist.' % (NEXMDir)
        sys.exit()
    header = header('%s/header' % (NEXMDir))

    ## Adding + 1 to include zeroth time-step ##
    header.n_class_steps = header.n_class_steps + 1

    ## Ask user to delete unfinished trajectories up to user-defined number of classical time-steps ##
    contq = input('Trajectories less than %d classical time-steps will be deleted.\nContinue? Answer yes [1] or no [0]: ' % (header.n_class_steps - 1))
    if contq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if contq == 0:
        sys.exit()

    ## Clean directories ##
    clnflag = 0
    traj = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'Path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        dirlist = open('%s/dirlist' % (NEXMD),'w')
        for dir in dirlist1:
            if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                clnflag = 1
            else:
                data = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                data = data.readlines()
                tsteps = len(data)
                if tsteps != header.n_class_steps:
                    clnflag = 1
            if clnflag == 1:
                files = next(os.walk('%s/%04d' % (NEXMD,dir)))[1]
                for file in files:
                    shutil.rmtree('%s/%04d/%s' % (NEXMD,dir,file))
                files = next(os.walk('%s/%04d' % (NEXMD,dir)))[2]
                for file in files:
                    if file != 'input.ceon':
                        os.remove('%s/%04d/%s' % (NEXMD,dir,file))
                if not os.path.exists('%s/%04d/input.ceon' % (NEXMD,dir)):
                    print 'path %s%04d/input.ceon does not exist.' % (NEXMD,dir)
                    sys.exit()
                old_inputfile = open('%s/%04d/input.ceon' % (NEXMD,dir),'r')
                new_inputfile = open('%s/%04d/ninput.ceon' % (NEXMD,dir),'w')
                for line in old_inputfile:
                    if 'n_class_steps' in line:
                        new_inputfile.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (header.n_class_steps - 1))
                    else:
                        new_inputfile.write(line)
                os.rename('%s/%04d/ninput.ceon' % (NEXMD,dir), '%s/%04d/input.ceon' % (NEXMD,dir))
                clnflag = 0
                traj += 1
                print >> dirlist, '%04d' % (dir)
                print '%s%04d' % (NEXMD,dir)
    print 'The contents of %d trajectories have been deleted.' % (traj)
