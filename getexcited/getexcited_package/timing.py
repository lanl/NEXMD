#/usr/bin/python

'''

This function collects timings from all trajectories.

If this function is requested, the timings located at the end of the
standard output files (i.e. md.out) are outputted to a file called
'timing.out'.  The first column is directory of the trajectory,
followed by its total CPU time, and timings for the ground state,
excited states, adiabatic forces, and non-adiabatic derivative coupling,
respectively.  These timings, averaged over all trajectories, are also
printed to screen when this function is executed.  An error file called
'timing.err' will be generated if any problems occur such as
non-existent or incomplete files.

'''

import numpy as np
import os
import sys
import subprocess
import shlex
import glob

cwd = os.getcwd()

def timing(pathtotime):

    print('Collecting timings from trajectories.')

    ## Directory names ##
    NEXMDir = input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print('Path %s does not exist.' % (NEXMDir))
        sys.exit()
    
    ## Collect and check timings ##
    print('Collecting timings. please wait ...')
    if not os.path.exists('%s/getexcited_package/collectime.sh' % (pathtotime)):
        print('The script, collectime.sh, must be in the getexcited_package.')
        sys.exit()
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print('There are no NEXMD folders in %s.' % (NEXMDir))
        sys.exit()
    error = open('%s/timing.err' % (cwd),'w')
    errflag = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/%s/dirlist1' % (cwd,NEXMD)):
            print('Path %s/%sdirlist1 does not exist.' % (cwd,NEXMD))
            sys.exit()
        dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            if not os.path.exists('%s/%s/%04d' % (cwd,NEXMD,dir)):
                print('%s%04d' % (NEXMD,dir), 'does not exist', file=error)
                errflag = 1
                continue
            os.chdir('%s/%s/%04d' % (cwd,NEXMD,dir))
            if not os.path.exists('%s/%s/%04d/md.out' % (cwd,NEXMD,dir)):
                print('%s%04d/md.out' % (NEXMD,dir), 'does not exist', file=error)
                errflag = 1
                continue
            subprocess.call(shlex.split('sh %s/getexcited_package/collectime.sh' % (pathtotime)))
            if not os.path.exists('%s/%s/%04d/timing.out' % (cwd,NEXMD,dir)):
                print('%s/%04d/timing.out' % (NEXMD,dir), 'does not exist', file=error)
                errflag = 1
                continue
            with open('%s/%s/%04d/timing.out' % (cwd,NEXMD,dir),'r') as data:
                data = data.readlines()
                if len(data) != 6 or 'MD total CPU time' not in data[0]:
                    print('%s%04d/timing.out' % (NEXMD,dir), 'is incomplete', file=error)
                    errflag = 1
            print('%s%04d' % (NEXMD,dir))
    if errflag == 1:
        print('One or more trajectories did not finish, check timing.err.')
        contq = eval(input('Continue? Answer yes [1] or no [0]: '))
        if contq not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if contq == 0:
            sys.exit()
    else:
        os.remove('%s/timing.err' % (cwd))

    ## Extract and combine timings ##

    print(' Collecting time....')
    timing = open('%s/timing.out' % (cwd),'w')
    times = np.zeros(5) ## change 5 to 1 for old code (NAESMD)
    traj = 0
    for NEXMD in NEXMDs:
        dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            try: 
                data = open('%s/%s/%04d/timing.out' % (cwd,NEXMD,dir),'r')
            except:
                print(('%s/%s/%04d/timing.out does not exist' % (cwd,NEXMD,dir)))
                continue
            data = data.readlines()
            data = np.delete(data, (1), axis = 0) ## comment out for old code (NAESMD)
            tarray = np.array([])
            index = 0
            fail = 0
            for line in data:
                val = line.split()
                if index == 0:
                    try: 
                        tarray = np.append(tarray, np.float(val[5]))
                    except:
                        print((' %s/%s/%04d did not finish' %(cwd,NEXMD,dir)))
                        fail = 1
                        traj -= 1
                        break
                else: ## comment out for old code (NAESMD)
                    tarray = np.append(tarray, np.float(val[0]))
                index += 1
            if fail == 0:
                times += tarray
                print('%s%04d' % (NEXMD,dir), ' '.join(str('%06d' % (x)) for x in tarray), file=timing)
            os.remove('%s/%s/%04d/timing.out' % (cwd,NEXMD,dir))
            print('%s%04d' % (NEXMD,dir))
            traj += 1
    times = times/traj
    print('Total number of trajectories: %d' %(traj))
    print('Mean total cpu [s]:', '%06d' % (times[0]))
    ## comment all below for old code (NAESMD) ##
    print('Mean ground state [s]:', '%06d' % (times[1]))
    print('Mean excited states [s]:', '%06d' % (times[2]))
    print('Mean adiabatic forces [s]:', '%06d' % (times[3]))
    print('Mean non-adiabatic derivatives [s]:', '%06d' % (times[4]))
