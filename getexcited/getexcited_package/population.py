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

cwd = os.getcwd()

def population():

    print 'Collecting populations.'

    ## directory names ##
    NEXMDir = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print 'Path %s does not exist.' % (NEXMDir)
        sys.exit()

    ## User-defined length of analysis and initialize arrays ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print 'path %s/header does not exist.' % (NEXMDir)
        sys.exit()
    header = open('%s/header' % (NEXMDir),'r')
    for line in header:
        if 'n_exc_states_propagate' in line:
            nstates = np.int(line.split()[0][len('n_exc_states_propagate='):-1])
        if 'time_init' in line:
            tinith = np.float(line.split()[0][len('time_init='):-1])
        if 'time_step' in line:
            dt = np.float(line.split()[0][len('time_step='):-1])
        if 'n_class_steps' in line:
            tsmax = np.int(line.split()[0][len('n_class_steps='):-1]) + 1
        if 'out_data_steps' in line:
            odata = np.int(line.split()[0][len('out_data_steps='):-1])
    tcoll = input('Collect populations up to what time in femtoseconds: ')
    if isinstance(tcoll, int) == false and isinstance(tcoll, float) == false:
        print 'time must be integer or float.'
        sys.exit()
    if tcoll < 0:
        print 'Time must be integer or float greater than zero.'
        sys.exit()
    tcoll = np.float(tcoll)
    if tcoll > (tsmax - 1)*dt:
        tcoll = (tsmax - 1)*dt
    tscol = 0
    while tscol*dt*odata <= tcoll:
        tscol += 1
    times = np.around(np.linspace(tinith, tcoll, tscol), decimals = 3)
    fpoph = np.zeros((tscol,nstates))
    fpopc = np.zeros((tscol,nstates))

    ## Collect populations ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDir)
        sys.exit()
    output = open('%s/pop.out' % (cwd),'w')
    error = open('%s/pop.err' % (cwd),'w')
    ttraj = 0
    ctraj = 0
    etraj = 0
    errflag = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(dirlist1,int) == true:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                print >> error, '%s%04d/coeff-n.out' % (NEXMD,dir), 'does not exist'
                errflag = 1
                ttraj += 1
                continue
            data = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
            data = data.readlines()
            tsteps = len(data)
            tflag = 0
            if tsteps >= tscol:
                poph = np.zeros((tscol,nstates))
                popc = np.zeros((tscol,nstates))
                index = 0
                for line in data[0:tscol:1]:
                    val = line.split()
                    pes = np.int(val[0])
                    time = np.around(np.float(val[1]), decimals = 3)
                    if time != times[index]:
                        print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index])
                        tflag = 1
                        errflag = 1
                        break
                    poph[index][pes-1] = 1.0
                    popc[index] = np.float_(val[2:2+nstates])
                    index += 1
                if tflag == 0:
                    fpoph += poph
                    fpopc += popc
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ctraj += 1
                    if tsteps == tsmax:
                        etraj += 1
            else:
                print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                print >> error, '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                errflag = 1
            ttraj += 1
    if ctraj == 0:
        print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)),tcoll)
        os.remove('%s/pop.out' % (cwd))
    else:
        fpoph = fpoph/ctraj
        fpopc = fpopc/ctraj
        print 'Total trajectories:', '%04d' % (ttraj)
        print 'Completed trajectories:', '%04d' % (ctraj)
        print 'Excellent trajectories:', '%04d' % (etraj)
        print >> output, 'Total trajectories: ', '%04d' % (ttraj)
        print >> output, 'Completed trajectories: ', '%04d' % (ctraj)
        print >> output, 'Excellent trajectories: ', '%04d' % (etraj)
        for tstep in np.arange(tscol):
            print >> output, '%0*.2f' % (len(str((tsmax))) + 2,dt*tstep), ' '.join(str('%.3f' % (x)) for x in fpoph[tstep]), '%.3f' % (np.sum(fpoph[tstep])), ' '.join(str('%.3f' % (x)) for x in fpopc[tstep]), '%.3f' % (np.sum(fpopc[tstep]))
    if errflag == 1:
        print 'One or more trajectories have experienced an error, check pop.err.'
    else:
        os.remove('%s/pop.err' % (cwd))

