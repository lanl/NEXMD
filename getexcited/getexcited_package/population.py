#/usr/bin/python

'''

This function collects NEXMD populations from surface hopping
and populations from Schrodinger equation.

Two files are generated in the current working directory if
collecting NEXMD populations are requested, 'pop.out' and
'pop.err'.  In 'pop.out', first column is time in fs, followed by
the average population on each PES, followed by the sum of all
PES populations (should be = 1.0), followed by the average
populations from quantum coefficients, followed by the sum of
quantum populations (should be = 1.0).  In 'pop.err', first
column is directory of trajectory, third column is the time at
which the trajectory has ended in fs, or 'does not exist'.  The
'pop.err' file is not generated if all trajectories are complete.              

'Completed' trajectories are trajectories that have completed
within the time defined by user while executing this function.
'Excellent' trajectories are trajectories that have completed
within the time defined in 'header' located in the NEXMD
directory.                                                        

Output Files:
- pop_[type].out, where [type] = mean_ensemble

Error Files:
- pop_[type].err, where [type] = mean_ensemble

'''

import numpy as np
import os
import sys
import glob

cwd = os.getcwd()

def population(header):

    print('Collecting populations.')

    ## Directory names ##
    NEXMDir = input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print('Path %s does not exist.' % (NEXMDir))
        sys.exit()
    ## Check if NEXMD folders exist ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print('There are no NEXMD folders in %s.' % (NEXMDir))
        sys.exit()

    ## Information from header ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print('Path %s/header does not exist.' % (NEXMDir))
        sys.exit()
    header = header('%s/header' % (NEXMDir))
        
    ## Adding + 1 to include zeroth time-step ##
    header.n_class_steps = header.n_class_steps + 1

    ## Collection time ##
    tcoll = eval(input('Calculate populations up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: '))
    if isinstance(tcoll, int) == False and isinstance(tcoll, float) == False:
        print('Time must be integer or float.')
        sys.exit()
    if tcoll < 0:
        print('Time must be integer or float greater than zero.')
        sys.exit()
    tcoll = np.float(tcoll)
    if tcoll > (header.n_class_steps - 1)*header.time_step:
        tcoll = (header.n_class_steps -1)*header.time_step

    ## Number of classical time-steps ##
    tscol = 0
    while tscol*header.time_step*header.out_data_steps <= tcoll:
        tscol += 1
    
    ## Collection time array ##
    times = np.around(np.linspace(header.time_init, tcoll, tscol), decimals = 3)
    fpoph = np.zeros((tscol,header.n_exc_states_propagate))
    fpopc = np.zeros((tscol,header.n_exc_states_propagate))
    
    ## Collect populations ##
    output = open('%s/pop_mean_ensemble.out' % (cwd),'w')
    error = open('%s/pop_mean_ensemble.err' % (cwd),'w')
    ttraj = 0
    ctraj = 0
    etraj = 0
    errflag = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print('Path %sdirlist1 does not exist.' % (NEXMD))
            sys.exit()
        dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                errflag = 1
                ttraj += 1
                continue
            data = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
            data = data.readlines()
            tsteps = len(data)
            tflag = 0
            if tsteps >= tscol:
                poph = np.zeros((tscol,header.n_exc_states_propagate))
                popc = np.zeros((tscol,header.n_exc_states_propagate))
                index = 0
                for line in data[0:tscol:1]:
                    val = line.split()
                    pes = np.int(val[0])
                    time = np.around(np.float(val[1]), decimals = 3)
                    if time != times[index]:
                        print('There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index]), file=error)
                        tflag = 1
                        errflag = 1
                        break
                    poph[index][pes-1] = 1.0
                    popc[index] = np.float_(val[2:2+header.n_exc_states_propagate])
                    index += 1
                if tflag == 0:
                    fpoph += poph
                    fpopc += popc
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ctraj += 1
                    if tsteps == header.n_class_steps:
                        etraj += 1
            else:
                print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step), file=error)
                errflag = 1
            ttraj += 1
    if ctraj == 0:
        print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)),tcoll))
        os.remove('%s/pop_mean_ensemble.out' % (cwd))
    else:
        fpoph = fpoph/ctraj
        fpopc = fpopc/ctraj
        print('Total trajectories:', '%04d' % (ttraj))
        print('Completed trajectories:', '%04d' % (ctraj))
        print('Excellent trajectories:', '%04d' % (etraj))
        print('Total trajectories: ', '%04d' % (ttraj), file=output)
        print('Completed trajectories: ', '%04d' % (ctraj), file=output)
        print('Excellent trajectories: ', '%04d' % (etraj), file=output)
        for tstep in np.arange(tscol):
            print('%0*.2f' % (len(str((header.n_class_steps))) + 2,header.time_step*tstep), ' '.join(str('%.3f' % (x)) for x in fpoph[tstep]), '%.3f' % (np.sum(fpoph[tstep])), ' '.join(str('%.3f' % (x)) for x in fpopc[tstep]), '%.3f' % (np.sum(fpopc[tstep])), file=output)
    if errflag == 1:
        print('One or more trajectories have experienced an error, check pop_mean_ensemble.err.')
    else:
        os.remove('%s/pop_mean_ensemble.err' % (cwd))
