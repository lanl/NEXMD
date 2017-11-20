#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares restart input files for NEXMD.             |
|                                                                   |
| If restart input files are requested, the function searches for   |
| trajectories that did not complete up to the user-defined         |
| number of classical steps and prepares input files with the name  |
| 'input.ceon' in their respective directories.  The input file is  |
| prepared with coordinates, velocities, quantum coefficients, and  |
| last-residing surface taken from the 'restart.out' file.  These   |
| are data from the last-generated time-step.  Within each 'NEXMD#' |
| folder, with directories that contain restart input files, a new  |
| 'dirlist' will be generated with a list of restart directories.   |
| An error file, 'restart.err', lists the directories that either   |
| (1) do not contain 'energy-ev.out' files, from which the last     |
| time-step is determined, or (2) have incomplete 'restart.out'     |
| files.  The former generally means the trajectory did not start   |
| when NEXMD was first attempted.  The 'restart.err' file is not    |
| generated if there are no such trajectories.  During every        |
| iteration of requesting restart input files, a file containing    |
| the random seeds is generated with the name 'rseedslist#'.  Part  |
| of this function deletes the data at all time-steps after the     |
| time-step of the 'restart.out' file.  The purpose of this is to   |
| have continuous set of data along the trajectory with no repeated |
| information.  An error file called 'delextra.out' is generated if |
| one or more output files do not exist.                            |
|___________________________________________________________________|

'''


import numpy as np
import os
import sys
import re
import glob
import filecmp
import subprocess
import shlex

cwd = os.getcwd()

def extract(file):
    num = re.findall('\d+$', file)
    return (np.int(num[0]) if num else 0, file)

def restart(pathtodel):

    print 'Preparing restart input files for NEXMD.'

    ## Directory names ##
    NEXMDir = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print 'Path %s does not exist.' % (NEXMDir)
        sys.exit()

    ## Choose classical time-steps, get number of quantum steps and verbosity ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print 'Path %s/header does not exist.' % (NEXMDir)
        sys.exit()
    header = open('%s/header' % (NEXMDir),'r')
    header = header.readlines()
    num = 0
    tline = len(header)
    verb = none
    for line in header:
        if 'bo_dynamics_flag' in line:
            boflag = np.int(line.split()[0][len('bo_dynamics_flag='):-1])
        if 'time_step' in line:
            dt = np.float(line.split()[0][len('time_step='):-1])
        if 'n_class_steps' in line:
            tsmax = np.int(line.split()[0][len('n_class_steps='):-1])
        if 'n_quant_steps' in line:
            nqstep = np.int(line.split()[0][len('n_quant_steps='):-1])
            if nqstep == 0:
                nqstep = 1
        if '&moldyn' in line:
            tline = num
        if 'verbosity' in line and num > tline and verb is none:
            verb = np.int(line.split()[0][len('verbosity='):-1])
        if 'out_data_steps' in line:
            odata = np.int(line.split()[0][len('out_data_steps='):-1])
            if odata == 0:
                print 'No data has been printed to files because out_data_steps = 0 in header.'
                sys.exit()
        num += 1
    print 'Currently, trajectories are set to run for %d classical steps with a time-step of %.2f fs.\nthis is a total of %.2f fs.' % (tsmax,dt,tsmax*dt)
    tsmaxq = input('Keep this trajectory length? Answer yes [1] or no [0]: ')
    if tsmaxq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if tsmaxq == 0:
        ntsmax = input('Enter new number of classical time-steps: ')
        if isinstance(ntsmax, int) == false:
            print 'Answer must be integer.'
            sys.exit()
        if ntsmax <= tsmax:
            print 'Answer must be greater than or equal to the previous number of classical steps used, which was %d.\nTo reduce number of classical steps past %d, simply change n_class_steps in header.' % (tsmax,tsmax)
            sys.exit()
        nheader = open('%s/nheader' % (NEXMDir),'w')
        for line in header:
            if 'n_class_steps' in line:
                nheader.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (ntsmax))
            else:
                nheader.write(line)
        nheader.close()
        tsmax = ntsmax
        os.rename('%s/nheader' % (NEXMDir), '%s/header' % (NEXMDir))
    
    ## Choose random seeds ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDir)
        sys.exit()
    ntraj = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'Path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        data = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(data,int) == true:
            data = np.array([data])
        ntraj += len(data)
    randq = input('New random seeds? Answer yes [1] or no [0]: ')
    if randq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    data = glob.glob('%s/rseedslist*' % (NEXMDir))
    if len(data) == 0:
        print 'There are no rseedslists in %s.' % (NEXMDir)
        sys.exit()
    maxdir = np.int(extract(max(data,key = extract))[0])
    if randq == 0:
        rseeds = raw_input('Path to random-seeds list (** must be different from past random seeds **): ')
        if not os.path.exists(rseeds):
            print 'Path %s does not exist.' % (rseeds)
            sys.exit()
        rseeds = np.int_(np.genfromtxt('%s' % (rseeds)))
        if isinstance(rseeds,int) == true:
            rseeds = np.array([rseeds])
        len = len(rseeds)
        if len < ntraj:
            print 'Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories is %d.' % (len,ntraj)
            sys.exit()
        for rseed in rseeds:
            if rseed < 0:
                print 'A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (rseed)
                sys.exit()

    ## Prepare NEXMD restart input files ##
    print 'Searching for incomplete trajectories and preparing input files.'
    header = open('%s/header' % (NEXMDir),'r')
    header = header.readlines()
    rseedslist = open('%s/rseedslist%d' % (NEXMDir,maxdir + 1),'w')
    error = open('%s/restart.err' % (cwd),'w')
    rtimes = np.array([])
    rstflag = 0
    traj = 0
    for NEXMD in NEXMDs:
        dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(dirlist1,int) == true:
            dirlist1 = np.array([dirlist1])
        dirlist = open('%s/dirlist' % (NEXMD),'w')
        for dir in dirlist1:
            if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                print >> error, '%s%04d/energy-ev.out' % (NEXMD,dir), 'does not exist'
                print >> rseedslist, '%d' % (-123456789)
                rstflag = 1
                traj += 1
                continue
            data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
            data = data.readlines()
            tsteps = len(data) - 2
            if tsteps != tsmax:
                if not os.path.exists('%s/%04d/restart.out' % (NEXMD,dir)):
                    print >> error, '%s/%04d/restart.out' % (NEXMD,dir), 'does not exist'
                    print >> rseedslist, '%d' % (-123456789)
                    rstflag = 1
                    traj += 1
                    continue
                data = open('%s/%04d/restart.out' % (NEXMD,dir),'r')
                data = data.readlines()
                index = 0
                array = np.array([])
                for line in data:
                    if 'time' in line:
                        time = np.around(np.float(line.split()[-1]), decimals = 3)
                    if 'state' in line:
                        state = np.int(line.split()[-1])
                    if 'seed' in line:
                        rseed = np.int(line.split()[-1])
                    if '$coord' in line:
                        array = np.append(array,index)
                    if '$endcoord' in line:
                        array = np.append(array,index)
                    if '$veloc' in line:
                        array = np.append(array,index)
                    if '$endveloc' in line:
                        array = np.append(array,index)
                    if '$coeff' in line:
                        array = np.append(array,index)
                    if '$endcoeff' in line:
                        array = np.append(array,index)
                    index += 1
                array = np.int_(array)
                if len(array) != 6:
                    print >> error, '%s%04d/restart.out' % (NEXMD,dir), 'is incomplete'
                    print >> rseedslist, '%d' % (-123456789)
                    rstflag = 1
                    traj += 1
                    continue
                coords = data[array[0]:array[1]+1:1]
                velocs = data[array[2]:array[3]+1:1]
                coeffs = data[array[4]:array[5]+1:1]
                ## start renormalize coefficients ##
                ncoeffs = np.zeros(((len(coeffs) - 2), 2))
                index = 0
                for line in coeffs[1:len(coeffs) - 1:1]:
                    ncoeffs[index] = np.float_(line.split())
                    index += 1
                ncoeffs[:,0] = ncoeffs[:,0]/np.sum(ncoeffs[:,0])
                ## end renormalize coefficients ##
                data = glob.glob('%s/view*' % (NEXMDir))
                data = [ x[:-10] for x in data ]
                if len(data) != 0:
                    max = np.int(extract(max(data,key = extract))[0])
                else:
                    max = 0
                input = open('%s/%04d/input.ceon' % (NEXMD,dir),'w')
                for line in header:
                    if 'rnd_seed' in line:
                        input.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseed if randq == 1 else rseeds[traj]))
                    else:
                        if 'exc_state_init_flag' in line:
                            input.write('   exc_state_init=%d, ! initial excited state (0 - ground state) [0]\n' % (state))
                        else:
                            if 'time_init' in line:
                                input.write('   time_init=%.1f, ! initial time, fs [0.00]\n' % (time))
                            else:
                                if 'n_class_steps' in line:
                                    input.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (np.int(tsmax-time/dt)))
                                else:
                                    if 'out_count_init' in line:
                                        input.write('   out_count_init=%d, ! initial count for output files [0]\n' % (max))
                                    else:
                                        if 'nucl_coord_veloc' in line:
                                            for line in coords:
                                                input.write(line)
                                            input.write('\n')
                                            for line in velocs:
                                                input.write(line)
                                        else:
                                            if 'quant_amp_phase' in line:
                                                input.write('$coeff\n')
                                                for line in ncoeffs:
                                                    input.write('  %.10f  %.10f\n' % (line[0],line[1]))
                                                input.write('$endcoeff\n')
                                            else:
                                                input.write(line)
                rtimes = np.append(rtimes,time)
                print >> dirlist, '%04d' % (dir)
                print '%s%04d' % (NEXMD,dir)
                print >> rseedslist, '%d' % (rseed if randq == 1 else rseeds[traj])
            else:
                print >> rseedslist, '%d' % (-123456789)
            traj += 1
        dirlist.close()
    if filecmp.cmp('%s/rseedslist%d' % (NEXMDir,maxdir), '%s/rseedslist%d' % (NEXMDir,maxdir + 1)):
        os.remove('%s/rseedslist%d' % (NEXMDir,maxdir + 1))
    if rstflag == 1:
        print 'One or more trajectories cannot be restarted, check restart.err.'
    else:
        os.remove('%s/restart.err' % (cwd))

    ## Delete extraneous data in output files ##
    if rstflag == 1:
        contq = input('Continue to delete extraneous data in output files? Answer yes [1] or no [0]: ')
        if contq not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if contq == 0:
            sys.exit()
    print 'Deleting extraneous data in output files. please wait ...'
    if boflag == 0:
        files = np.array(['coeff-n.out', 'energy-ev.out', 'nact.out', 'order.out', 'pes.out', 'temperature.out', 'transition-densities.out'])
        ofiles = np.array(['hops.out', 'hops-trial.out'])
    else:
        files = np.array(['energy-ev.out', 'temperature.out', 'pes.out','transition-densities.out'])
        ofiles = np.array([])
    error = open('%s/delextra.err' % (cwd),'w')
    rstflag = 0
    traj = 0
    for NEXMD in NEXMDs:
        if os.stat('%s/%s/dirlist' % (cwd,NEXMD)).st_size == 0:
            continue
        else:
            dirlist = np.int_(np.genfromtxt('%s/%s/dirlist' % (cwd,NEXMD)))
        if isinstance(dirlist,int) == true:
            dirlist = np.array([dirlist])
        for dir in dirlist:
            if not os.path.exists('%s/%s/%04d' % (cwd,NEXMD,dir)):
                print >> error, '%s/%s%04d' % (cwd,NEXMD,dir), 'does not exist'
                rstflag = 1
                traj += 1
                continue
            os.chdir('%s/%s/%04d' % (cwd,NEXMD,dir))
            for index in np.arange(len(files)):
                if not os.path.exists('%s/%s/%04d/%s' % (cwd,NEXMD,dir,files[index])):
                    print >> error, '%s/%s%04d/%s' % (cwd,NEXMD,dir,files[index]), 'does not exist'
                    rstflag = 1
                    continue
                ## Derivation of the following algorithm is provided at the end of this script ##
                data = subprocess.check_output(['tail','-1','%s' % (files[index])])
                ltime = np.around(np.float(np.fromstring(data,dtype=float,sep=' ')[1 if index == 0 and boflag == 0 else 0]), decimals = 3)
                if rtimes[traj] > ltime + odata*dt:
                    print >> error, 'last time-step in', '%s/%s%04d/%s' % (cwd,NEXMD,dir,'restart.out'), 'exceeds last time-step in', '%s/%s%04d/%s' % (cwd,NEXMD,dir,files[index])
                    rstflag = 1
                    continue
                ncsteps = 0
                while rtimes[traj] + ncsteps*(odata*dt) <= ltime:
                    ncsteps += 1
                if verb == 3 and index in [2,4]:
                    lctime = rtimes[traj] + (ncsteps - 1)*(odata*dt)
                    nqsteps = 0
                    while lctime*nqstep + nqsteps*dt <= ltime*nqstep:
                        nqsteps += 1
                    if ltime == lctime:
                        nlines = (ncsteps - 1)*(odata*(nqstep - 1) + 1) + (nqsteps - 1) + 1
                    else:
                        ncsteps = 0
                        while rtimes[traj] + ncsteps*dt <= ltime:
                            ncsteps += 1
                        nlines = (ncsteps - 1)*(odata*(nqstep - 1) + 1) + (nqsteps - 1) - (ncsteps - 1) + 1
                else:
                    nlines = (ncsteps - 1) + 1
                if boflag == 0:
                    subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtodel,(nlines - 1 if index in [2,3,5] else nlines),files[index])))
                else:
                    subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtodel,(nlines - 1 if index == 1 else nlines),files[index])))
                os.rename('%s.restart' % (files[index]), '%s' % (files[index]))
            for index in np.arange(len(ofiles)):
                if not os.path.exists('%s/%s/%04d/%s' % (cwd,NEXMD,dir,ofiles[index])):
                    print >> error, '%s/%s%04d/%s' % (cwd,NEXMD,dir,ofiles[index]), 'does not exist'
                    rstflag = 1
                    continue
                ltime = 1000000
                nlines = 0
                while ltime >= rtimes[traj]:
                    data = subprocess.check_output(['tail','%d' % (-(nlines + 1)),'%s' % (ofiles[index])])
                    ltime = np.float(np.fromstring(data,dtype=float,sep=' ')[0])
                    nlines += 1
                subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtodel, nlines - 1, ofiles[index])))
                os.rename('%s.restart' % (ofiles[index]), '%s' % (ofiles[index]))
            print '%s%04d' % (NEXMD,dir)
            traj += 1
    if rstflag == 1:
        print 'One or more trajectories cannot be restarted properly, check delextra.err.'
    else:
        os.remove('%s/delextra.err' % (cwd))
    print 'A total of %d trajectories have been prepared for restart.' % (traj)

'''
    
    dt     = time-step
    verb   = verbosity of molecular dynamics data
    odata  = data are printed every odata classical time-steps
    ltime  = last printed time-step (may be classical or quantum)
    lctime = last printed classical time-step
    rtime  = restart time-step (always classical)
    nqstep = number of quantum steps per classical step
    
    a = (lctime - rtime)/dt                                      = total # of classical steps computed
    b = (lctime - rtime)/(odata x dt)                            = total # of classical steps printed
    c = (lctime - rtime)/dt x ((odata - 1)/odata)                = total # of classical steps not printed
    d = (nqstep - 1) x (lctime - rtime)/dt x ((odata - 1)/odata) = total # of quantum steps from unprinted classical steps
    e = nqstep x (lctime - rtime)/(odata x dt)                   = total # of classical and quantum steps from printed classical steps
    
    if verb != 3:
        total # of lines to delete = b
    if verb == 3 and file not in [nact.out, pes.out]:
        total # of lines to delete = b
    if verb == 3 and file in [nact.out, pes.out]:
        total # of lines to delete = d + e + resdiual quantum steps after lctime
        if ltime == lctime:
            residual quantum steps after lctime = nqstep x (ltime - lctime)/dt
            total # of lines to delete = d + e + nqstep x (ltime - lctime)/dt + 1 ( + 1 to delete the restart time-step )
            total # of lines to delete = (lctime - rtime)*(odata*(nqstep - 1) + 1)/(odata*dt) + nqstep*(ltime - lctime)/dt + 1
            total # of lines to delete = (ncsteps - 1)*(odata*(nqstep - 1) + 1) + (nqsteps - 1) + 1
            *** see in code how to calculate ncsteps and nqsteps ***
        if ltime != lctime:
            residual quantum steps after lctime = nqstep x (ltime - lctime)/dt
            number of unprinted classical steps after lctime = ncsteps - 1
            total # of lines to delete = d + e + nqstep x (ltime - lctime)/dt - (ncsteps - 1) + 1 ( + 1 to delete the restart time-step)
            total # of lines to delete = (lctime - rtime)*(odata*(nqstep - 1) + 1)/(odata*dt) + nqstep*(ltime - lctime)/dt - (ncsteps - 1) + 1
            total # of lines to delete = (ncsteps - 1)*(odata*(nqstep - 1) + 1) + (nqsteps - 1) - (ncsteps - 1) + 1
            *** see in code how to calculate ncsteps and nqsteps ***
            
'''
