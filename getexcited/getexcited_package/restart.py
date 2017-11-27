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

def restart(pathtopack,header):

    print 'Preparing restart input files for NEXMD.'

    ## Type of calculation and directory check ##
    dynq = input('Restart a single trajectory or an ensemble of trajectories?\nAnswer one [1] or ensemble [0]: ')
    if dynq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if dynq == 0: ## ensemble
        NEXMDir = raw_input('Ensemble directory [e.g. NEXMD]: ')
        if not os.path.exists(NEXMDir):
            print 'Path %s does not exist.' % (NEXMDir)
            sys.exit()
        ## Check if NEXMD folders exist ##
        NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
        NEXMDs.sort()
        if len(NEXMDs) == 0:
            print 'There are no NEXMD folders in %s.' % (NEXMDir)
            sys.exit()
    if dynq == 1: ## single trajectory
        NEXMDir = raw_input('Single trajectory directory: ')
        if not os.path.exists(NEXMDir):
            print 'Path %s does not exist.' % (NEXMDir)
            sys.exit()

    ## Information from header ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print 'Path %s/header does not exist.' % (NEXMDir)
        sys.exit()
    header = header('%s/header' % (NEXMDir))

    ## Check output data ##
    if header.out_data_steps == 0:
        print 'No data have been printed to files because out_data_steps = 0 in %s/header.' % (NEXMDir)
        sys.exit()

    ## Redefine number of quantum steps per classical step if set to 0 in header ##
    if header.n_quant_steps == 0:
        header.n_quant_steps = 1

## Figure out what this part is about
    print 'Currently, trajectories are set to run for %d classical steps with a time-step of %.2f fs.\nThis is a total of %.2f fs.' % (header.n_class_steps,header.time_step,header.n_class_steps*header.time_step)
    tsmaxq = input('Keep this trajectory length? Answer yes [1] or no [0]: ')
    if tsmaxq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if tsmaxq == 0:
        ntsmax = input('Enter new number of classical time-steps: ')
        if isinstance(ntsmax, int) == False:
            print 'Answer must be integer.'
            sys.exit()
        if ntsmax <= header.n_class_steps:
            print 'Answer must be greater than or equal to the previous number of classical steps used, which was %d.\nTo reduce number of classical steps past %d, simply change n_class_steps in header.' % (header.n_class_steps,header.n_class_steps)
            sys.exit()
        nheader = open('%s/nheader' % (NEXMDir),'w')
        for line in header:
            if 'n_class_steps' in line:
                nheader.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (ntsmax))
            else:
                nheader.write(line)
        nheader.close()
        header.n_class_steps = ntsmax
        os.rename('%s/nheader' % (NEXMDir), '%s/header' % (NEXMDir))
    
    ## Choose random seeds ##
    if dynq == 0: ## ensemble
        ntraj = 0
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            data = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(data,int) == True:
                data = np.array([data])
            ntraj += len(data)
    if dynq == 1: ## single trajectory
        ntraj = 1
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
        if isinstance(rseeds,int) == True:
            rseeds = np.array([rseeds])
        lenrseeds = len(rseeds)
        if lenrseeds < ntraj:
            print 'Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories is %d.' % (lenrseeds,ntraj)
            sys.exit()
        for rseed in rseeds:
            if rseed < 0:
                print 'A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (rseed)
                sys.exit()

    ## Prepare NEXMD restart input files ##
    if dynq == 0:
        print 'Searching for incomplete trajectories and preparing input files.'
        #header = open('%s/header' % (NEXMDir),'r')
        #header = header.readlines()
        rseedslist = open('%s/rseedslist%d' % (NEXMDir,maxdir + 1),'w')
        error = open('%s/restart.err' % (cwd),'w')
        rtimes = np.array([])
        rstflag = 0
        traj = 0
        for NEXMD in NEXMDs:
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1,int) == True:
                dirlist1 = np.array([dirlist1])
            dirlist = open('%s/dirlist' % (NEXMD),'w')
            for dir in dirlist1:
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir)
                    print >> rseedslist, '%d' % (-123456789)
                    rstflag = 1
                    traj += 1
                    continue
                data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                data = data.readlines()
                tsteps = len(data) - 2
                if tsteps != header.n_class_steps:
                    if not os.path.exists('%s/%04d/restart.out' % (NEXMD,dir)):
                        print >> error, 'Path %s/%04d/restart.out does not exist.' % (NEXMD,dir)
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
                        if 'State' in line:
                            state = np.int(line.split()[-1])
                        if 'Seed' in line:
                            rseed = np.int(line.split()[-1])
                        if '$COORD' in line:
                            array = np.append(array,index)
                        if '$ENDCOORD' in line:
                            array = np.append(array,index)
                        if '$VELOC' in line:
                            array = np.append(array,index)
                        if '$ENDVELOC' in line:
                            array = np.append(array,index)
                        if '$COEFF' in line:
                            array = np.append(array,index)
                        if '$ENDCOEFF' in line:
                            array = np.append(array,index)
                        index += 1
                    array = np.int_(array)
                    if len(array) != 6:
                        print >> error, 'Path %s%04d/restart.out is incomplete.' % (NEXMD,dir)
                        print >> rseedslist, '%d' % (-123456789)
                        rstflag = 1
                        traj += 1
                        continue
                    coords = data[array[0]:array[1]+1:1]
                    velocs = data[array[2]:array[3]+1:1]
                    coeffs = data[array[4]:array[5]+1:1]
                    ## Start renormalize coefficients ##
                    ncoeffs = np.zeros(((len(coeffs) - 2), 2))
                    index = 0
                    for line in coeffs[1:len(coeffs) - 1:1]:
                        ncoeffs[index] = np.float_(line.split())
                        index += 1
                    if np.sum(ncoeffs[:,0]) != 0:
                        ncoeffs[:,0] = ncoeffs[:,0]/np.sum(ncoeffs[:,0])
                    ## Find the maximum view file ##
                    data = glob.glob('%s/%04d/view*' % (NEXMD,dir))
                    data = [ x[:-10] for x in data ]
                    if len(data) != 0:
                        maxview = np.int(extract(max(data,key = extract))[0])
                    else:
                        maxview = 0
                    ## Make new input file ##
                    inputfile = open('%s/%04d/input.ceon' % (NEXMD,dir),'w')
                    for line in header:
                        if 'rnd_seed' in line:
                            inputfile.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseed if randq == 1 else rseeds[traj]))
                        else:
                            if 'exc_state_init_flag' in line:
                                inputfile.write('   exc_state_init=%d, ! initial excited state (0 - ground state) [0]\n' % (state))
                            else:
                                if 'time_init' in line:
                                    inputfile.write('   time_init=%.1f, ! initial time, fs [0.00]\n' % (time))
                                else:
                                    if 'n_class_steps' in line:
                                        inputfile.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (np.int(header.n_class_steps-time/header.time_step)))
                                    else:
                                        if 'out_count_init' in line:
                                            inputfile.write('   out_count_init=%d, ! initial count for output files [0]\n' % (maxview))
                                        else:
                                            if 'nucl_coord_veloc' in line:
                                                for line in coords:
                                                    inputfile.write(line)
                                                inputfile.write('\n')
                                                for line in velocs:
                                                    inputfile.write(line)
                                            else:
                                                if 'quant_amp_phase' in line:
                                                    inputfile.write('$coeff\n')
                                                    for line in ncoeffs:
                                                        inputfile.write('  %.10f  %.10f\n' % (line[0],line[1]))
                                                    inputfile.write('$endcoeff\n')
                                                else:
                                                    inputfile.write(line)
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

    if dynq == 1: ## single trajectory
        print 'Searching for incomplete trajectory and preparing input file.'
        #header = open('%s/header' % (NEXMDir),'r')
        #header = header.readlines()
        rseedslist = open('%s/rseedslist%d' % (NEXMDir,maxdir + 1),'w')
        rtimes = np.array([])
        rstflag = 0
        traj = 0
        if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
            print 'Path %s/energy-ev.out does not exist.' % (NEXMDir)
            os.remove('%s/rseedslist%d' % (NEXMDir, maxdir + 1))
            sys.exit()
        data = open('%s/energy-ev.out' % (NEXMDir),'r')
        data = data.readlines()
        tsteps = len(data) - 2
        if tsteps != header.n_class_steps:
            if not os.path.exists('%s/restart.out' % (NEXMDir)):
                print >> error, 'Path %s/restart.out does not exist.' % (NEXMDir)
                os.remove('%s/rseedslist%d' % (NEXMDir, maxdir + 1))
                sys.exit()
            data = open('%s/restart.out' % (NEXMDir),'r')
            data = data.readlines()
            index = 0
            array = np.array([])
            for line in data:
                if 'time' in line:
                    time = np.around(np.float(line.split()[-1]), decimals = 3)
                if 'State' in line:
                    state = np.int(line.split()[-1])
                if 'Seed' in line:
                    rseed = np.int(line.split()[-1])
                if '$COORD' in line:
                    array = np.append(array,index)
                if '$ENDCOORD' in line:
                    array = np.append(array,index)
                if '$VELOC' in line:
                    array = np.append(array,index)
                if '$ENDVELOC' in line:
                    array = np.append(array,index)
                if '$COEFF' in line:
                    array = np.append(array,index)
                if '$ENDCOEFF' in line:
                    array = np.append(array,index)
                index += 1
            array = np.int_(array)
            if len(array) != 6:
                print 'Path %s/restart.out is incomplete.' % (NEXMDir)
                os.remove('%s/rseedslist%d' % (NEXMDir, maxdir + 1))
            coords = data[array[0]:array[1]+1:1]
            velocs = data[array[2]:array[3]+1:1]
            coeffs = data[array[4]:array[5]+1:1]
            ## Start renormalize coefficients ##
            ncoeffs = np.zeros(((len(coeffs) - 2), 2))
            index = 0
            for line in coeffs[1:len(coeffs) - 1:1]:
                ncoeffs[index] = np.float_(line.split())
                index += 1
            if np.sum(ncoeffs[:,0]) != 0:
                ncoeffs[:,0] = ncoeffs[:,0]/np.sum(ncoeffs[:,0])
            ## Find the maximum view file ##
            data = glob.glob('%s/view*' % (NEXMDir))
            data = [ x[:-10] for x in data ]
            if len(data) != 0:
                maxview = np.int(extract(max(data,key = extract))[0])
            else:
                maxview = 0
            ## Make new input file ##
            inputfile = open('%s/input.ceon' % (NEXMDir),'w')
            for line in header:
                if 'rnd_seed' in line:
                    inputfile.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseed if randq == 1 else rseeds[traj]))
                else:
                    if 'exc_state_init_flag' in line:
                        inputfile.write('   exc_state_init=%d, ! initial excited state (0 - ground state) [0]\n' % (state))
                    else:
                        if 'time_init' in line:
                            inputfile.write('   time_init=%.1f, ! initial time, fs [0.00]\n' % (time))
                        else:
                            if 'n_class_steps' in line:
                                inputfile.write('   n_class_steps=%d, ! number of classical steps [1]\n' % (np.int(header.n_class_steps-time/header.time_step)))
                            else:
                                if 'out_count_init' in line:
                                    inputfile.write('   out_count_init=%d, ! initial count for output files [0]\n' % (maxview))
                                else:
                                    if 'nucl_coord_veloc' in line:
                                        for line in coords:
                                            inputfile.write(line)
                                        inputfile.write('\n')
                                        for line in velocs:
                                            inputfile.write(line)
                                    else:
                                        if 'quant_amp_phase' in line:
                                            inputfile.write('&coeff\n')
                                            for line in ncoeffs:
                                                inputfile.write('  %.10f  %.10f\n' % (line[0],line[1]))
                                            inputfile.write('&endcoeff\n')
                                        else:
                                            inputfile.write(line)
            rtimes = np.append(rtimes,time)
            print '%s' % (NEXMDir)
            print >> rseedslist, '%d' % (rseed if randq == 1 else rseeds[traj])
        else:
            print 'Trajectory has completed.'
            sys.exit()
        traj += 1
        if filecmp.cmp('%s/rseedslist%d' % (NEXMDir,maxdir), '%s/rseedslist%d' % (NEXMDir,maxdir + 1)):
            os.remove('%s/rseedslist%d' % (NEXMDir,maxdir + 1))

    ## Determine whether or not to delete extraneous data in output files ##
    contq = input('Continue to delete extraneous data in output files? Answer yes [1] or no [0]: ')
    if contq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if contq == 0:
        sys.exit()
    print 'Deleting extraneous data in output files. please wait ...'
    if header.bo_dynamics_flag == 0:
        files = np.array(['coeff-n.out', 'energy-ev.out', 'nact.out', 'order.out', 'pes.out', 'temperature.out', 'transition-densities.out'])
        ofiles = np.array(['hops.out', 'hops-trial.out'])
    else:
        files = np.array(['energy-ev.out', 'temperature.out', 'pes.out','transition-densities.out'])
        ofiles = np.array([])

    ## Delete extraneous data ##
    if dynq == 0: ## ensemble
        error = open('%s/delextra.err' % (cwd),'w')
        rstflag = 0
        traj = 0
        for NEXMD in NEXMDs:
            if os.stat('%s/%s/dirlist' % (cwd,NEXMD)).st_size == 0:
                continue
            else:
                dirlist = np.int_(np.genfromtxt('%s/%s/dirlist' % (cwd,NEXMD)))
            if isinstance(dirlist,int) == True:
                    dirlist = np.array([dirlist])
            for dir in dirlist:
                if not os.path.exists('%s/%s/%04d' % (cwd,NEXMD,dir)):
                    print >> error, 'Path %s/%s%04d does not exist.' % (cwd,NEXMD,dir)
                    rstflag = 1
                    traj += 1
                    continue
                os.chdir('%s/%s/%04d' % (cwd,NEXMD,dir))
                for index in np.arange(len(files)):
                    if not os.path.exists('%s/%s/%04d/%s' % (cwd,NEXMD,dir,files[index])):
                        print >> error, 'Path %s/%s%04d/%s does not exist.' % (cwd,NEXMD,dir,files[index])
                        rstflag = 1
                        continue
                    ## Derivation of the following algorithm is provided at the end of this script ##
                    data = subprocess.check_output(['tail','-1','%s' % (files[index])])
                    ltime = np.around(np.float(np.fromstring(data,dtype=float,sep=' ')[1 if index == 0 and header.bo_dynamics_flag == 0 else 0]), decimals = 3)
                    if rtimes[traj] > ltime + header.out_data_steps*header.time_step:
                        print >> error, 'last time-step in', '%s/%s%04d/%s' % (cwd,NEXMD,dir,'restart.out'), 'exceeds last time-step in', '%s/%s%04d/%s' % (cwd,NEXMD,dir,files[index])
                        rstflag = 1
                        continue
                    ncsteps = 0
                    while rtimes[traj] + ncsteps*(header.out_data_steps*header.time_step) <= ltime:
                        ncsteps += 1
                    if header.moldyn_verbosity == 3 and index in [2,4]:
                        lctime = rtimes[traj] + (ncsteps - 1)*(header.out_data_steps*header.time_step)
                        nqsteps = 0
                        while lctime*header.n_quant_steps + nqsteps*header.time_step <= ltime*header.n_quant_steps:
                            nqsteps += 1
                        if ltime == lctime:
                            nlines = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps - 1) + 1) + (nqsteps - 1) + 1
                        else:
                            ncsteps = 0
                            while rtimes[traj] + ncsteps*header.time_step <= ltime:
                                ncsteps += 1
                            nlines = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps - 1) + 1) + (nqsteps - 1) - (ncsteps - 1) + 1
                    else:
                        nlines = (ncsteps - 1) + 1
                    if header.bo_dynamics_flag == 0:
                        subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack,(nlines - 1 if index in [2,3,5] else nlines),files[index])))
                    else:
                        subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack,(nlines - 1 if index == 1 else nlines),files[index])))
                    os.rename('%s.restart' % (files[index]), '%s' % (files[index]))
                for index in np.arange(len(ofiles)):
                    if not os.path.exists('%s/%s/%04d/%s' % (cwd,NEXMD,dir,ofiles[index])):
                        print >> error, 'Path %s/%s%04d/%s does not exist.' % (cwd,NEXMD,dir,ofiles[index])
                        rstflag = 1
                        continue
                    ltime = 1000000
                    nlines = 0
                    while ltime >= rtimes[traj]:
                        data = subprocess.check_output(['tail','%d' % (-(nlines + 1)),'%s' % (ofiles[index])])
                        ltime = np.float(np.fromstring(data,dtype=float,sep=' ')[0])
                        nlines += 1
                    subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack, nlines - 1, ofiles[index])))
                    os.rename('%s.restart' % (ofiles[index]), '%s' % (ofiles[index]))
                print '%s%04d' % (NEXMD,dir)
                traj += 1
        if rstflag == 1:
            print 'One or more trajectories cannot be restarted properly, check delextra.err.'
        else:
            os.remove('%s/delextra.err' % (cwd))
        print 'A total of %d trajectories have been prepared for restart.' % (traj)

    ## Delete extraneous data ##
    if dynq == 1: ## single trajectory
        error = open('%s/delextra.err' % (cwd),'w')
        rstflag = 0
        traj = 0
        os.chdir('%s/%s' % (cwd,NEXMDir))
        for index in np.arange(len(files)):
            if not os.path.exists('%s/%s/%s' % (cwd,NEXMDir,files[index])):
                print >> error, '%s/%s/%s does not exist.' % (cwd,NEXMDir,files[index])
                rstflag = 1
                continue
            ## Derivation of the following algorithm is provided at the end of this script ##
            data = subprocess.check_output(['tail','-1','%s' % (files[index])])
            ltime = np.around(np.float(np.fromstring(data,dtype=float,sep=' ')[1 if index == 0 and header.bo_dynamics_flag == 0 else 0]), decimals = 3)
            if rtimes[traj] > ltime + header.out_data_steps*header.time_step:
                print >> error, 'last time-step in', '%s/%s/%s' % (cwd,NEXMDir,'restart.out'), 'exceeds last time-step in', '%s/%s/%s' % (cwd,NEXMDir,files[index])
                rstflag = 1
                continue
            ncsteps = 0
            while rtimes[traj] + ncsteps*(header.out_data_steps*header.time_step) <= ltime:
                ncsteps += 1
            if header.moldyn_verbosity == 3 and index in [2,4]:
                lctime = rtimes[traj] + (ncsteps - 1)*(header.out_data_steps*header.time_step)
                nqsteps = 0
                while lctime*header.n_quant_steps + nqsteps*header.time_step <= ltime*header.n_quant_steps:
                    nqsteps += 1
                if ltime == lctime:
                    nlines = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps - 1) + 1) + (nqsteps - 1) + 1
                else:
                    ncsteps = 0
                    while rtimes[traj] + ncsteps*header.time_step <= ltime:
                        ncsteps += 1
                    nlines = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps - 1) + 1) + (nqsteps - 1) - (ncsteps - 1) + 1
            else:
                nlines = (ncsteps - 1) + 1
            if header.bo_dynamics_flag == 0:
                subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack,(nlines - 1 if index in [2,3,5] else nlines),files[index])))
            else:
                subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack,(nlines - 1 if index == 1 else nlines),files[index])))
            os.rename('%s.restart' % (files[index]), '%s' % (files[index]))
        for index in np.arange(len(ofiles)):
            if not os.path.exists('%s/%s/%s' % (cwd,NEXMDir,ofiles[index])):
                print >> error, '%s/%s/%s' % (cwd,NEXMDir,ofiles[index]), 'does not exist'
                rstflag = 1
                continue
            ltime = 1000000
            nlines = 0
            while ltime >= rtimes[traj]:
                data = subprocess.check_output(['tail','%d' % (-(nlines + 1)),'%s' % (ofiles[index])])
                ltime = np.float(np.fromstring(data,dtype=float,sep=' ')[0])
                nlines += 1
            subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (pathtopack, nlines - 1, ofiles[index])))
            os.rename('%s.restart' % (ofiles[index]), '%s' % (ofiles[index]))
        print '%s' % (NEXMDir)
        traj += 1
        if rstflag == 1:
            print 'Trajectory cannot be restarted properly, check delextra.err.'
        else:
            os.remove('%s/delextra.err' % (cwd))
        print 'A total of %d trajectory has been prepared for restart.' % (traj)

'''
    
    header.time_step     = time-step
    header.moldyn_verbosity   = verbosity of molecular dynamics data
    header.out_data_steps  = data are printed every out_data_steps classical time-steps
    ltime  = last printed time-step (may be classical or quantum)
    lctime = last printed classical time-step
    rtime  = restart time-step (always classical)
    header.n_quant_steps= number of quantum steps per classical step
    
    a = (lctime - rtime)/header.time_step                                                                                   = total # of classical steps computed
    b = (lctime - rtime)/(header.out_data_steps x header.time_step)                                                         = total # of classical steps printed
    c = (lctime - rtime)/header.time_step x ((header.out_data_steps - 1)/header.out_data_steps)                             = total # of classical steps not printed
    d = (header.n_quant_steps- 1) x (lctime - rtime)/header.time_step x ((header.out_data_steps - 1)/header.out_data_steps) = total # of quantum steps from unprinted classical steps
    e = header.n_quant_stepsx (lctime - rtime)/(header.out_data_steps x header.time_step)                                   = total # of classical and quantum steps from printed classical steps
    
    if header.moldyn_verbosity != 3:
        total # of lines to delete = b
    if header.moldyn_verbosity == 3 and file not in [nact.out, pes.out]:
        total # of lines to delete = b
    if header.moldyn_verbosity == 3 and file in [nact.out, pes.out]:
        total # of lines to delete = d + e + resdiual quantum steps after lctime
        if ltime == lctime:
            residual quantum steps after lctime = header.n_quant_stepsx (ltime - lctime)/header.time_step
            total # of lines to delete = d + e + header.n_quant_stepsx (ltime - lctime)/header.time_step + 1 ( + 1 to delete the restart time-step )
            total # of lines to delete = (lctime - rtime)*(header.out_data_steps*(header.n_quant_steps- 1) + 1)/(header.out_data_steps*header.time_step) + header.n_quant_steps*(ltime - lctime)/header.time_step + 1
            total # of lines to delete = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps- 1) + 1) + (nqsteps - 1) + 1
            *** see in code how to calculate ncsteps and nqsteps ***
        if ltime != lctime:
            residual quantum steps after lctime = header.n_quant_stepsx (ltime - lctime)/header.time_step
            number of unprinted classical steps after lctime = ncsteps - 1
            total # of lines to delete = d + e + header.n_quant_stepsx (ltime - lctime)/header.time_step - (ncsteps - 1) + 1 ( + 1 to delete the restart time-step)
            total # of lines to delete = (lctime - rtime)*(header.out_data_steps*(header.n_quant_steps- 1) + 1)/(header.out_data_steps*header.time_step) + header.n_quant_steps*(ltime - lctime)/header.time_step - (ncsteps - 1) + 1
            total # of lines to delete = (ncsteps - 1)*(header.out_data_steps*(header.n_quant_steps- 1) + 1) + (nqsteps - 1) - (ncsteps - 1) + 1
            *** see in code how to calculate ncsteps and nqsteps ***
            
'''
