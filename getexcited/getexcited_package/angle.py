#/usr/bin/python

'''

This function calculates a user-specified angle in time.

The angle is calculated between three atoms at every time-step
and may be tracked along a single trajectory or an ensemble of
trajectories during adiabatic or non-adiabatic dynamics.  The
user must supply the line numbers of three atoms, where the first atom
between '&coord' and '&endcoord' in 'input.ceon' is always labeled with
line number = 0.  The second atom specified is the vertex.

Type of calculation:

[1] Single Trajectory
< collection time (fs)
> time (fs), angle (degrees)

[2] Ensemble of Trajectories

[2a] Mean
< collection time (fs)
> time (fs), angle (degrees), standard deviation (degrees)

[2b] All time-steps
> trajectory directory, angle (degrees) at all time-steps and
trajectories

Output Files:
- angle_[type].out, where [type] = single, mean_ensemble, raw_ensemble

Error Files:
- angle_[type].err, where [type] = single, mean_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob
import fileinput
import math

cwd = os.getcwd()

def angle(header):

    print('Calculating angle between bonds as a function of time.')

    ## Type of calculation and directory check ##
    dynq = input('Calculate angle between bonds along one trajectory or an ensemble of trajectories?\nAnswer one [1] or ensemble [0]: ')
    if dynq not in [1,0]:
        print('Answer must be 1 or 0.')
        sys.exit()
    if dynq == 0: ## ensemble
        NEXMDir = raw_input('Ensemble directory [e.g. NEXMD]: ')
        if not os.path.exists(NEXMDir):
            print('Path %s does not exist.' % (NEXMDir))
            sys.exit()
        ## Check if NEXMD folders exist ##
        NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
        NEXMDs.sort()
        if len(NEXMDs) == 0:
            print('There are no NEXMD folders in %s.' % (NEXMDir))
            sys.exit()
        ## Determine mean or all ##
        typeq = input('Output mean angle in time or output angle at all time-steps and trajectories?\nAnswer mean [0] or all [1]: ')
        if typeq not in [0,1]:
            print('Answer must be 0 or 1.')
            sys.exit()
    if dynq == 1: ## single trajectory
        typeq = 0
        NEXMDir = raw_input('Single trajectory directory: ')
        if not os.path.exists(NEXMDir):
            print('Path %s does not exist.' % (NEXMDir))
            sys.exit()

    ## Information from header ##
    if dynq == 0: ## ensemble
        if not os.path.exists('%s/header' % (NEXMDir)):
            print('Path %s/header does not exist.' % (NEXMDir))
            sys.exit()
        header = header('%s/header' % (NEXMDir))
    if dynq == 1: ## single trajectory
        if not os.path.exists('%s/input.ceon' % (NEXMDir)):
            print('Path %s/input.ceon does not exist.' % (NEXMDir))
            sys.exit()
        header = header('%s/input.ceon' % (NEXMDir))

    ## Adding + 1 to include zeroth time-step ##
    header.n_class_steps = header.n_class_steps + 1
    
    ## Collection time ##
    if typeq == 0: ## mean angle
        if dynq == 0: ## ensemble
            tcoll = input('Calculate angle up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: ')
        if dynq == 1: ## single trajectory
            tcoll = input('Calculate angle up to what time in femtoseconds? ')
        if isinstance(tcoll, int) == False and isinstance(tcoll, float) == False:
            print('Time must be integer or float.')
            sys.exit()
        if tcoll < 0:
            print('Time must be integer or float greater than zero.')
            sys.exit()
        tcoll = np.float(tcoll)
        if tcoll > (header.n_class_steps - 1)*header.time_step:
            tcoll = (header.n_class_steps - 1)*header.time_step
    if typeq == 1: ## all bla
        tcoll = (header.n_class_steps - 1)*header.time_step

    ## Number of classical steps ##
    tscol = 0
    while tscol*header.time_step*header.out_data_steps <= tcoll:
        tscol += 1

    ## Number of time-steps for coordinates ##
    ccoll = 0
    num = 0
    while ccoll <= tcoll:
        ccoll += header.time_step*header.out_data_steps*header.out_coords_steps
        num += 1

    ## Collection time ##
    times = np.linspace(header.time_init, ccoll - header.time_step*header.out_data_steps*header.out_coords_steps, num)

    ## Three unique atoms defined by user ##
    lines = input('Input the line numbers labeling the three atoms in the form A(atom)-V(vertex atom)-A(atom).\nInput an array of the form [[atom1, atom2, atom3], .. ]: ')
    for line in lines:
        if isinstance(line, list) == False:
            print('Subarray must be of the form [atom1, atom2, atom3], where atom# = line number of atom#.')
            sys.exit()
        if len(line) != 3:
            print('Subarray must contain three elements labeling the line numbers of three atoms.')
            sys.exit()
        index = 0
        for atom in line:
            if isinstance(atom, int) == False:
                print('Element number %d of subarray must be integer.\nUser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1]))
                sys.exit()
            if atom < 0:
                print('Element number %d of subarray must be a positive integer.\nUser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1]))
                sys.exit()
            if atom > header.natoms - 1: # -1 for python indexing
                print('Element number %d of subarray must be less than the max number of atoms (-1).\nuser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1]))
                sys.exit()
            index += 1
        if len(np.unique(line)) != 3:
            print('All elements of subarray must be unique.\nUser inputted [%s, %s], which is not allowed.' % (line[0], line[1]))
            sys.exit()
    nangles = len(lines)

    ## Calculate angle along a single trajectory ##
    if dynq == 1: ## single trajectory
        print('Collecting angle along single trajectory.  Please wait ...')
        ## genrate output file ##
        output = open('%s/angle_single.out' % (cwd),'w')
        etraj = 0
        ## Determine completed number of time-steps ##
        if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
            print('Path %s/energy-ev.out does not exist.' % (NEXMDir))
            sys.exit()
        tsteps = np.genfromtxt('%s/energy-ev.out' % (NEXMDir), usecols=[0], skip_header=1).size
        ## Generate array with indices of the coordinate blocks along a trajectory ##
        if tsteps >= tscol:
            if not os.path.exists('%s/coords.xyz' % (NEXMDir)):
                print('Path %s/coords.xyz does not exist.' % (NEXMDir))
                sys.exit()
            data = open('%s/coords.xyz' % (NEXMDir),'r')
            data = data.readlines()
            lenc = len(data)
            ncoords = 0
            cindex = 0
            tflag1 = 0
            tflag2 = 0
            array = np.array([])
            for line in data:
                if 'time' in line:
                    if ncoords == 0:
                        tinit = np.float(line.split()[-1])
                        if tinit != header.time_init:
                            tflag1 = 1
                            break
                    else:
                        time = np.around(np.float(line.split()[-1]), decimals = 3)
                        if time > tcoll:
                            tflag3 = 1
                            break
                        if time != times[ncoords]:
                            tflag2 = 1
                            break
                    ncoords += 1
                    array = np.append(array,cindex)
                cindex += 1
            if tflag1 == 1:
                print('Initial time in %s/coords.xyz does not match time_init in %s/input.ceon.' % (NEXMDir,NEXMDir))
                sys.exit()
            if tflag2 == 1:
                print('There is an inconsistency in time-step in %s/coords.xyz.' % (NEXMDir))
                sys.exit()
            ## Append lines for last coordinate set ##
            if tflag3 == 1:
                array = np.append(array,cindex)
            else:
                array = np.append(array, lenc + 1)
            array = np.int_(array)
            ## Checks to ensure angle calculation ##
            if ncoords == 0:
                print('No coordinates were found in %s/coords.xyz' % (NEXMDir))
                sys.exit()
            if ncoords == 1:
                print('Only initial coordinates, at %.2f fs, were found in %s/coords.xyz.' % (tinit,NEXMDir))
                sys.exit()
            ## Calculate angle along a single trajectory ##
            sangle = np.zeros((ncoords,nangles))
            for ncoord in np.arange(ncoords):
                coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                index = 0
                for line in lines:
                    vec0 = np.float_(coords[line[0,0]].split()[1:])
                    vec1 = np.float_(coords[line[0,1]].split()[1:])
                    vec2 = np.float_(coords[line[0,2]].split()[1:])
                    a = np.subtract(vec0, vec1)
                    b = np.subtract(vec2, vec1)
                    cosine = np.arccos(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)))
                    sangle[ncoord,index] = cosine
                    index += 1
            print('%s' % (NEXMDir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
            ctraj = 1
            if tsteps == header.n_class_steps:
                etraj = 1
        else:
            print('%s' % (NEXMDir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
        ttraj = 1
        ## Summary of results ##
        if ctraj == 0:
            print('No trajectories completed within %0*.2f.' % (len(str(header.n_class_steps)),tcoll))
        else:
            print('Total trajectories:', '%04d' % (ttraj))
            print('Completed trajectories:', '%04d' % (ctraj))
            print('Excellent trajectories:', '%04d' % (etraj))
            output.wirte( 'Total trajectories: ', '%04d' % (ttraj))
            output.wirte( 'Completed trajectories: ', '%04d' % (ctraj))
            output.wirte( 'Excellent trajectories: ', '%04d' % (etraj))
            for ncoord in np.arange(ncoords):
                output.wirte( '%0*.2f' % (len(str((header.n_class_steps))) + 2,header.time_step*header.out_data_steps*header.out_coords_steps*ncoord), ' '.join('%08.3f' % (np.degrees(angle)) for angle in sangle[ncoord]))

    ## Calculate angle along an ensemble of trajectories ##
    if dynq == 0 and typeq == 0: ## mean from ensemble
        print('Collecting mean angle from ensemble.  Please wait ...')
        ## Determine total number of trajectories in ensemble ##
        with open('%s/totdirlist' % (NEXMDir),'w') as data:
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('Path %sdirlist1 does not exist.' % (NEXMD))
                    sys.exit()
                inputdata = fileinput.input('%s/dirlist1' % (NEXMD))
                data.writelines(inputdata)
        dirlist1 = np.int_(np.genfromtxt('%s/totdirlist' % (NEXMDir)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        os.remove('%s/totdirlist' % (NEXMDir))
        ## Generate output and error files ##
        output = open('%s/angle_mean_ensemble.out' % (cwd),'w')
        error = open('%s/angle_mean_ensemble.err' % (cwd),'w')
        ## Generate angle array for final results ##
        fangle = np.zeros((len(times), nangles))
        eangle = np.zeros((len(times), len(dirlist1), nangles))
        ttraj = 0
        ctraj = 0
        etraj = 0
        errflag = 0
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print('Path %dirlist1 does not exist.' % (NEXMD))
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1, int) == True:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine completed number of time-steps ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                    error.wirte( 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                ## Generate array with indices of the coordinate blocks along trajectory ##
                if tsteps >= tscol:
                    if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,dir)):
                        error.wirte( 'Path %s%04d/coords.xyz does not exist.' % (NEXMD,dir))
                        errflag = 1
                        ttraj += 1
                        continue
                    data = open('%s/%04d/coords.xyz' % (NEXMD,dir),'r')
                    data = data.readlines()
                    lenc = len(data)
                    ncoords = 0
                    cindex = 0
                    tflag1 = 0
                    tflag2 = 0
                    tflag3 = 0
                    array = np.array([])
                    for line in data:
                        if 'time' in line:
                            if ncoords == 0:
                                tinit = np.float(line.split()[-1])
                                if tinit != header.time_init:
                                    tflag1 = 1
                                    continue
                            else:
                                time = np.around(np.float(line.split()[-1]), decimals = 3)
                                if time > tcoll:
                                    tflag3 = 1
                                    continue
                                if time != times[ncoords]:
                                    tflag2 = 1
                                    continue
                            ncoords += 1
                            array = np.append(array,cindex)
                        cindex += 1
                    if tflag1 == 1:
                        error.wirte( 'Initial time in %s%04d/coords.xyz does not match time_init in %s/header.' % (NEXMD,dir,NEXMDir))
                        errflag = 1
                        ttraj += 1
                        continue
                    if tflag2 == 1:
                        error.wirte( 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (NEXMD,dir))
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Append lines for last coordinate set ##
                    if tflag3 == 1:
                        array = np.append(array,cindex)
                    else:
                        array = np.append(array, lenc + 1)
                    array = np.int_(array)
                    ## Checks to ensure angle calculation ##
                    if ncoords == 0:
                        error.wirte( 'No coordinates were found in %s%04d/coords.xyz' % (NEXMD,dir))
                        errflag = 1
                        ttraj += 1
                        continue
                    if ncoords == 1:
                        error.wirte( 'Only initial coordinates, at %.2f fs, were found in %s%04d/coords.xyz.' % (tinit,NEXMD,dir))
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Calculate angle along a single trajectory ##
                    sangle = np.zeros((ncoords, nangles))
                    for ncoord in np.arange(ncoords):
                        coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                        index = 0
                        for line in lines:
                            vec0 = np.float_(coords[line[0]].split()[1:])
                            vec1 = np.float_(coords[line[1]].split()[1:])
                            vec2 = np.float_(coords[line[2]].split()[1:])
                            a = np.subtract(vec0, vec1)
                            b = np.subtract(vec2, vec1)
                            cosine = np.arccos(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)))
                            sangle[ncoord,index] = cosine
                            eangle[ncoord,ctraj,index] = cosine
                            index += 1
                    fangle += sangle
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ctraj += 1
                    if tsteps == header.n_class_steps:
                        etraj += 1
                else:
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    error.wirte( '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    errflag = 1
                ttraj += 1
        ## Summary of results ##
        if ctraj == 0:
            print('No trajectories completed within %0*.2f.' % (len(str(header.n_class_steps)),tcoll))
        else:
            ## Mean and standard deviation for angle ##
            eangle = np.delete(eangle, np.arange(ctraj, ttraj), axis = 1)
            eangle = np.std(eangle, axis = 1)
            fangle = fangle/ctraj
            print('Total trajectories:', '%04d' % (ttraj))
            print('Completed trajectories:', '%04d' % (ctraj))
            print('Excellent trajectories:', '%04d' % (etraj))
            output.wirte( 'Total trajectories: ', '%04d' % (ttraj))
            output.wirte( 'Completed trajectories: ', '%04d' % (ctraj))
            output.wirte( 'Excellent trajectories: ', '%04d' % (etraj))
            for ncoord in np.arange(ncoords):
                output.wirte( '%0*.2f' % (len(str((header.n_class_steps))) + 2,header.time_step*header.out_data_steps*header.out_coords_steps*ncoord), ' '.join('%08.3f' % (np.degrees(angle)) for angle in fangle[ncoord]), ' '.join('%07.3f' % (np.degrees(angle)) for angle in eangle[ncoord]))
        if errflag == 1:
            print('One or more trajectories have experienced an error, check angle_mean_ensemble.err.')
        else:
            os.remove('%s/angle_mean_ensemble.err' % (cwd))

    ## Calculate angle from ensemble of trajectories at all time-steps ##
    if dynq == 0 and typeq == 1: ## all from ensemble
        print('Collecting all angles from ensemble.  Please wait ...')
        ## Generate output and error files ##
        output = open('%s/angle_raw_ensemble.out' % (cwd),'w')
        error = open('%s/angle_raw_ensemble.err' % (cwd),'w')
        ttraj = 0
        etraj = 0
        errflag = 0
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print('Path %dirlist1 does not exist.' % (NEXMD))
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1, int) == True:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine completed number of time-steps ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                    error.wirte( 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                ## Generate array with indices of the coordinate blocks along trajectory ##
                if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,dir)):
                    error.wirte( 'Path %s%04d/coords.xyz does not exist.' % (NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/coords.xyz' % (NEXMD,dir),'r')
                data = data.readlines()
                lenc = len(data)
                ncoords = 0
                cindex = 0
                tflag1 = 0
                tflag2 = 0
                tflag3 = 0
                array = np.array([])
                for line in data:
                    if 'time' in line:
                        if ncoords == 0:
                            tinit = np.float(line.split()[-1])
                            if tinit != header.time_init:
                                tflag1 = 1
                                continue
                        else:
                            time = np.around(np.float(line.split()[-1]), decimals = 3)
                            if time > tcoll:
                                tflag3 = 1
                                continue
                            if time != times[ncoords]:
                                tflag2 = 1
                                continue
                        ncoords += 1
                        array = np.append(array,cindex)
                    cindex += 1
                if tflag1 == 1:
                    error.wirte( 'Initial time in %s%04d/coords.xyz does not match time_init in %s/header.' % (NEXMD,dir,NEXMDir))
                    errflag = 1
                    ttraj += 1
                    continue
                if tflag2 == 1:
                    error.wirte( 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                ## Append lines for last coordinate set ##
                if tflag3 == 1:
                    array = np.append(array,cindex)
                else:
                    array = np.append(array, lenc + 1)
                array = np.int_(array)
                ## Checks to ensure angle calculation ##
                if ncoords == 0:
                    error.wirte( 'No coordinates were found in %s%04d/coords.xyz' % (NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                if ncoords == 1:
                    error.wirte( 'Only initial coordinates, at %.2f fs, were found in %s%04d/coords.xyz.' % (tinit,NEXMD,dir))
                    errflag = 1
                    ttraj += 1
                    continue
                ## Calculate angle along a single trajectory ##
                for ncoord in np.arange(ncoords):
                    coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                    sangle = np.zeros(nangles)
                    index = 0
                    for line in lines:
                        vec0 = np.float_(coords[line[0]].split()[1:])
                        vec1 = np.float_(coords[line[1]].split()[1:])
                        vec2 = np.float_(coords[line[2]].split()[1:])
                        a = np.subtract(vec0, vec1)
                        b = np.subtract(vec2, vec1)
                        cosine = np.arccos(np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b)))
                        sangle[index] = cosine
                        index += 1
                    output.wirte( '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2,header.time_step*header.out_data_steps*header.out_coords_steps*ncoord), ' '.join('%08.3f' % (np.degrees(angle)) for angle in sangle))
                print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                if tsteps == header.n_class_steps:
                    etraj += 1
                else:
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    error.wirte( '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    errflag = 1
                ttraj += 1
        ## Summary of results ##
        if ttraj == 0:
            print('No trajectories completed within %0*.2f.' % (len(str(header.n_class_steps)),tcoll))
        else:
            print('Total trajectories:', '%04d' % (ttraj))
            print('Excellent trajectories:', '%04d' % (etraj))
        if errflag == 1:
            print('One or more trajectories have experienced an error, check angle_raw_ensemble.err.')
        else:
            os.remove('%s/angle_raw_ensemble.err' % (cwd))
