#/usr/bin/python

'''

This function calculates excited-state permanent dipole
moment in time.

The dipole moment at every time-step is determined by the
dipole moment of the occupied state according to the
adiabatic or nonadiabatic dynamics and may be tracked
along a single trajectory or an ensemble of trajectories.

Type of calculation:

[1] Single Trajectory
< collection time (fs)
> time (fs), dipole moment (Debye)

[2] Ensemble of Trajectories

[2a] Mean
< collection time (fs)
> time (fs), dipole moment (Debye), standard deviation (Debye)

[2b] All time-steps
> trajectory directory, dipole moment (Debye) at all times and 
trajectories

Other options include:

[1] Relative direction of dipole
< line numbers of two atoms in the system, as presented
in 'input.ceon', which will be used to construct a vector and 
dotted with the dipole moment
> for option [1]: angle (degrees)
> for option [2a]: angle (degrees), standard deviation (degrees)
> for options [2b]: angle (degrees)

Output Files:
- permdipole_[type].out, where [type] = single, mean_ensemble, raw_ensemble

Error Files:
- permdipole_[type].err, where [type] = single, mean_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob
import subprocess
import shlex
import fileinput

cwd = os.getcwd()

def permdipole(pathtodipole):

    print 'Calculating excited-state permanent dipole moment as a function of time.'

    ## Type of calculation and directory check ##
    dynq = input('Calculate dipole moment along one trajectory or an ensemble of trajectories?\nAnswer one [1] or ensemble [0]: ')
    if dynq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if dynq == 0: ## ensemble
        NEXMDir = raw_input('Ensemble directory [e.g. NEXMD]: ')
        if not os.path.exists(NEXMDir):
            print 'Path %s does not exist.' % (NEXMDir)
            sys.exit()
        ## check if NEXMD folders exist ##
        NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
        NEXMDs.sort()
        if len(NEXMDs) == 0:
            print 'There are no NEXMD folders in %s.' % (NEXMDir)
            sys.exit()
        ## Determine mean or all ##
        typeq = input('Output mean dipole in time or output dipoles at all time-steps and trajectories?\nAnswer mean [0] or all [1]: ')
        if typeq not in [0,1]:
            print 'Answer must be 0 or 1.'
            sys.exit()
    if dynq == 1: ## single trajectory
        typeq = 0
        NEXMDir = raw_input('Single trajectory directory: ')
        if not os.path.exists(NEXMDir):
            print 'Path %s does not exist.' % (NEXMDir)
            sys.exit()

    ## User-defined length of analysis and initialize arrays ##
    if dynq == 0: ## ensemble
        if not os.path.exists('%s/header' % (NEXMDir)):
            print 'Path %s/header does not exist.' % (NEXMDir)
            sys.exit()
        header = open('%s/header' % (NEXMDir),'r')
        header = header.readlines()
    if dynq == 1: ## single trajectory
        if not os.path.exists('%s/input.ceon' % (NEXMDir)):
            print 'Path %s/input.ceon does not exist.' % (NEXMDir)
            sys.exit()
        header = open('%s/input.ceon' % (NEXMDir),'r')
        header = header.readlines()
    stateinit = none
    for line in header:
        if 'bo_dynamics_flag' in line:
            boflag = np.int(line.split()[0][len('bo_dynamics_flag='):-1])
        if 'exc_state_init=' in line:
            stateinit = np.int(line.split()[0][len('exc_state_init='):-1])
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
        if 'out_coords_steps' in line:
            cdata = np.int(line.split()[0][len('out_coords_steps='):-1])
        if 'natoms' in line:
            natoms = np.int(line.split()[0][len('natoms='):-1])
    if boflag == 1 and stateinit == none:
        print 'Dynamics are set to born-oppenheimer, but the initial state is not set.\nPlease check bo_dynamics_flag and exc_state_init in header.'
        sys.exit()

    ## Collection time ##
    if typeq == 0: ## mean dipole
        if dynq == 0: ## ensemble
            tcoll = input('Calculate dipole up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: ')
        if dynq == 1: ## single trajectory
            tcoll = input('Calculate dipole up to what time in femtoseconds? ')
        if isinstance(tcoll, int) == false and isinstance(tcoll, float) == false:
            print 'Time must be integer or float.'
            sys.exit()
        if tcoll < 0:
            print 'Time must be integer or float greater than zero.'
            sys.exit()
        tcoll = np.float(tcoll)
        if tcoll > (tsmax - 1)*dt:
            tcoll = (tsmax - 1)*dt
    if typeq == 1: ## all dipoles
        tcoll = (tsmax - 1)*dt

    ## Determine direction of dipole ##
    dotq = input('Find the angle between the dipole and a user-defined vector on the molecule?\nAnswer yes [1] or no [0]: ')
    if dotq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if dotq == 0:
        if natoms < 2:
            print 'Number of atoms set under natoms is less than zero.\nPlease check header (for ensemble) or input.ceon (for single trajectory).'
            sys.exit()
        else:
            lines = [0,1]
    if dotq == 1:
        lines = input('Input an array of the form [atom1, atom2], where atom# = line number of atom (0 is the first line).\nThese two atoms will be used to construct a vector: ')
        if isinstance(lines, list) == false:
            print 'Input must be an array of the form [atom 1, atom2], where atom# = line number of atom (0 is the first line).'
            sys.exit()
        if len(lines) != 2:
            print 'Input must be an array with two elements labeling the line numbers of two atoms.'
            sys.exit()
        index = 0
        for i in lines:
            if isinstance(i, int) == false:
                print 'Element number %d of input array must be integer.\nUser inputted [%s, %s], which is not allowed.' % (index + 1, lines[0], lines[1])
                sys.exit()
            if i < 0:
                print 'Element number %d of input array must be a positive integer.\nUser inputted [%s, %s], which is not allowed.' % (index + 1, lines[0], lines[1])
                sys.exit()
            if i > natoms - 1:
                print 'Element number %d of input array must be less than the max number of atoms (-1).\nUser inputted [%s, %s], which is not allowed.' % (index + 1, lines[0], lines[1])
                sys.exit()
            index += 1
        if len(np.unique(lines)) != 2:
            print 'All elements of input array must be unique.\nUser inputted [%s, %s], which is not allowed.' % (lines[0], lines[1])
            sys.exit()

    ## Number of classical time-steps ##
    tscol = 0
    while tscol*dt*odata <= tcoll:
        tscol += 1
    ccoll = 0
    
    ## Number of time-steps for coordinates ##
    num = 0
    while ccoll <= tcoll:
        ccoll += dt*odata*cdata
        num += 1
    
    ## Number of dipoles ##
    edipoles = ccoll

    ## Collection time array ##
    times = np.linspace(tinith, ccoll - dt*odata*cdata, num)

    ## Grep permanent dipole moments and states from ensembles ##
    if dynq == 0: ## ensemble
        print 'Checking permanent dipole moments and states.  Please wait ...'
        ## Checks to make sure scripts are available ##
        if not os.path.exists('%s/getexcited_package/collectdipline.sh' % (pathtodipole)):
            print 'The script, collectdipline.sh, must be in the getexcited_package.'
            sys.exit()
        if not os.path.exists('%s/getexcited_package/collectpermdipole.sh' % (pathtodipole)):
            print 'The script, collectpermdipole.sh, must be in the getexcited_package.'
            sys.exit()
        ## Generation of error file ##
        error = open('%s/dipole_collection_ensemble.err' % (cwd),'w')
        errflag = 0
        for NEXMD in NEXMDs:
            ## Check and open list of directories ##
            if not os.path.exists('%s/%s/dirlist1' % (cwd,NEXMD)):
                print 'Path %s/%s/dirlist1 does not exist.' % (cwd,NEXMD)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
            if isinstance(dirlist1,int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Check if directory exists ##
                if not os.path.exists('%s/%s/%04d' % (cwd,NEXMD,dir)):
                    print >> error, 'Path %s%04d does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                ## Go to directory ##
                os.chdir('%s/%s/%04d' % (cwd,NEXMD,dir))
                ## Check if standard output exists ##
                if not os.path.exists('%s/%s/%04d/md.out' % (cwd,NEXMD,dir)):
                    print >> error, 'Path %s%04d/md.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                ## Grep line number of classical step and classical step ##
                subprocess.call(shlex.split('sh %s/getexcited_package/collectdipline.sh %d' % (pathtodipole, odata*cdata)))
                if not os.path.exists('%s/%s/%04d/dipline.out' % (cwd,NEXMD,dir)):
                    print >> error, 'Path %s%04d/dipline.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                print '%s%04d' % (NEXMD,dir), 'dipole lines in md.out found'
                ## data = [line number of classical step, classical step] ##
                data = np.genfromtxt('%s/%s/%04d/dipline.out' % (cwd,NEXMD,dir))
                tdipoles = len(data)
                ## Check to ensure dipole calculation ##
                if np.array_equal(np.around(data[1:edipoles:1,1]*dt, decimals = 3), times[1:edipoles:1]) == false:
                    print >> error, 'There is an inconsistency in time-step in %s%04d/dipline.out.' % (NEXMD,dir)
                    errflag = 1
                    continue
                ## Delete previous permdipole.out if exists ##
                if os.path.exists('%s/%s/%04d/permdipole.out' % (cwd,NEXMD,dir)):
                    os.remove('%s/%s/%04d/permdipole.out' % (cwd,NEXMD,dir))
                ## Grep dipoles from standard output ##
                subprocess.call(shlex.split('sh %s/getexcited_package/collectpermdipole.sh %d' % (pathtodipole, nstates + 2)))
                if not os.path.exists('%s/%s/%04d/permdipole.out' % (cwd,NEXMD,dir)):
                    print >> error, 'Path %s%04d/permdipole.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                print '%s%04d' % (NEXMD,dir), 'dipoles in md.out extracted'
                ## Another check to ensure dipole calculation ##
                with open('%s/%s/%04d/permdipole.out' % (cwd,NEXMD,dir),'r') as data:
                    if len(data.readlines()) != tdipoles*(nstates + 3):
                        print >> error, 'Path %s%04d/permdipole.out is incomplete.' % (NEXMD,dir)
                        errflag = 1
                        os.remove('%s/%s/%04d/permdipole.out' % (cwd,NEXMD,dir))
                        continue
                ## Delete previous pop.out if exists ##
                if os.path.exists('%s/%s/%04d/pop.out' % (cwd,NEXMD,dir)):
                    os.remove('%s/%s/%04d/pop.out' % (cwd,NEXMD,dir))
                ## Get states ##
                if boflag == 1: ## adiabatic
                    states = np.int_(np.ones(edipoles)*statinit)
                if boflag == 0: ## nonadiabatic
                    ## Check coefficient file exists ##
                    if not os.path.exists('%s/%s/%04d/coeff-n.out' % (cwd,NEXMD,dir)):
                        print >> error, 'Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir)
                        errflag = 1
                        continue
                    data = open('%s/%s/%04d/coeff-n.out' % (cwd,NEXMD,dir),'r')
                    data = data.readlines()
                    states = np.zeros(edipoles)
                    index = 0
                    for line in data[0:tscol:odata*cdata]:
                        val = line.split()
                        pes = np.int(val[0])
                        time = np.around(np.float(val[1]), decimals = 3)
                        ## another check to ensure dipole calculation ##
                        if time != times[index]:
                            print >> error, 'There is an inconsistency in time-step in %s%04d/coeff-n.out at %.3f fs' % (NEXMD,dir,times[index])
                            errflag = 1
                            break
                        states[index] = pes
                        index += 1
                    ## if simulation becomes adiabatic after nonadiabatic ##
                    if len(data[0:tscol:odata*cdata]) < edipoles:
                        states[index::] = pes
                ## save populations to file ##
                np.savetxt('pop.out', np.transpose([times,states]), fmt=['%10.5e','%d'])
        if errflag == 1:
            print 'One or more trajectories have experienced an error, check dipole_collection_ensemble.err.'
            contq = input('Continue? Answer yes [1] or no [0]: ')
            if contq not in [1,0]:
                print 'Answer must to be 1 or 0.'
                sys.exit()
            if contq == 0:
                sys.exit()
        else:
            os.remove('%s/dipole_collection_ensemble.err' % (cwd))

    ## Grep permanent dipole moments and states from single trajectory ##
    if dynq == 1: ## single trajectory
        print 'Checking permanent dipole moments and states.  Please wait ...'
        ## Checks to make sure scripts are available ##
        if not os.path.exists('%s/getexcited_package/collectdipline.sh' % (pathtodipole)):
            print 'The script, collectdipline.sh, must be in the getexcited_package.'
            sys.exit()
        if not os.path.exists('%s/getexcited_package/collectpermdipole.sh' % (pathtodipole)):
            print 'The script, collectpermdipole.sh, must be in the getexcited_package.'
            sys.exit()
        ## Go to directory ##
        os.chdir('%s/%s' % (cwd,NEXMDir))
        ## Check if standard output exists ##
        if not os.path.exists('%s/%s/md.out' % (cwd,NEXMDir)):
            print 'Path %s/md.out does not exist.' % (NEXMDir)
            sys.exit()
        ## Grep line number of classical step and classical step ##
        subprocess.call(shlex.split('sh %s/getexcited_package/collectdipline.sh %d' % (pathtodipole, odata*cdata)))
        if not os.path.exists('%s/%s/dipline.out' % (cwd,NEXMDir)):
            print 'Path %s/dipline.out does not exist.' % (NEXMDir)
            sys.exit()
        print '%s' % (NEXMDir), 'dipole lines in md.out found'
        ## data = [line number of classical step, classical step] ##
        data = np.genfromtxt('%s/%s/dipline.out' % (cwd,NEXMDir))
        tdipoles = len(data)
        ## Check to ensure dipole calculation ##
        if np.array_equal(np.around(data[1:edipoles:1,1]*dt, decimals = 3), times[1:edipoles:1]) == false:
            print 'There is an inconsistency in time-step in %s/dipline.out.' % (NEXMDir)
            sys.exit()
        ## Delete previous permdipole.out if exists ##
        if os.path.exists('%s/%s/permdipole.out' % (cwd,NEXMDir)):
            os.remove('%s/%s/permdipole.out' % (cwd,NEXMDir))
        ## Grep dipoles from standard output ##
        subprocess.call(shlex.split('sh %s/getexcited_package/collectpermdipole.sh %d' % (pathtodipole, nstates + 2)))
        if not os.path.exists('%s/%s/permdipole.out' % (cwd,NEXMDir)):
            print 'Path %s/permdipole.out does not exist.' % (NEXMDir)
            sys.exit()
        print '%s dipoles in md.out extracted' % (NEXMDir)
        ## Another check to ensure dipole calculation ##
        with open('%s/%s/permdipole.out' % (cwd,NEXMDir),'r') as data:
            if len(data.readlines()) != tdipoles*(nstates + 3):
                print 'Path %s/permdipole.out is incomplete.' % (NEXMDir)
                sys.exit()
        ## Delete previous pop.out if exists ##
        if os.path.exists('%s/%s/pop.out' % (cwd,NEXMDir)):
            os.remove('%s/%s/pop.out' % (cwd,NEXMDir))
        ## Get states ##
        if boflag == 1: ## adiabatic
            states = np.int_(np.ones(edipoles)*statinit)
        if boflag == 0: ## nonadiabatic
            ## Check coefficient file exists ##
            if not os.path.exists('%s/%s/coeff-n.out' % (cwd,NEXMDir)):
                print 'Path %s/coeff-n.out does not exist.' % (NEXMDir)
                sys.exit()
            data = open('%s/%s/coeff-n.out' % (cwd,NEXMDir),'r')
            data = data.readlines()
            states = np.zeros(edipoles)
            index = 0
            for line in data[0:tscol:odata*cdata]:
                val = line.split()
                pes = np.int(val[0])
                time = np.around(np.float(val[1]), decimals = 3)
                ## Another check to ensure dipole calculation ##
                if time != times[index]:
                    print 'There is an inconsistency in time-step in %s/coeff-n.out at %.3f fs.' % (NEXMDir,times[index])
                    sys.exit()
                states[index] = pes
                index += 1
            ## If simulation becomes adiabatic after nonadiabatic ##
            if len(data[0:tscol:odata*cdata]) < edipoles:
                states[index::] = pes
        ## Save populations to file ##
        np.savetxt('pop.out', np.transpose([times,states]), fmt=['%10.5e','%d'])
    
    ## Collect user-defined vector in the molecule to determine dipole direction for ensembles ##
    if dynq == 0: ## ensemble
        os.chdir('%s' % (cwd))
        print 'Checking coordinate files for vector generation.  Please wait ...'
        ## Generation of error file ##
        error = open('%s/uservec_collection.err' % (cwd),'w')
        errflag = 0
        for NEXMD in NEXMDs:
            ## Check and open list of directories ##
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %s/dirlist1 does not exist.' % (NEXMD)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1,int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Check if trajectory directory exists ##
                if not os.path.exists('%s/%04d' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                ## Check if coordinate file exists ##
                if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d/coords.xyz does not exist.' % (NEXMD,dir)
                    errflag = 1
                    continue
                ## Find geometries ##
                data = open('%s/%04d/coords.xyz' % (NEXMD,dir),'r')
                data = data.readlines()
                lenc = len(data)
                ncoords = 0
                index = 0
                tflag1 = 0
                tflag2 = 0
                tflag3 = 0
                array = np.array([])
                for line in data:
                    if 'time' in line:
                        if ncoords == 0:
                            tinit = np.float(line.split()[-1])
                            if tinit != tinith:
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
                        array = np.append(array,index)
                    index += 1
                if tflag1 == 1:
                    print >> error, 'Initial time in %s%04d/coords.xyz does not match time_init in %s%04d/input.ceon.' % (NEXMD,dir)
                    errflag = 1
                    continue
                if tflag2 == 1:
                    print >> error, 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (NEXMD,dir)
                    errflag = 1
                    continue
                if tflag3 == 1:
                    array = np.append(array,index)
                else:
                    array = np.append(array, lenc + 1)
                array = np.int_(array)
                ## Checks to ensure user-defined vector calculation ##
                if ncoords == 0:
                    print >> error, 'No coordinates were found in %s%04d.' % (NEXMD,dir)
                    errflag = 1
                    continue
                if ncoords == 1:
                    print >> error, 'Only initial coordinates, at %.2f fs, were found in %s%04d.' % (tinit,NEXMD,dir)
                    errflag = 1
                    continue
                ## Generation of user-defined vector file ##
                output = open('%s/%04d/uservec.out' % (NEXMD,dir),'w')
                ## Print user-defined vector to file ##
                for ncoord in np.arange(ncoords):
                    coords = data[array[ncoord] + 1: array[ncoord + 1] - 1:1]
                    vec0 = np.float_(coords[lines[0]].split()[1:])
                    vec1 = np.float_(coords[lines[1]].split()[1:])
                    uvec = np.subtract(vec1, vec0)
                    output.write('{:>12}  {:>12}  {:>12}  {:>12}\n'.format(times[ncoord],uvec[0],uvec[1],uvec[2]))
                print '%s%04d user-defined vector from coords.xyz extracted' % (NEXMD,dir)
        if errflag == 1:
            print 'One or more trajectories have experienced an error, check uservec_collection.err.'
            contq = input('Continue? Answer yes [1] or no [0]: ')
            if contq not in [1,0]:
                print 'Answer must to be 1 or 0.'
                sys.exit()
            if contq == 0:
                sys.exit()
        else:
            os.remove('%s/uservec_collection.err' % (cwd))

    ## Collect user-defined vector in the molecule to determine dipole direction for single trajectory ##
    if dynq == 1: ## single trajectory
        os.chdir('%s' % (cwd))
        print 'Checking coordinate files for vector generation.  Please wait ...'
        ## Check if trajectory directory exists ##
        if not os.path.exists('%s' % (NEXMDir)):
            print 'Path %s does not exist.' % (NEXMDir)
            sys.exit()
        ## Check if coordinate file exists ##
        if not os.path.exists('%s/coords.xyz' % (NEXMDir)):
            print 'Path %s/coords.xyz does not exist.' % (NEXMDir)
            sys.exit()
        ## Find geometries ##
        data = open('%s/coords.xyz' % (NEXMDir),'r')
        data = data.readlines()
        lenc = len(data)
        ncoords = 0
        index = 0
        tflag1 = 0
        tflag2 = 0
        tflag3 = 0
        array = np.array([])
        for line in data:
            if 'time' in line:
                if ncoords == 0:
                    tinit = np.float(line.split()[-1])
                    if tinit != tinith:
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
                array = np.append(array,index)
            index += 1
        if tflag1 == 1:
            print >> error, 'Initial time in %s/coords.xyz does not match time_init in %s/input.ceon.' % (NEXMDir)
            sys.exit()
        if tflag2 == 1:
            print >> error, 'There is an inconsistency in time-step in %s/coords.xyz.' % (NEXMDir)
            sys.exit()
        if tflag3 == 1:
            array = np.append(array,index)
        else:
            array = np.append(array, lenc + 1)
        array = np.int_(array)
        ## Checks to ensure user-defined vector calculation ##
        if ncoords == 0:
            print >> error, 'No coordinates were found in %s/coords.xyz.' % (NEXMDir)
            sys.exit()
        if ncoords == 1:
            print >> error, 'Only initial coordinates, at %.2f fs, were found in %s/coords.xyz.' % (tinit,NEXMDir)
            sys.exit()
        ## Generation of user-defined vector file ##
        output = open('%s/uservec.out' % (NEXMDir),'w')
        ## Print user-defined vector to file ##
        for ncoord in np.arange(ncoords):
            coords = data[array[ncoord] + 1: array[ncoord + 1] - 1:1]
            vec0 = np.float_(coords[lines[0]].split()[1:])
            vec1 = np.float_(coords[lines[1]].split()[1:])
            uvec = np.subtract(vec1, vec0)
            output.write('{:>12}  {:>12}  {:>12}  {:>12}\n'.format(times[ncoord],uvec[0],uvec[1],uvec[2]))
        print '%s user-defined vector from coords.xyz extracted' % (NEXMDir)
            
    ## Collect dipole moment along a single trajectory ##
    if dynq == 1: ## single trajectory
        os.chdir('%s' % (cwd))
        print 'Collecting permanent dipole moment along a single trajectory.  Please wait ...'
        ## Generate output file ##
        output = open('%s/permdipole_single.out' % (cwd),'w')
        ttraj = 0
        ctraj = 0
        etraj = 0
        ## Determine completed number of time-steps ##
        if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
            print 'Path %s/energy-ev.out does not exist.' % (NEXMDir)
            sys.exit()
        data = open('%s/energy-ev.out' % (NEXMDir),'r')
        data = data.readlines()
        tsteps = len(data) - 1
        ## Generate array with indices of the dipole blocks along trajectory ##
        if tsteps >= tscol:
            if not os.path.exists('%s/permdipole.out' % (NEXMDir)):
                print 'Path %s/permdipole.out does not exist.' % (NEXMDir)
                sys.exit()
            data = open('%s/permdipole.out' % (NEXMDir),'r')
            data = data.readlines()
            lend = len(data)
            ndipoles = 0
            tflag = 0
            index = 0
            array = np.array([])
            for line in data:
                if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in line:
                    if ndipoles == 0:
                        time = tinith
                    else:
                        time += dt*odata*cdata
                        if time > tcoll:
                            tflag = 1
                            break
                    ndipoles += 1
                    array = np.append(array,index)
                index += 1
            ## Another check to ensure dipole calculation ##
            if ndipoles != edipoles:
                print 'Number of dipoles detected in %s/permdipole.out, %d, does not match the expected %d.' % (NEXMDir,ndipoles,edipoles)
                sys.exit()
            ## Append lines for last dipole set ##
            if tflag == 1:
                array = np.append(array, index)
            else:
                array = np.append(array, lend + 1)
            array = np.int_(array)
            ## Another check to ensure dipole calculation ##
            if ndipoles == 0:
                print 'No dipoles were found in %s/permdipole.out' % (NEXMDir)
                sys.exit()
            ## Open the user-defined vector file = [time, vx, vy, vz] ##
            if not os.path.exists('%s/uservec.out' % (NEXMDir)):
                print 'Path %s/uservec.out does not exist.' % (NEXMDir)
                sys.exit()
            uservec = np.genfromtxt('%s/uservec.out' % (NEXMDir))
            ## Open population data = [time, state] ##
            states = np.genfromtxt('%s/pop.out' % (NEXMDir))
            ## Collect dipole along a single trajectory ##
            sdipole = np.zeros((ndipoles,2))
            for ndipole in np.arange(ndipoles):
                dipoles = data[array[ndipole] + 1:array[ndipole + 1]:1]
                vdipole = np.float_(dipoles[np.int(states[ndipole, 1])].split()[2:5:1]) ## extracts dipole vector
                mdipole = np.float(dipoles[np.int(states[ndipole, 1])].split()[5]) ## extracts dipole magnitude
                cosine = np.arccos(np.dot(vdipole, uservec[ndipole,1:4:1])/(np.linalg.norm(vdipole)*np.linalg.norm(uservec[ndipole,1:4:1]))) if dotq == 1 else 0 ## inverse cosine of the dot product
                sdipole[ndipole] = np.array([mdipole, cosine])
            ## Delete extraneous data ##
            os.remove('%s/permdipole.out' % (NEXMDir))
            os.remove('%s/pop.out' % (NEXMDir))
            os.remove('%s/uservec.out' % (NEXMDir))
            print '%s' % (NEXMDir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
            ctraj += 1
            if tsteps == tsmax:
                etraj += 1
        else:
            print '%s' % (NEXMDir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
        ttraj += 1
        ## Summary of results ##
        if ctraj == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(tsmax)), tcoll)
        else:
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Completed trajectories:', '%04d' % (ctraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
            print >> output, 'Total trajectories: ', '%04d' % (ttraj)
            print >> output, 'Completed trajectories: ', '%04d' % (ctraj)
            print >> output, 'Excellent trajectories: ', '%04d' % (etraj)
            for ndipole in np.arange(ndipoles):
                print >> output, '%0*.2f' % (len(str((tsmax))) + 2, dt*odata*cdata*ndipole), '%d' % (states[ndipole, 1]), '%03.6f' % (sdipole[ndipole,0]), '%.6f' % (np.degrees(sdipole[ndipole,1]))

    ## Calculate mean dipole from ensemble of trajectories ##
    if dynq == 0 and typeq == 0: ## mean from ensemble
        os.chdir('%s' % (cwd))
        print 'Collecting mean permanent dipole moment from ensemble.  Please wait ...'
        ## Determine total number of trajectories in ensemble ##
        with open('%s/totdirlist' % (NEXMDir),'w') as data:
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                input = fileinput.input('%s/dirlist1' % (NEXMD))
                data.writelines(input)
        dirlist1 = np.int_(np.genfromtxt('%s/totdirlist' % (NEXMDir)))
        if isinstance(dirlist1,int) == true:
            dirlist1 = np.array([dirlist1])
        os.remove('%s/totdirlist' % (NEXMDir))
        ## Generate output and error files ##
        output = open('%s/permdipole_mean_ensemble.out' % (cwd),'w')
        error = open('%s/permdipole_mean_ensemble.err' % (cwd),'w')
        ## Generate dipole arrays for final results ##
        fdipole = np.zeros(len(times)) ## magnitude
        fcosine = np.zeros(len(times)) ## direction
        edipole = np.zeros((edipoles, len(dirlist1)))
        ecosine = np.zeros((edipoles, len(dirlist1)))
        ttraj = 0
        ctraj = 0
        etraj = 0
        errflag = 0
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1, int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine completed number of time-steps ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                data = data.readlines()
                tsteps = len(data) - 1
                ## Generate array with indices of the dipole blocks along a single trajectory ##
                if tsteps >= tscol:
                    if not os.path.exists('%s/%04d/permdipole.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/permdipole.out does not exist.' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    data = open('%s/%04d/permdipole.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    lend = len(data)
                    ndipoles = 0
                    tflag = 0
                    index = 0
                    array = np.array([])
                    for line in data:
                        if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in line:
                            if ndipoles == 0:
                                time = tinith
                            else:
                                time += dt*odata*cdata
                                if time > tcoll:
                                    tflag = 1
                                    break
                            ndipoles += 1
                            array = np.append(array,index)
                        index += 1
                    ## Check to ensure dipole calculation ##
                    if ndipoles != edipoles:
                        print >> error, 'Number of dipoles detected in %s%04d/permdipole.out, %d, does not match the expected %d.' % (NEXMD,dir,ndipoles,edipoles)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Append lines for last dipole set ##
                    if tflag == 1:
                        array = np.append(array,index)
                    else:
                        array = np.append(array, lend + 1)
                    array = np.int_(array)
                    ## Another check to ensure dipole calculation ##
                    if ndipoles == 0:
                        print >> error, 'No dipoles were found in %s%04d/permdipole.out' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Open the user-defined vector file = [time, vx, vy, vz] ##
                    if not os.path.exists('%s/%04d/uservec.out' % (NEXMD,dir)):
                        print >> error, 'Path %s/%04d/uservec.out' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    uservec = np.genfromtxt('%s/%04d/uservec.out' % (NEXMD,dir))
                    ## Open population data = [time, state] ##
                    states = np.genfromtxt('%s/%04d/pop.out' % (NEXMD,dir))
                    ## Collect dipole along a single trajectory ##
                    sdipole = np.zeros(ndipoles)
                    scosine = np.zeros(ndipoles)
                    for ndipole in np.arange(ndipoles):
                        dipoles = data[array[ndipole] + 1:array[ndipole + 1]:1]
                        vdipole = np.float_(dipoles[np.int(states[ndipole, 1])].split()[2:5:1]) ## extracts dipole vector
                        mdipole = np.float(dipoles[np.int(states[ndipole, 1])].split()[5]) ## extracts dipole magnitude
                        cosine = np.arccos(np.dot(vdipole, uservec[ndipole,1:4:1])/(np.linalg.norm(vdipole)*np.linalg.norm(uservec[ndipole,1:4:1]))) if dotq == 1 else 0 ## inverse cosine of the dot product
                        sdipole[ndipole] = mdipole
                        scosine[ndipole] = cosine
                        edipole[ndipole,ctraj] = sdipole[ndipole]
                        ecosine[ndipole,ctraj] = scosine[ndipole]
                    fdipole += sdipole
                    fcosine += scosine
                    ## Delete extraneous data ##
                    os.remove('%s/%04d/permdipole.out' % (NEXMD,dir))
                    os.remove('%s/%04d/pop.out' % (NEXMD,dir))
                    os.remove('%s/%04d/uservec.out' % (NEXMD,dir))
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ctraj += 1
                    if tsteps == tsmax:
                        etraj += 1
                else:
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                ttraj += 1
        ## Summary of results ##
        if ctraj == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(tsmax)), tcoll)
        else:
            ## Mean and standard deviation for dipole magnitude ##
            edipole = np.delete(edipole, np.arange(ctraj, ttraj), axis = 1)
            edipole = np.std(edipole, axis = 1)
            fdipole = fdipole/ctraj
            ## Mean and standard deviation for dipole direction ##
            ecosine = np.delete(ecosine, np.arange(ctraj, ttraj), axis = 1)
            ecosine = np.std(ecosine, axis = 1)
            fcosine = fcosine/ctraj
            ## Summary of results ##
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Completed trajectories:', '%04d' % (ctraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
            print >> output, 'Total trajectories:', '%04d' % (ttraj)
            print >> output, 'Completed trajectories:', '%04d' % (ctraj)
            print >> output, 'Excellent trajectories:', '%04d' % (etraj)
            for ndipole in np.arange(ndipoles):
                print >>  output, '%0*.2f' % (len(str((tsmax))) + 2, dt*odata*cdata*(ndipole)), '%03.6f' % (fdipole[ndipole]), '%03.6f' % (edipole[ndipole]), '%.6f' % (np.degrees(fcosine[ndipole])), '%.6f' % (np.degrees(ecosine[ndipole]))
        if errflag == 1:
            print 'One or more trajectories have experienced an error, check permdipole_mean_ensemble.err.'
        else:
            os.remove('%s/permdipole_mean_ensemble.err' % (cwd))
    
    ## Collect dipoles from ensemble of trajectories at all time-steps ##
    if dynq == 0 and typeq == 1: ## all from ensemble
        os.chdir('%s' % (cwd))
        print 'Collecting all permanent dipole moments from ensemble.  Please wait ...'
        ## Generate output and error files ##
        output = open('%s/permdipole_raw_ensemble.out' % (cwd),'w')
        error = open('%s/permdipole_raw_ensemble.err' % (cwd),'w')
        ttraj = 0
        etraj = 0
        errflag = 0
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(dirlist1, int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine number of time-steps completed ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                data = data.readlines()
                tsteps = len(data) - 1
                ## Generate array with indices of the dipole blocks along trajectory ##
                if not os.path.exists('%s/%04d/permdipole.out' % (NEXMD,dir)):
                    print >> error, 'Path %s%04d/permdipole.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/permdipole.out' % (NEXMD,dir),'r')
                data = data.readlines()
                lend = len(data)
                ndipoles = 0
                tflag = 0
                index = 0
                array = np.array([])
                for line in data:
                    if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in line:
                        if ndipoles == 0:
                            time = tinith
                        else:
                            time += dt*odata*cdata
                        ndipoles += 1
                        array = np.append(array,index)
                    index += 1
                ## Another check to ensure dipole calculation ##
                if ndipoles != edipoles:
                    print >> error, 'Number of dipoles detected in %s%04d/permdipole.out, %d, does not match the expected %d.' % (NEXMD,dir,ndipoles,edipoles)
                    errflag = 1
                    ttraj += 1
                    continue
                ## Append lines for last dipole set ##
                array = np.append(array, lend + 1)
                array = np.int_(array)
                ## Another check to ensure the dipole calculation ##
                if ndipoles == 0:
                    print >> error, 'No dipoles were found in %s%04d/permdipole.out.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                ## Open the user-defined vector file = [time, vx, vy, vz] ##
                if not os.path.exists('%s/%04d/uservec.out' % (NEXMD,dir)):
                    print >> error, 'Path %s/%04d/uservec.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                uservec = np.genfromtxt('%s/%04d/uservec.out' % (NEXMD,dir))
                ## Open population data = [time, state] ##
                if not os.path.exists('%s/%04d/pop.out' % (NEXMD,dir)):
                    print >> error, 'Path %s/%04d/pop.out does not exist.' % (NEXMD,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                states = np.genfromtxt('%s/%04d/pop.out' % (NEXMD,dir))
                ## Collect dipole along a single trajectory ##
                for ndipole in np.arange(ndipoles):
                    dipoles = data[array[ndipole] + 1:array[ndipole + 1]:1]
                    vdipole = np.float_(dipoles[np.int(states[ndipole, 1])].split()[2:5:1]) ## extracts dipole vector
                    mdipole = np.float(dipoles[np.int(states[ndipole, 1])].split()[5]) ## extracts dipole magnitude
                    cosine = np.arccos(np.dot(vdipole, uservec[ndipole,1:4:1])/(np.linalg.norm(vdipole)*np.linalg.norm(uservec[ndipole,1:4:1]))) if dotq == 1 else 0 ## inverse cosine of the dot product
                    print >> output,  '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, dt*odata*cdata*ndipole), '%d' % (states[ndipole, 1]), '%03.6f' % (mdipole), '%.6f' % (np.degrees(cosine))
                ## Delete extraneous data ##
                os.remove('%s/%04d/permdipole.out' % (NEXMD,dir))
                os.remove('%s/%04d/pop.out' % (NEXMD,dir))
                os.remove('%s/%04d/uservec.out' % (NEXMD,dir))
                print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                if tsteps == tsmax:
                    etraj += 1
                ttraj += 1
        ## Summary of results ##
        if ttraj == 0:
            print 'No trajectories completed with %0*.2f.' % (len(str(tsmax)), tcoll)
        else:
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
        if errflag == 1:
            print 'One or more trajectories have experienced an error, check permdipole_raw_ensemble.err.'
        else:
            os.remove('%s/permdipole_raw_ensemble.err' % (cwd))
