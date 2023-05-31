#/usr/bin/python

'''
    
This function uses the diagonal elements of the transition
density in time to calculate fraction of induced charge and
width of exciton center of mass.

The user must supply how the molecule is to be divided.  This is 
done by defining atom numbers, according to how they are ordered 
in a general input file.  The atom numbers are to be written in 
a user-generated file, where the keyword 'break' is to be inserted 
wherever a new fragment of the molecule is defined.

Type of calculation:

[1] Single Trajectory
< sample input.ceon file with molecule, fragment file, collection time (fs)
> time (fs), occupied state (if non-adiabatic), occupancies in same order as fragment file (sums to 1), width of exciton center of mass

[2] Ensemble of Trajectories

[2a] Mean
< sample input.ceon file with molecule, fragment file, collection time (fs)
> time (fs), occupancies in same order as fragment file (sums to 1), width of exciton center of mass

[2b] All time-steps
< sample input.ceon file with molecule, fragment file, collection time (fs)
> trajectory directory, time (fs), occupied state (if non-adiabatic), occupancies in same order as fragment file (sums to 1), width of exciton center of mass

[2c] User-defined time
< sample input.ceon file with molecule, fragment file, collection time (fs)
> trajectory directory, time (fs), occupied state (if non-adiabatic), occupancies in same order as fragment file (sums to 1), width of exciton center of mass

Output Files:
- td_[type].out, where [type] = single, mean_ensemble, raw_ensemble

Error Files:
- td_[type].out, where [type] = single, mean_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob
import math

cwd = os.getcwd()

def tdiagonal(header):


    print('Calculating occupation of excitation as a function of time according to transition densities.')

    ## Type of calculation and directory check ##
    dynq = eval(input('Calculate occupation along one trajectory or an ensemble of trajectories?\nAnswer one [1] or ensemble [0]: '))
    if dynq not in [1,0]:
        print('Answer must be 1 or 0.')
        sys.exit()
    if dynq == 0: ## ensemble
        NEXMDir = input('Ensemble directory [e.g. NEXMD]: ')
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
        typeq = eval(input('Output mean occupation in time, at all time-steps and trajectories, or up to some user-defined time?\nAnswer mean [0], all [1], or user-defined [2]: '))
        if typeq not in [0,1,2]:
            print('Answer must be 0, 1, or 2.')
            sys.exit()
    if dynq == 1: ## single trajectory
        typeq = eval(input('Output occupation at all time-steps and trajectories, or up to some user-defined time?\nAnswer all [1], or user-defined [2]: '))
        if typeq not in [1,2]:
            print('Answer must be 1 or 2.')
            sys.exit()
        NEXMDir = input('Single trajectory directory: ')
        if not os.path.exists(NEXMDir):
            print('Path %s does not exist.' % (NEXMDir))
            sys.exit()
    
    ## Number and type of coordinates ##
    coordq = input('Directory with an input.ceon file for coordinates, do not include input.ceon in the path: ')
    if not os.path.exists('%s/input.ceon' % (coordq)):
        print('Path %s/input.ceon does not exist.' % (coordq))
        sys.exit()
    coords = open('%s/input.ceon' % (coordq),'r')
    coords = coords.readlines()
    index = 0
    for line in coords:
        if '&coord' in line:
            low = index
        if '&endcoord' in line:
            high = index
        index += 1
    natoms = (high - 1) - (low + 1) + 1
    orbitals = np.array([],dtype=np.int64)
    for line in coords[low + 1:high:1]:
        atype = np.int(line.split()[0])
        if atype == 1:
            orbitals = np.append(orbitals,1)
        else:
            orbitals = np.append(orbitals,4)
    ## Number of orbitals ##
    norbits = np.int(np.sum(orbitals))
    ## Indices to split the orbitals ##
    orbitals = np.cumsum(orbitals)
    ## Fragments of the molecule ##
    fragq = input('Directory with fragment file, include name of file in the path: ')
    if not os.path.exists(fragq):
        print('Path %s does not exist.' % (fragq))
        sys.exit()
    frag = open('%s' % (fragq),'r')
    frag = frag.readlines()
    barray = np.array([],dtype=np.int64)
    fragments = np.array([],dtype=np.int64)
    index = 0
    for line in frag:
        if 'break' in line:
            barray = np.append(barray, index - len(barray))
        else:
            fragments = np.append(fragments, np.int(line.split()[0])-1) ## - 1 for python indexing
        index += 1
    if len(fragments) != natoms:
        print('Number of atoms in the fragment file is not consistent with the coordinate file.')
    fragments = np.split(fragments, barray)
    nfrag = len(fragments)
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

    ## Check output data ##
    if header.out_data_steps == 0:
        print('No data have been printed to files because out_data_steps = 0 in %s/header.' % (NEXMDir))
        sys.exit()

    ## Check state is set for BO dynamics ##
    try:
        header.exc_state_init
        state_set = 1
    except AttributeError:
        state_set = 0
    if header.bo_dynamics_flag == 1 and state_set == 0:
        print('Dynamics are set to Born-Oppenheimer (bo_dynamics_flag = 1), but the initial state is not set.\nPlease check bo_dynamics_flag and exc_state_init in header.')
        sys.exit()

    ## Collection time ##
    if typeq == 0: ## mean occupation
        if dynq == 0: ## ensemble
            tcoll = eval(input('Calculate occupation up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: '))
        if dynq == 1: ## single trajectory
            tcoll = eval(input('Calculate occupation up to what time in femtoseconds? '))
        if isinstance(tcoll, int) == False and isinstance(tcoll, float) == False:
            print('Time must be integer or float.')
            sys.exit()
        if tcoll < 0:
            print('Time must be integer or float greater than zero.')
            sys.exit()
        tcoll = np.float(tcoll)
        nsteps = 1
        while nsteps*header.time_step <= tcoll:
            nsteps += 1
        tcoll = (nsteps - 1)*header.time_step
        if tcoll > (header.n_class_steps - 1)*header.time_step:
            tcoll = (header.n_class_steps - 1)*header.time_step
    if typeq == 1: ## all occupations
        tcoll = (header.n_class_steps - 1)*header.time_step
    if typeq == 2: ## all up to user-defined time
        tcoll = eval(input('Calculate occupation up to what time in femtoseconds? '))
        if isinstance(tcoll, int) == False and isinstance(tcoll, float) == False:
            print('Time must be integer or float.')
            sys.exit()
        if tcoll < 0:
            print('Time must be integer or float greater than zero.')
            sys.exit()
        tcoll = np.float(tcoll)
        nsteps = 1
        while nsteps*header.time_step <= tcoll:
            nsteps += 1
        tcoll = (nsteps - 1)*header.time_step
        if tcoll > (header.n_class_steps - 1)*header.time_step:
            tcoll = (header.n_class_steps - 1)*header.time_step

    x = (round(tcoll/(header.time_step*header.out_data_steps),0) - round(tcoll,3)/ round((header.time_step*header.out_data_steps),3))
    while(x  >=  header.time_step*0.1 or x <= - header.time_step*0.1):
        tcoll = tcoll - header.time_step
        x = (round(tcoll/(header.time_step*header.out_data_steps),0) - round(tcoll,3)/ round((header.time_step*header.out_data_steps),3))

    ## Number of classical steps ##
    tscol = 0
    while round(tscol*header.time_step*header.out_data_steps,3) <= round(tcoll,3):
        tscol += 1
    ## Line numbers and collection time array ##
    linenums = np.arange(tscol)
    times = np.around(np.linspace(header.time_init, tcoll, tscol), decimals  = 3)
    ## Single adiabatic trajectory ##
    if dynq == 1 and header.bo_dynamics_flag == 1:
        ## All time-steps - single adiabatic trajectory ##
        if typeq == 1:
            ## Generate output file ##
            if os.path.exists('%s/td_single.out' % (cwd)):
                os.remove('%s/td_single.out' % (cwd))
            ## Begin looping over trajectory ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            ## Determine completed number of time-steps ##
            if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
                print('Path %s/energy-ev.out does not exist.' % (NEXMDir))
                sys.exit()
            tsteps = np.genfromtxt('%s/energy-ev.out' % (NEXMDir), usecols=[0], skip_header=1).size
            ## Determine if transition density file exists ##
            if not os.path.exists('%s/transition-densities.out' % (NEXMDir)):
                print('Path %s/transition-densities.out does not exist.' % (NEXMDir))
                sys.exit()
            ## Check times ##
            times_td = np.around(np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=[0]), decimals = 3)
            if np.array_equal(times_td, times) == False:
                print('There is an inconsistency in time-step in %s/transition-densities.out.' % (NEXMDir))
                sys.exit()
            ## Collect data ##
            tds = np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=np.arange(1, norbits + 1))
            tds = np.hsplit(tds,orbitals)[0:-1:1]
            single_occ = np.zeros((tscol,natoms))
            aindex = 0
            for atom in tds:
                tindex = 0
                for time in atom:
                    single_occ[tindex, aindex] = np.sum(time)
                    tindex += 1
                aindex += 1
            single_occ = np.abs(single_occ)
            single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
            exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
            ## Generate fragment array ##
            frag_occ = np.zeros((tscol, nfrag))
            ## Split data by time ##
            index = 0
            for line in single_occ:
                ## Split data by fragments ##
                findex = 0
                for fragment in fragments:
                    frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                    findex += 1
                index += 1
            ## Output data ##
            dir_name = []
            dir_name.extend(['%s' % (NEXMDir)] * len(times_td))
            dtype = [('dir_name', list),('times_td', float)]
            for i in np.arange(nfrag + 2):
                dtype.append(('var%d' % (i), float))
            data = np.zeros(times_td.size, dtype=dtype)
            data['dir_name'] = dir_name
            data['times_td'] = times_td
            for i in np.arange(nfrag):
                data['var%d' % (i)] = frag_occ[:,i]
            data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
            data['var%d' % (nfrag + 1)] = exciton_com
            with open('%s/td_single.out' % (cwd),'a') as output:
                np.savetxt(output, data, fmt='%10s %07.2f '+ '%.3f ' * (nfrag + 2))
            ctraj += 1
            if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                etraj += 1
            print('%s' % (NEXMDir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
            ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
        
        ## User-defined time-steps - single adiabatic trajectory ##
        if typeq == 2:
            ## Generate output file ##
            if os.path.exists('%s/td_single.out' % (cwd)):
                os.remove('%s/td_single.out' % (cwd))
            ## Begin looping over trajectory ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            ## Determine completed number of time-steps ##
            if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
                print('Path %s/energy-ev.out does not exist.' % (NEXMDir))
                sys.exit()
            tsteps = np.genfromtxt('%s/energy-ev.out' % (NEXMDir), usecols=[0], skip_header=1).size
            ## Determine if transition density file exists ##
            if not os.path.exists('%s/transition-densities.out' % (NEXMDir)):
                print('Path %s/transition-densities.out does not exist.' % (NEXMDir))
                sys.exit()
            ## Check times ##
            times_td = np.around(np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=[0]), decimals = 3)
            if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                print('There is an inconsistency in time-step in %s/transition-densities.out.' % (NEXMDir))
                sys.exit()
            ## Compare completed time-steps to collection time-steps and collect data ##
            alength = round((tcoll/(header.time_step*header.out_data_steps)+1),0)
            if tsteps >= tscol:
                tds = np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=np.arange(1, norbits + 1))
                tds = np.hsplit(tds,orbitals)[0:-1:1]
                single_occ = np.zeros((tscol,natoms))
                aindex = 0
                for atom in tds:
                    tindex = 0
                    for time in atom:
                        if tindex > tscol-1:
                            break
                        single_occ[tindex, aindex] = np.sum(time)
                        tindex += 1
                    aindex += 1
                single_occ = np.abs(single_occ)
                single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                ## Generate fragment array ##
                frag_occ = np.zeros((tscol, nfrag))
                ## Split data by time ##
                index = 0
                for line in single_occ:
                    ## Split data by fragments ##
                    findex = 0
                    for fragment in fragments:
                        frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                        findex += 1
                    index += 1
                ## Output data ##
                dir_name = []
                dir_name.extend(['%s' % (NEXMDir)] * (int(alength)))
                dtype = [('dir_name', list),('times_td', float)]
                for i in np.arange(nfrag + 2):
                    dtype.append(('var%d' % (i), float))
                data = np.zeros(alength, dtype=dtype)
                data['dir_name'] = dir_name
                timearray = np.linspace(0,round(tcoll,2),alength)
                data['times_td'] = timearray
                vari1 = data['var%d' % (nfrag)]
                vari2= np.sum(frag_occ, axis = 1)
                for i in np.arange(nfrag):
                    data['var%d' % (i)] = frag_occ[:,i]
                data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                data['var%d' % (nfrag + 1)] = exciton_com
                with open('%s/td_single.out' % (cwd),'a') as output:
                    np.savetxt(output, data, fmt='%10s %07.2f '+ '%.3f ' * (nfrag + 2))
                ctraj += 1
                if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                    etraj += 1
            print('%s' % (NEXMDir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
            ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))

    ## Single non-adiabatic trajectory ##
    if dynq == 1 and header.bo_dynamics_flag == 0:
        ## All time-steps - single non-adiabatic trajectory ##
        if typeq == 1:
            ## Generate output file ##
            if os.path.exists('%s/td_single.out' % (cwd)):
                os.remove('%s/td_single.out' % (cwd))
            ## Begin looping over trajectory ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            ## Determine completed number of time-steps ##
            if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
                print('Path %s/energy-ev.out does not exist.' % (NEXMDir))
                sys.exit()
            tsteps = np.genfromtxt('%s/energy-ev.out' % (NEXMDir), usecols=[0], skip_header=1).size
            ## Determine if coefficient files exist and open it ##
            if not os.path.exists('%s/coeff-n.out' % (NEXMDir)):
                print('Path %s/coeff-n.out does not exist.' % (NEXMDir))
                sys.exit()
            hops = np.int_(np.genfromtxt('%s/coeff-n.out' % (NEXMDir),usecols=[0]))
            if np.where(linenums > hops.size)[0].size:
                hops = np.append(hops, hops[-1]*np.ones(np.where(linenums > hops.size)[0].size))
            else:
                hops = hops[linenums]
            ## Determine if transition density file exists ##
            if not os.path.exists('%s/transition-densities.out' % (NEXMDir)):
                print('Path %s/transition-densities.out does not exist.' % (NEXMDir))
                sys.exit()
            ## Check times ##
            times_td = np.around(np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=[0]), decimals = 3)
            if np.array_equal(times_td, times) == False:
                print('There is an inconsistency in time-step in %s/transition-densities.out.' % (NEXMDir))
                sys.exit()
            ## Collect data ##
            tds = np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=np.arange(1, norbits + 1))
            tds = np.hsplit(tds,orbitals)[0:-1:1]
            single_occ = np.zeros((tscol,natoms))
            aindex = 0
            for atom in tds:
                tindex = 0
                for time in atom:
                    single_occ[tindex, aindex] = np.sum(time)
                    tindex += 1
                aindex += 1
            single_occ = np.abs(single_occ)
            single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
            exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
            ## Generate fragment array ##
            frag_occ = np.zeros((tscol, nfrag))
            ## Split data by time ##
            index = 0
            for line in single_occ:
                ## Split data by fragments ##
                findex = 0
                for fragment in fragments:
                    frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                    findex += 1
                index += 1
            ## Output data ##
            dir_name = []
            dir_name.extend(['%s' % (NEXMDir)] * len(times_td))
            dtype = [('dir_name', list),('times_td', float), ('hops', float)]
            for i in np.arange(nfrag + 2):
                dtype.append(('var%d' % (i), float))
            data = np.zeros(times_td.size, dtype=dtype)
            data['dir_name'] = dir_name
            data['times_td'] = times_td
            data['hops'] = hops
            for i in np.arange(nfrag):
                data['var%d' % (i)] = frag_occ[:,i]
            data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
            data['var%d' % (nfrag + 1)] = exciton_com
            with open('%s/td_single.out' % (cwd),'a') as output:
                np.savetxt(output, data, fmt='%10s %07.2f %d '+ '%.3f ' * (nfrag + 2))
            ctraj += 1
            if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                etraj += 1
            print('%s' % (NEXMDir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
            ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))

        ## User-defined time-steps - single non-adiabatic trajectory ##
        if typeq == 2:
            ## Generate output files ##
            if os.path.exists('%s/td_single.out' % (cwd)):
                os.remove('%s/td_single.out' % (cwd))
            ## Begin looping over trajectory ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            ## Determine completed number of time-steps ##
            if not os.path.exists('%s/energy-ev.out' % (NEXMDir)):
                print('Path %s/energy-ev.out does not exist.' % (NEXMDir))
                sys.exit()
            tsteps = np.genfromtxt('%s/energy-ev.out' % (NEXMDir), usecols=[0], skip_header=1).size
            ## Determine if coefficient files exist and open it ##
            if not os.path.exists('%s/coeff-n.out' % (NEXMDir)):
                print('Path %s/coeff-n.out does not exist.' % (NEXMDir))
                sys.exit()
            hops = np.int_(np.genfromtxt('%s/coeff-n.out' % (NEXMDir),usecols=[0]))
            if np.where(linenums > hops.size)[0].size:
                hops = np.append(hops, hops[-1]*np.ones(np.where(linenums > hops.size)[0].size))
            else:
                hops = hops[linenums]
            ## Determine if transition density file exists ##
            if not os.path.exists('%s/transition-densities.out' % (NEXMDir)):
                print('Path %s/transition-densities.out does not exist.' % (NEXMDir))
                sys.exit()
            ## Check times ##
            times_td = np.around(np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=[0]), decimals = 3)
            if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                print('There is an inconsistency in time-step in %s/transition-densities.out.' % (NEXMDir))
                sys.exit()
            ## Compare completed time-steps to collection time-steps and collect data ##
            alength = int(round((tcoll/(header.time_step*header.out_data_steps)+1),0))
            if tsteps >= tscol:
                tds = np.genfromtxt('%s/transition-densities.out' % (NEXMDir),usecols=np.arange(1, norbits + 1))
                tds = np.hsplit(tds,orbitals)[0:-1:1]
                single_occ = np.zeros((tscol,natoms))
                aindex = 0
                for atom in tds:
                    tindex = 0
                    for time in atom:
                        if tindex > tscol-1:
                            break
                        single_occ[tindex, aindex] = np.sum(time)
                        tindex += 1
                    aindex += 1
                single_occ = np.abs(single_occ)
                single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                ## Generate fragment array ##
                frag_occ = np.zeros((tscol, nfrag))
                ## Split data by time ##
                index = 0
                for line in single_occ:
                    ## split data by fragments ##
                    findex = 0
                    for fragment in fragments:
                        frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                        findex += 1
                    index += 1
                ## Output data ##
                dir_name = []
                dir_name.extend(['%s' % (NEXMDir)] * (int(alength)))
                dtype = [('dir_name', list),('times_td', float), ('hops', float)]
                for i in np.arange(nfrag + 2):
                    dtype.append(('var%d' % (i), float))
                data = np.zeros(alength, dtype=dtype)
                timearray = np.linspace(0,round(tcoll,2),alength)
                data['dir_name'] = dir_name
                data['times_td'] = timearray
                data['hops'] = hops
                vari1 = data['var%d' % (nfrag)]
                vari2= np.sum(frag_occ, axis = 1)
                for i in np.arange(nfrag):
                    data['var%d' % (i)] = frag_occ[:,i]
                data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                data['var%d' % (nfrag + 1)] = exciton_com
                with open('%s/td_single.out' % (cwd),'a') as output:
                    np.savetxt(output, data, fmt='%10s %07.2f %d '+ '%.3f ' * (nfrag + 2))
                ctraj += 1
                if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                    etraj += 1
            print('%s' % (NEXMDir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
            ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))

    ## Adiabatic ensemble ##
    if dynq == 0 and header.bo_dynamics_flag == 1:
        ## All time-steps - adiabatic ensemble ##
        if typeq == 1:
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))
            error = open('%s/td_raw_ensemble.err' % (cwd),'w')
            ## Begin looping over trajectories ##
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
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if transition density file exists ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Check times ##
                    times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                    if np.array_equal(times_td, times) == False:
                        print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Collect data ##
                    tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1, norbits + 1))
                    tds = np.hsplit(tds,orbitals)[0:-1:1]
                    single_occ = np.zeros((tscol,natoms))
                    aindex = 0
                    for atom in tds:
                        tindex = 0
                        for time in atom:
                            single_occ[tindex, aindex] = np.sum(time)
                            tindex += 1
                        aindex += 1
                    single_occ = np.abs(single_occ)
                    single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                    exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                    ## Generate fragment array ##
                    frag_occ = np.zeros((tscol, nfrag))
                    ## Split data by time ##
                    index = 0
                    for line in single_occ:
                        ## Split data by fragments ##
                        findex = 0
                        for fragment in fragments:
                            frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                            findex += 1
                        index += 1
                    ## Output data ##
                    dir_name = []
                    dir_name.extend(['%s%04d' % (NEXMD,dir)] * len(times_td))
                    dtype = [('dir_name', list),('times_td', float)]
                    for i in np.arange(nfrag + 2):
                        dtype.append(('var%d' % (i), float))
                    data = np.zeros(times_td.size, dtype=dtype)
                    data['dir_name'] = dir_name
                    data['times_td'] = times_td
                    for i in np.arange(nfrag):
                        data['var%d' % (i)] = frag_occ[:,i]
                    data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                    data['var%d' % (nfrag + 1)] = exciton_com
                    with open('%s/td_raw_ensemble.out' % (cwd),'a') as output:
                        np.savetxt(output, data, fmt='%10s %07.2f '+ '%.3f ' * (nfrag + 2))
                    ctraj += 1
                    if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                        etraj += 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories have experienced an error, check td_raw_ensemble.err.')
            else:
                os.remove('%s/td_raw_ensemble.err' % (cwd))

        ## User-defined time-steps - adiabatic ensemble ##
        if typeq == 2:
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))
            error = open('%s/td_raw_ensemble.err' % (cwd),'w')
            ## Begin looping over trajectories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('path %sdirlist1 does not exist.' % (NEXMD))
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1, int) == True:
                    dirlist1 = np.array([dirlist1])
                for dir in dirlist1:
                    alength = round((tcoll/(header.time_step*header.out_data_steps)+1),0)
                    timearray = np.linspace(0,round(tcoll,2),alength)
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if transition density file exists ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Check times ##
                    times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                    if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                        print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1, norbits + 1))
                        tds = np.hsplit(tds,orbitals)[0:-1:1]
                        single_occ = np.zeros((tscol,natoms))
                        aindex = 0
                        for atom in tds:
                            tindex = 0
                            for time in atom:
                                if tindex > tscol-1:
                                    break
                                single_occ[tindex, aindex] = np.sum(time)
                                tindex += 1
                            aindex += 1
                        single_occ = np.abs(single_occ)
                        single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                        exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                        ## Generate fragment array ##
                        frag_occ = np.zeros((tscol, nfrag))
                        ## Split data by time ##
                        index = 0
                        for line in single_occ:
                            ## Split data by fragments ##
                            findex = 0
                            for fragment in fragments:
                                frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                                findex += 1
                            index += 1
                        ## Output data ##
                        dir_name = []
                        dir_name.extend(['%s%04d' % (NEXMD,dir)] * (int(alength)))
                        dtype = [('dir_name', list),('times_td', float)]
                        for i in np.arange(nfrag + 2):
                            dtype.append(('var%d' % (i), float))
                        data = np.zeros(alength, dtype=dtype)
                        data['dir_name'] = dir_name
                        data['times_td'] = timearray
                        len(data['var%d' % (i)])
                        vari1 = data['var%d' % (nfrag)]
                        vari2= np.sum(frag_occ, axis = 1)
                        for i in np.arange(nfrag):
                            data['var%d' % (i)] = frag_occ[:,i]
                        data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                        data['var%d' % (nfrag + 1)] = exciton_com
                        with open('%s/td_raw_ensemble.out' % (cwd),'a') as output:
                            np.savetxt(output, data, fmt='%10s %07.2f '+ '%.3f ' * (nfrag + 2))
                        ctraj += 1
                        if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                            etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step)*header.out_data_steps, file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories have experienced an error, check td_raw_ensemble.err.')
            else:
                os.remove('%s/td_raw_ensemble.err' % (cwd))

        ## Calculate mean - adiabatic ensemble ##
        if typeq == 0:
            ## Generate occupancy array for final results ##
            final_occ = np.zeros((tscol,natoms))
            final_exciton_com = np.zeros(tscol)
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))
            error = open('%s/td_mean_ensemble.err' % (cwd),'w')
            output = open('%s/td_mean_ensemble.out' % (cwd),'w')
            ## Begin looping over trajectories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('Path %sdirlist1 does not exist.' % (NEXMD))
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1, int) == True:
                    dirlist1 = np.array([dirlist1])
                for dir in dirlist1:
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if transition density files exist ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                        if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                            print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1,norbits + 1))
                        tds = np.hsplit(tds,orbitals)[0:-1:1]
                        single_occ = np.zeros((tscol,natoms))
                        aindex = 0
                        for atom in tds:
                            tindex = 0
                            for time in atom:
                                if tindex > tscol-1:
                                    break
                                single_occ[tindex, aindex] = np.sum(time)
                                tindex += 1
                            aindex += 1
                        single_occ = np.abs(single_occ)
                        single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                        exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                        final_occ += single_occ
                        final_exciton_com += exciton_com
                        ctraj += 1
                        if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                            etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                ## Average occupation ##
                final_occ = final_occ/ctraj
                final_exciton_com = final_exciton_com/ctraj
                ## Generate fragment array ##
                frag_occ = np.zeros((tscol, nfrag))
                ## Split data by atoms ##
                index = 0
                for line in final_occ:
                    ## Split data by fragments ##
                    findex = 0
                    for fragment in fragments:
                        frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                        findex += 1
                    index += 1
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
                print('Total trajectories:', '%04d' % (ttraj), file=output)
                print('Completed trajectories:', '%04d' % (ctraj), file=output)
                print('Excellent trajectories:', '%04d' % (etraj), file=output)
                for tstep in np.arange(tscol):
                    print('%0*.2f' % (len(str((header.n_class_steps))) + 2, tstep*header.time_step*header.out_data_steps), ' '.join(str('%.3f' % (x)) for x in frag_occ[tstep]), '%.3f' % (np.sum(frag_occ[tstep])), '%.3f' % (final_exciton_com[tstep]), file=output)
            if errflag == 1:
                print('One or more trajectories have experienced an error, check td_mean_ensemble.err')
            else:
                os.remove('%s/td_mean_ensemble.err' % (cwd))
                
    ## Non-adiabtic ensemble ##
    if dynq == 0 and header.bo_dynamics_flag == 0:
        ## All time-steps - non-adiabatic ensemble ##
        if typeq == 1:
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))

            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            error = open('%s/td_raw_ensemble.err' % (cwd),'w')
            ## Begin looping over trajectories ##
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
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if coefficient files exist and open them ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = np.int_(np.genfromtxt('%s/%04d/coeff-n.out' % (NEXMD,dir),usecols=[0]))
                    if np.where(linenums > hops.size)[0].size:
                        hops = np.append(hops, hops[-1]*np.ones(np.where(linenums > hops.size)[0].size))
                    else:
                        hops = hops[linenums]
                    ## Determine if transition density file exists ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Check times ##
                    times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                    if np.array_equal(times_td, times) == False:
                        print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Collect data ##
                    tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1, norbits + 1))
                    tds = np.hsplit(tds,orbitals)[0:-1:1]
                    single_occ = np.zeros((tscol,natoms))
                    aindex = 0
                    for atom in tds:
                        tindex = 0
                        for time in atom: 
                            single_occ[tindex, aindex] = np.sum(time)
                            tindex += 1
                        aindex += 1
                    single_occ = np.abs(single_occ)
                    single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                    exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                    ## Generate fragment array ##
                    frag_occ = np.zeros((tscol, nfrag))
                    ## Split data by time ##
                    index = 0
                    for line in single_occ:
                        ## Split data by fragments ##
                        findex = 0
                        for fragment in fragments:
                            frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                            findex += 1
                        index += 1
                    ## Output data ##
                    dir_name = []
                    dir_name.extend(['%s%04d' % (NEXMD,dir)] * len(times_td))
                    dtype = [('dir_name', list ),('times_td', float), ('hops', float)]
                    for i in np.arange(nfrag + 2):
                        dtype.append(('var%d' % (i), float))
                    data = np.zeros(times_td.size, dtype=dtype)
                    data['dir_name'] = dir_name
                    data['times_td'] = times_td
                    data['hops'] = hops
                    for i in np.arange(nfrag):
                        data['var%d' % (i)] = frag_occ[:,i]
                    data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                    data['var%d' % (nfrag + 1)] = exciton_com
                    with open('%s/td_raw_ensemble.out' % (cwd),'a') as output:
                        np.savetxt(output, data, fmt='%10s %07.2f %d '+ '%.3f ' * (nfrag + 2))
                    ctraj += 1
                    if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                        etraj += 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories have experienced an error, check td_raw_ensemble.err.')
            else:
                os.remove('%s/td_raw_ensemble.err' % (cwd))

        ## User-defined time-steps - non-adiabatic ensemble ##
        if typeq == 2:
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))

            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            error = open('%s/td_raw_ensemble.err' % (cwd),'w')
            ## Begin looping over trajectories ##
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
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if coefficient files exist and open them ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = np.int_(np.genfromtxt('%s/%04d/coeff-n.out' % (NEXMD,dir),usecols=[0]))
                    if np.where(linenums > hops.size)[0].size:
                        hops = np.append(hops, hops[-1]*np.ones(np.where(linenums > hops.size)[0].size))
                    else:
                        hops = hops[linenums]
                    ## Determine if transition density file exists ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Check times ##
                    times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                    if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                        print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1, norbits + 1))
                        tds = np.hsplit(tds,orbitals)[0:-1:1]
                        single_occ = np.zeros((tscol,natoms))
                        aindex = 0
                        for atom in tds:
                            tindex = 0
                            for time in atom:
                                if tindex > tscol-1:
                                    break
                                single_occ[tindex, aindex] = np.sum(time)
                                tindex += 1
                            aindex += 1
                        single_occ = np.abs(single_occ)
                        single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                        exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                        ## Generate fragment array ##
                        frag_occ = np.zeros((tscol, nfrag))
                        ## Split data by time ##
                        index = 0
                        for line in single_occ:
                            ## Split data by fragments ##
                            findex = 0
                            for fragment in fragments:
                                frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                                findex += 1
                            index += 1
                        ## Output data ##
                        alength = int(round((tcoll/(header.time_step*header.out_data_steps)+1),0))
                        dir_name = []
                        dir_name.extend(['%s%04d' % (NEXMD,dir)] * int(alength))
                        dtype = [('dir_name', list),('times_td', float), ('hops', float)]
                        for i in np.arange(nfrag + 2):
                            dtype.append(('var%d' % (i), float))
                        data = np.zeros(alength, dtype=dtype)
                        timearray = np.linspace(0,round(tcoll,2),alength)
                        data['dir_name'] = dir_name
                        data['times_td'] = timearray
                        data['hops'] = hops
                        vari1 = data['var%d' % (nfrag)]
                        vari2= np.sum(frag_occ, axis = 1)
                        len(data['var%d' % (i)])
                        for i in np.arange(nfrag):
                            data['var%d' % (i)] = frag_occ[:,i]
                        data['var%d' % (nfrag)] = np.sum(frag_occ, axis = 1)
                        data['var%d' % (nfrag + 1)] = exciton_com
                        with open('%s/td_raw_ensemble.out' % (cwd),'a') as output:
                            np.savetxt(output, data, fmt='%10s %07.2f %d '+ '%.3f ' * (nfrag + 2))
                        ctraj += 1
                        if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                            etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One of more trajectories have experienced an error, check td_raw_ensemble.err.')
            else:
                os.remove('%s/td_raw_ensemble.err' % (cwd))
        
        ## Calculate mean - non-adiabatc ensemble ##
        if typeq == 0:
            ## Generate occupancy array for final results ##
            final_occ = np.zeros((tscol,natoms))
            final_exciton_com = np.zeros(tscol)
            ## Generate output/error files ##
            if os.path.exists('%s/td_raw_ensemble.out' % (cwd)):
                os.remove('%s/td_raw_ensemble.out' % (cwd))
            if os.path.exists('%s/td_raw_ensemble.err' % (cwd)):
                os.remove('%s/td_raw_ensemble.err' % (cwd))

            output = open('%s/td_mean_ensemble.out' % (cwd),'w')
            error = open('%s/td_mean_ensemble.err' % (cwd),'w')
            ## Begin looping over trajectories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('path %sdirlist1 does not exist.' % (NEXMD))
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1, int) == True:
                    dirlist1 = np.array([dirlist1])
                for dir in dirlist1:
                    ## Determine completed number of time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    tsteps = np.genfromtxt('%s/%04d/energy-ev.out' % (NEXMD,dir), usecols=[0], skip_header=1).size
                    ## Determine if transition density files exist ##
                    if not os.path.exists('%s/%04d/transition-densities.out' % (NEXMD,dir)):
                        print('Path %s%04d/transition-densities.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        times_td = np.around(np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=[0]), decimals = 3)
                        if len(times_td) != math.floor((header.n_class_steps)/(header.out_data_steps)):
                            print('There is an inconsistency in time-step in %s%04d/transition-densities.out.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        tds = np.genfromtxt('%s/%04d/transition-densities.out' % (NEXMD,dir),usecols=np.arange(1,norbits + 1))
                        tds = np.hsplit(tds,orbitals)[0:-1:1]
                        single_occ = np.zeros((tscol,natoms))
                        aindex = 0
                        for atom in tds:
                            tindex = 0
                            for time in atom:
                                if tindex > tscol-1:
                                    break
                                single_occ[tindex, aindex] = np.sum(time)
                                tindex += 1
                            aindex += 1
                        single_occ = np.abs(single_occ)
                        single_occ = single_occ/single_occ.sum(axis = 1, keepdims = True)
                        exciton_com = 1.0/np.sum(np.square(single_occ), axis = 1)
                        final_occ += single_occ
                        final_exciton_com += exciton_com
                        ctraj += 1
                        if tsteps == math.floor((header.n_class_steps)/(header.out_data_steps)):
                            etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str(header.n_class_steps)) + 2, (tsteps - 1)*header.time_step*header.out_data_steps))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
            else:
                ## Average occupation ##
                final_occ = final_occ/ctraj
                final_exciton_com = final_exciton_com/ctraj
                ## Generate fragment array ##
                frag_occ = np.zeros((tscol, nfrag))
                ## Split data by atoms ##
                index = 0
                for line in final_occ:
                    ## Split data by fragments ##
                    findex = 0
                    for fragment in fragments:
                        frag_occ[index, findex] = np.sum(line[np.int_(fragment)])
                        findex += 1
                    index += 1    
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
                print('Total trajectories:', '%04d' % (ttraj), file=output)
                print('Completed trajectories:', '%04d' % (ctraj), file=output)
                print('Excellent trajectories:', '%04d' % (etraj), file=output)
                for tstep in np.arange(tscol):
                    print('%0*.2f' % (len(str((header.n_class_steps))) + 2, tstep*header.time_step*header.out_data_steps), ' '.join(str('%.3f' % (x)) for x in frag_occ[tstep]), '%.3f' % (np.sum(frag_occ[tstep])), '%.3f' % (final_exciton_com[tstep]), file=output)
            if errflag == 1:
                print('One or more trajectories have experienced an error, check td_mean_ensemble.err')
            else:
                os.remove('%s/td_mean_ensemble.err' % (cwd))
                
