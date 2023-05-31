#/usr/bin/python

'''
    
This function collects pess and nacts from NEXMD during adiabatic or 
non-adiabatic dynamics.

A maximum total of four output files will be generated in the current 
working directory if collecting pess and nacts is requested.  In 
'pes_raw_ensemble.out' and 'nact_raw_ensemble.out', the pess and nacts
at all time-steps or up to a time defined by the user are shown, 
respectively.  In 'pes_raw_ensemble.out', columns from left to right are: 
directory of trajectory, current state, new state, followed by all pess.  
Likewise, 'nact_raw_ensemble.out' is in similar format.  The columns 
showing nacts are consecutive rows of the nact matrix, same as that 
shown in 'nact.out', located in the directory of each trajectory.  
In 'pes_hop_ensemble.out' and 'nact_hop_ensemble.out' are same data as 
'...raw_ensemble.out', but only at time-steps where hops occur.  If the 
BO flag (i.e. bo_dynamics_flag) in 'header' is set to '1', the 
simulation is adiabatic and only 'pes_raw_ensemble.out' will be 
generated.  An error file will be generated if certain files do not 
exist or if trajectories did not finish up to the user-defined length of 
analysis.

Type of calculation

[1] Ensemble of Trajectories

[1a] All time-steps
In 'pes_raw_ensemble.out':
> trajectory directory, current state, new state, pess
In 'nact_raw_ensemble.out':
> trajectory directory, current state, new state, nacts

[2a] All up to user-defined time
In 'pes_raw_ensemble.out':
> trajectory directory, current state, new state, pess
In 'nact_raw_ensemble.out':
> trajectory directory, current state, new state, nacts

[2b] Hops up to user-defined time
In 'pes_hop_ensemble.out':
> trajectory directory, current state, new state, pess
In 'nact_hop_ensemble.out':
> trajectory directory, current state, new state, nacts
Note: current state will not equal new state for ..._hops_...
output files

Output files:
- pes_[type].out, where [type] = hop_ensemble, raw_ensemble
- nact_[type].out, where [type] = hop_ensemble, raw_ensemble

Error files:
- pesnact_[type].out, where [type] = hop_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob

cwd = os.getcwd()

def pesnact(header):

    print('Collecting pess and/or nacts.')

    ## Type of calculaton and directory ##
    NEXMDir = input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print('Path %s does not exist.' % (NEXMDir))
        sys.exit()
    typeq = eval(input('Output pess and/or nacts at all time-steps and trajectories, up to some user-defined time, or hops only?\nanswer all [0], user-defined time [1], or hops [2]: '))
    if typeq not in [0,1,2]:
        print('Answer must be 0, 1, or 2.')
        sys.exit()

    ## Information from header ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print('Path %s/header does not exist.' % (NEXMDir))
        sys.exit()
    header = header('%s/header' % (NEXMDir))

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

    ## Check that dynamics are non-Born-Oppenheimer if data for hops are requested ##
    if header.bo_dynamics_flag == 1 and typeq == 2:
        print('Dynamics are set to Born-Oppenheimer (bo_dynamics_flag = 1), but hops occur during non-Born-Oppenheimer (bo_dynamics_flag = 0) dynamics.\nPlease check bo_dynamics_flag in header.')
        sys.exit()

    ## Redefine number of quantum steps per classical step if set to 0 in header ##
    if header.n_quant_steps == 0:
        header.n_quant_steps = 1

    ## Collection time ##
    if typeq == 0: ## all time-steps
        tcoll = (header.n_class_steps - 1)*header.time_step
    if typeq == 1 or typeq == 2: ## user-defined time or time-steps at hops only
        tcoll = eval(input('Collect data up to what time in femtoseconds: '))
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

    ## Number of classical steps ##
    tscol = 0
    while tscol*header.time_step*header.out_data_steps <= tcoll:
        tscol += 1

    ## Data type ##
    if header.bo_dynamics_flag == 0: ## non-adiabatic
        dtypeq = eval(input('Collect pess [1], nacts [2], or both [3]: '))
    if header.bo_dynamics_flag == 1: ## adiabatic
        dtypeq = 1
    if dtypeq not in [1,2,3]:
        print('Answer must be 1, 2, or 3.')
        sys.exit()

    ## Line number array ##
    if header.moldyn_verbosity == 3:
        ##if tstepq == 0 and header.n_quant_steps != 1:
        if  header.n_quant_steps != 1:
            linenums = np.arange(0, tscol*header.n_quant_steps - (header.n_quant_steps - 1), header.n_quant_steps)
        else:
            linenums = np.arange(0, tscol)
    else:
        linenums = np.arange(0, tscol)

    ## Time array ##
    if dtypeq == 1:
        times = np.around(np.linspace(header.time_init, tcoll, tscol), decimals = 3)
    if dtypeq == 2 or dtypeq == 3:
        times = np.around(np.linspace(header.time_init + header.time_step*header.out_data_steps, tcoll, tscol - 1), decimals = 3)

    ## Indices to cut extraneous nact data ##
    indices = np.array([])
    index = header.n_exc_states_propagate
    for term in np.split(np.arange(header.n_exc_states_propagate*header.n_exc_states_propagate),header.n_exc_states_propagate):
        indices = np.append(indices,term[-index::])
        index -= 1
    indices = np.int_(np.insert(indices + 1, 0, 0, 0))

    ## Check for NEXMD folders ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print('There are no NEXMD folders in %s.' % (NEXMDir))
        sys.exit()

    ### Adiabatic ###
    if header.bo_dynamics_flag == 1:
        ## All time-steps ##
        if typeq == 0:
            ## Generate output/error files ##
            pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
            error = open('%s/pesnact.err' % (cwd),'w')
            ## Begin looping over directories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('Path %dirlist1 does not exist.' % (NEXMD))
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
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## Determine if pes files exist and open them ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                        print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                    pes = pes.readlines()
                    lines = linenums[0:len(pes):1]
                    ## Collect data ##
                    tflag = 0
                    index = 0
                    for line in lines:
                        if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                            pess = np.float_(pes[int(line/header.n_quant_steps)].split())
                        else:
                            pess = np.float_(pes[line].split())

                        time = np.around(pess[0], decimals = 3)
                        if time != times[index]:
                            print('There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index]), file=error)
                            tflag = 1
                            errflag = 1
                            break
                        print('%s%04d' % (NEXMD,dir), '%d' % (header.exc_state_init), '%d' % (header.exc_state_init), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)
                        index += 1
                    if tflag == 0:
                        ctraj += 1
                        if tsteps == header.n_class_steps:
                            etraj += 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
                os.remove('%s/pes_raw_ensemble.out' % (cwd))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(header.n_class_steps)),tcoll))
            else:
                os.remove('%s/pesnact.err' % (cwd))

        ## User-defined time-steps ##
        if typeq == 1:
            ## Generate output/error files ##
            pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
            error = open('%s/pesnact.err' % (cwd),'w')
            ## Begin looping over directories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('Path %dirlist1 does not exist.' % (NEXMD))
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
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## Determine if pes files exist and open them ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                        print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                    pes = pes.readlines()
                    lines = linenums[0:len(pes):1]
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        tflag = 0
                        index = 0
                        for line in lines:
                            if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                                pess = np.float_(pes[int(line/header.n_quant_steps)].split())
                            else:
                                pess = np.float_(pes[line].split())

                            time = np.around(pess[0], decimals = 3)
                            if time != times[index]:
                                print('There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index]), file=error)
                                tflag = 1
                                errflag = 1
                                break
                            print('%s%04d' % (NEXMD,dir), '%d' % (header.exc_state_init), '%d' % (header.exc_state_init), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == header.n_class_steps:
                                etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
                os.remove('%s/pes_raw_ensemble.out' % (cwd))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(header.n_class_steps)),tcoll))
            else:
                os.remove('%s/pesnact.err' % (cwd))

    ### Non-adiabatic ###
    if header.bo_dynamics_flag == 0:
        ## All time-steps ##
        if typeq == 0:
            ## Generate output/error files ##
            if dtypeq == 1:
                pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
            if dtypeq == 2:
                nactall = open('%s/nact_raw_ensemble.out' % (cwd),'w')
            if dtypeq == 3:
                pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
                nactall = open('%s/nact_raw_ensemble.out' % (cwd),'w')
            error = open('%s/pesnact.err' % (cwd),'w')
            ## Begin looping over directories ##
            ttraj = 0
            ctraj = 0
            etraj = 0
            errflag = 0
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('path %sdirlist1 does not exist.' % (NEXMD))
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
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## Determine if pes/nact files exist and open them ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    hsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtypeq == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                    ## Collect data ##
                    tflag = 0
                    index = 0
                    cstate = np.int(hops[0].split()[0])
                    for line in lines:

                        if line <= hsteps-1:
                            nstate = np.int(hops[line].split()[0])

                        if dtypeq in [1,3]:
                            pess = np.float_(pes[line].split())
                            time = np.around(pess[0], decimals = 3)
                            if time != times[index]:
                                print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                tflag = 1
                                errflag = 1
                                break

                        if dtypeq in [2,3]:
                            headerSteps = (hsteps-1)
                            if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                                headerSteps = header.n_quant_steps*(hsteps-1)

                            if line <= headerSteps:
                                nacts = np.float_(nact[line - 1].split())[indices]
                                time = np.around(nacts[0], decimals = 3)
                                if time != times[index]:
                                    print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                    tflag = 1
                                    errflag = 1
                                    break

                        if dtypeq == 1:
                            print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)
                            cstate = nstate
                            index += 1
                            continue
                        if dtypeq == 2:
                            if line <= headerSteps:
                                print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nactall)
                                cstate = nstate
                                index += 1
                                continue
                        if dtypeq == 3:
                            print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)

                            if line <= headerSteps:
                                print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nactall)
                        cstate = nstate
                        index += 1
                    if tflag == 0:
                        ctraj += 1
                        if tsteps == header.n_class_steps:
                            etraj += 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
                if dtypeq == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One of more trajectories have experienced an error, check pesnact.err.') 
            else:
                os.remove('%s/pesnact.err' % (cwd))

        ## User-defined time-steps ##
        if typeq == 1:
            ## Generate output/error files ##
            if dtypeq == 1:
                pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
            if dtypeq == 2:
                nactall = open('%s/nact_raw_ensemble.out' % (cwd),'w')
            if dtypeq == 3:
                pesall = open('%s/pes_raw_ensemble.out' % (cwd),'w')
                nactall = open('%s/nact_raw_ensemble.out' % (cwd),'w')
            error = open('%s/pesnact.err' % (cwd),'w')
            ## Begin looping over directories ##
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
                    ## Determine number of completed time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir))
                        sys.exit()
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## Determine if pes/nact files exist ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    htsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtypeq == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        tflag = 0
                        index = 0
                        cstate = np.int(hops[0].split()[0])
                        for line in lines:

                            if line <= htsteps-1:
                                nstate = np.int(hops[line].split()[0])
                            if dtypeq in [1,3]:
                                pess = np.float_(pes[line].split())
                                time = np.around(pess[0], decimals = 3)
                                if time != times[index]:
                                    print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                    tflag = 1
                                    errflag = 1
                                    break
                            if dtypeq in [2,3]:
                                headerSteps = (htsteps-1)
                                if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                                    headerSteps = header.n_quant_steps*(htsteps-1)

                                if line <= headerSteps:
                                    nacts = np.float_(nact[(line-1)].split())[indices]
                                    time = np.around(nacts[0], decimals = 3)
                                    if time != times[index]:
                                        print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                        tflag = 1
                                        errflag = 1
                                        break
                            if dtypeq == 1:
                                print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)
                                cstate = nstate
                                index += 1
                                continue
                            if dtypeq == 2:
                                if line <= headerSteps:
                                    print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nactall)
                                    cstate = nstate
                                    index += 1
                                    continue
                            if dtypeq == 3:
                                print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=pesall)
                                if line <= headerSteps:
                                    print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nactall)
                            cstate = nstate
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == header.n_class_steps:
                                etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
                if dtypeq == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(header.n_class_steps)),tcoll))
            else:
                os.remove('%s/pesnact.err' % (cwd))

        ## Time-steps at hops only ##
        if typeq == 2:
            ## Generate output/error files ##
            if dtypeq == 1:
                peshop = open('%s/pes_hop_ensemble.out' % (cwd),'w')
            if dtypeq == 2:
                nacthop = open('%s/nact_hop_ensemble.out' % (cwd),'w')
            if dtypeq == 3:
                peshop = open('%s/pes_hop_ensemble.out' % (cwd),'w')
                nacthop = open('%s/nact_hop_ensemble.out' % (cwd),'w')
            error = open('%s/pesnact.err' % (cwd),'w')
            ## Begin looping over directories ##
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
                    ## determine number of completed time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print('Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir))
                        sys.exit()
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## determine if pes/nact files exist ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print('Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir), file=error)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    hsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtypeq == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print('Path %s%04d/pes.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print('Path %s%04d/nact.out does not exist.' % (NEXMD,dir), file=error)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                    ## Compare completed time-steps to collection time-steps and collect data ##
                    if tsteps >= tscol:
                        tflag = 0
                        cstate = np.int(hops[0].split()[0])
                        index = 0
                        for line in lines:
                            headerSteps = hsteps - 1
                            if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                                headerSteps = header.n_quant_steps*(hsteps-1)


                            if line <= headerSteps:
                                if(header.moldyn_verbosity ==3 and header.n_quant_steps >1):
                                    nstate = np.int(hops[int(line/header.n_quant_steps)].split()[0])
                                else:
                                    nstate = np.int(hops[int(line)].split()[0])
                            if dtypeq in [1,3]:
                                pess = np.float_(pes[line].split())
                                time = np.around(pess[0], decimals = 3)
                                if time != times[index]:
                                    print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                    tflag = 1
                                    errflag = 1
                                    break
                            if dtypeq in [2,3]:
                                if line <= headerSteps:
                                    nacts = np.float_(nact[line - 1].split())[indices]
                                    time = np.around(nacts[0], decimals = 3)
                                    if time != times[index]:
                                        print('There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index]), file=error)
                                        tflag = 1
                                        errflag = 1
                                        break
                            if nstate != cstate:
                                if dtypeq == 1:
                                    print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=peshop)
                                    cstate = nstate
                                    index += 1
                                    continue
                                if dtypeq == 2:
                                    if line <= headerSteps:
                                        print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nacthop)
                                        cstate = nstate
                                        index += 1
                                        continue
                                if dtypeq == 3:
                                    print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess), file=peshop)
                                    if line <= headerSteps:
                                        print('%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts), file=nacthop)
                            cstate = nstate
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == header.n_class_steps:
                                etraj += 1
                    else:
                        print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step), file=error)
                        errflag = 1
                    print('%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((header.n_class_steps))) + 2, (tsteps - 1)*header.time_step))
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print('No trajectories completed within %0*.2f fs.' % (len(str(header.n_class_steps)), tcoll))
                if dtypeq == 1:
                    os.remove('%s/pes_hop_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_hop_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_hop_ensemble.out' % (cwd))
                    os.remove('%s/nact_hop_ensemble.out' % (cwd))
            else:
                print('Total trajectories:', '%04d' % (ttraj))
                print('Completed trajectories:', '%04d' % (ctraj))
                print('Excellent trajectories:', '%04d' % (etraj))
            if errflag == 1:
                print('One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(header.n_class_steps)),tcoll))
            else:
                os.remove('%s/pesnact.err' % (cwd))
