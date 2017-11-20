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

def pesnact():

    print 'Collecting pess and/or nacts.'

    ## Type of calculaton and directory ##
    NEXMDir = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDir):
        print 'Path %s does not exist.' % (NEXMDir)
        sys.exit()
    typeq = input('Output pess and/or nacts at all time-steps and trajectories, up to some user-defined time, or hops only?\nanswer all [0], user-defined time [1], or hops [2]: ')
    if typeq not in [0,1,2]:
        print 'Answer must be 0, 1, or 2.'
        sys.exit()

    ## Information from header ##
    if not os.path.exists('%s/header' % (NEXMDir)):
        print 'Path %s/header does not exist.' % (NEXMDir)
        sys.exit()
    header = open('%s/header'% (NEXMDir),'r')
    header = header.readlines()
    num = 0
    tline = len(header)
    verb = none
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
    if boflag == 1 and stateinit == none:
        print 'Dynamics are set to Born-Oppenheimer, but the initial state is not set.\nPlease check bo_dynamics_flag and exc_state_init in header.'
        sys.exit()
    if boflag == 1 and typeq == 2:
        print 'Dynamics are set to Born-Oppenheimer. hops only occur during non-born-oppenheimer dynamics.\nPlease check bo_dynamics_flag in header.'
        sys.exit()

    ## Collection time ##
    if typeq == 0: ## all time-steps
        tcoll = (tsmax - 1)*dt
    if typeq == 1 or typeq == 2: ## user-defined time or time-steps at hops only
        tcoll = input('Collect data up to what time in femtoseconds: ')
        if isinstance(tcoll, int) == false and isinstance(tcoll, float) == false:
            print 'Time must be integer or float.'
            sys.exit()
        if tcoll < 0:
            print 'Time must be integer or float greater than zero.'
            sys.exit()
        tcoll = np.float(tcoll)
        nsteps = 1
        while nsteps*dt <= tcoll:
            nsteps += 1
        tcoll = (nsteps - 1)*dt
        if tcoll > (tsmax - 1)*dt:
            tcoll = (tsmax - 1)*dt

    ## Number of classical steps ##
    tscol = 0
    while tscol*dt*odata <= tcoll:
        tscol += 1

    ## Data type ##
    if boflag == 0: ## non-adiabatic
        dtypeq = input('Collect pess [1], nacts [2], or both [3]: ')
    if boflag == 1: ## adiabatic
        dtypeq = 1
    if dtypeq not in [1,2,3]:
        print 'Answer must be 1, 2, or 3.'
        sys.exit()

    ## Line number array ##
    if verb == 3:
        if tstepq == 0 and nqstep != 1:
            linenums = np.arange(nqstep, tscol*nqstep - (nqstep - 1), nqstep)
        else:
            linenums = np.arange(0, tscol)
    else:
        linenums = np.arange(0, tscol)

    ## Time array ##
    if dtypeq == 1:
        times = np.around(np.linspace(tinith, tcoll, tscol), decimals = 3)
    if dtypeq == 2 or dtypeq == 3:
        times = np.around(np.linspace(tinith + dt*odata, tcoll, tscol - 1), decimals = 3)

    ## Indices to cut extraneous nact data ##
    indices = np.array([])
    index = nstates
    for term in np.split(np.arange(nstates*nstates),nstates):
        indices = np.append(indices,term[-index::])
        index -= 1
    indices = np.int_(np.insert(indices + 1, 0, 0, 0))

    ## Check for NEXMD folders ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (NEXMDir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDir)
        sys.exit()

    ### Adiabatic ###
    if boflag == 1:
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
                    print 'Path %dirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1,int) == true:
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
                    ## Determine if pes files exist and open them ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
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
                        pess = np.float_(pes[line].split())
                        time = np.around(pess[0], decimals = 3)
                        if time != times[index]:
                            print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index])
                            tflag = 1
                            errflag = 1
                            break
                        print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (stateinit), '%d' % (stateinit), ' '.join(str('%.10f') % (x) for x in pess)
                        index += 1
                    if tflag == 0:
                        ctraj += 1
                        if tsteps == tsmax:
                            etraj += 1
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)), tcoll)
                os.remove('%s/pes_raw_ensemble.out' % (cwd))
            else:
                print 'Total trajectories:', '%04d' % (ttraj)
                print 'Completed trajectories:', '%04d' % (ctraj)
                print 'Excellent trajectories:', '%04d' % (etraj)
            if errflag == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(tsmax)),tcoll)
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
                    print 'Path %dirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1,int) == true:
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
                    ## Determine if pes files exist and open them ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
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
                            pess = np.float_(pes[line].split())
                            time = np.around(pess[0], decimals = 3)
                            if time != times[index]:
                                print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,dir,times[index])
                                tflag = 1
                                errflag = 1
                                break
                            print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (stateinit), '%d' % (stateinit), ' '.join(str('%.10f') % (x) for x in pess)
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == tsmax:
                                etraj += 1
                    else:
                        print >> error, '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                        errflag = 1
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)), tcoll)
                os.remove('%s/pes_raw_ensemble.out' % (cwd))
            else:
                print 'Total trajectories:', '%04d' % (ttraj)
                print 'Completed trajectories:', '%04d' % (ctraj)
                print 'Excellent trajectories:', '%04d' % (etraj)
            if errflag == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(tsmax)),tcoll)
            else:
                os.remove('%s/pesnact.err' % (cwd))

    ### Non-adiabatic ###
    if boflag == 0:
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
                    print 'path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1,int) == true:
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
                    ## Determine if pes/nact files exist and open them ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    hsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtypeq == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
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
                        if line <= hsteps - 1:
                            nstate = np.int(hops[line].split()[0])
                        if dtypeq in [1,3]:
                            pess = np.float_(pes[line].split())
                            time = np.around(pess[0], decimals = 3)
                            if time != times[index]:
                                print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                tflag = 1
                                errflag = 1
                                break
                        if dtypeq in [2,3]:
                            if line <= hsteps - 1:
                                nacts = np.float_(nact[line - 1].split())[indices]
                                time = np.around(nacts[0], decimals = 3)
                                if time != times[index]:
                                    print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                    tflag = 1
                                    errflag = 1
                                    break
                        if dtypeq == 1:
                            print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                            cstate = nstate
                            index += 1
                            continue
                        if dtypeq == 2:
                            if line <= hsteps - 1:
                                print >> nactall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                                cstate = nstate
                                index += 1
                                continue
                        if dtypeq == 3:
                            print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                            if line <= hsteps - 1:
                                print >> nactall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                        cstate = nstate
                        index += 1
                    if tflag == 0:
                        ctraj += 1
                        if tsteps == tsmax:
                            etraj += 1
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)), tcoll)
                if dtypeq == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
            else:
                print 'Total trajectories:', '%04d' % (ttraj)
                print 'Completed trajectories:', '%04d' % (ctraj)
                print 'Excellent trajectories:', '%04d' % (etraj)
            if errflag == 1:
                print 'One of more trajectories have experienced an error, check pesnact.err.' % (len(str(tsmax)), tcoll)
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
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1,int) == true:
                    dirlist1 = np.array([dirlist1])
                for dir in dirlist1:
                    ## Determine number of completed time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir)
                        sys.exit()
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## Determine if pes/nact files exist ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    htsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtype == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
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
                            if line <= htsteps - 1:
                                nstate = np.int(hops[line].split()[0])
                            if dtypeq in [1,3]:
                                pess = np.float_(pes[line].split())
                                time = np.around(pess[0], decimals = 3)
                                if time != times[index]:
                                    print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                    tflag = 1
                                    errflag = 1
                                    break
                            if dtypeq in [2,3]:
                                if line <= htsteps - 1:
                                    nacts = np.float_(nact[line-1].split())[indices]
                                    time = np.around(nacts[0], decimals = 3)
                                    if time != times[index]:
                                        print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                        tflag = 1
                                        errflag = 1
                                        break
                            if dtypeq == 1:
                                print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                                cstate = nstate
                                index += 1
                                continue
                            if dtypeq == 2:
                                if line <= htsteps - 1:
                                    print >> nactall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                                    cstate = nstate
                                    index += 1
                                    continue
                            if dtypeq == 3:
                                print >> pesall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                                if line <= htsteps - 1:
                                    print >> nactall, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                            cstate = nstate
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == tsmax:
                                etraj += 1
                    else:
                        print >> error, '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                        errflag = 1
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)), tcoll)
                if dtypeq == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (cwd))
                    os.remove('%s/nact_raw_ensemble.out' % (cwd))
            else:
                print 'Total trajectories:', '%04d' % (ttraj)
                print 'Completed trajectories:', '%04d' % (ctraj)
                print 'Excellent trajectories:', '%04d' % (etraj)
            if errflag == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(tsmax)),tcoll)
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
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(dirlist1,int) == true:
                    dirlist1 = np.array([dirlist1])
                for dir in dirlist1:
                    ## determine number of completed time-steps ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,dir)):
                        print 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,dir)
                        sys.exit()
                    data = open('%s/%04d/energy-ev.out' % (NEXMD,dir),'r')
                    data = data.readlines()
                    tsteps = len(data) - 1
                    ## determine if pes/nact files exist ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,dir)):
                        print >> error, 'Path %s%04d/coeff-n.out does not exist.' % (NEXMD,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    hops = open('%s/%04d/coeff-n.out' % (NEXMD,dir),'r')
                    hops = hops.readlines()
                    hsteps = len(hops)
                    if dtypeq == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[0:len(pes):1]
                    if dtypeq == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        nact = open('%s/%04d/nact.out' % (NEXMD,dir),'r')
                        nact = nact.readlines()
                        lines = linenums[1:len(nact) + 1:1]
                    if dtypeq == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/pes.out does not exist.' % (NEXMD,dir)
                            errflag = 1
                            ttraj += 1
                            continue
                        pes = open('%s/%04d/pes.out' % (NEXMD,dir),'r')
                        pes = pes.readlines()
                        lines = linenums[1:len(pes):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,dir)):
                            print >> error, 'Path %s%04d/nact.out does not exist.' % (NEXMD,dir)
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
                            if line <= hsteps - 1:
                                nstate = np.int(hops[line].split()[0])
                            if dtypeq in [1,3]:
                                pess = np.float_(pes[line].split())
                                time = np.around(pess[0], decimals = 3)
                                if time != times[index]:
                                    print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                    tflag = 1
                                    errflag = 1
                                    break
                            if dtypeq in [2,3]:
                                if line <= hsteps - 1:
                                    nacts = np.float_(nact[line - 1].split())[indices]
                                    time = np.around(nacts[0], decimals = 3)
                                    if time != times[index]:
                                        print >> error, 'There is an inconsistency in time-step in %s%04d at %.3f fs.' % (NEXMD,dir,times[index])
                                        tflag = 1
                                        errflag = 1
                                        break
                            if nstate != cstate:
                                if dtypeq == 1:
                                    print >> peshop, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                                    cstate = nstate
                                    index += 1
                                    continue
                                if dtypeq == 2:
                                    if line <= hsteps - 1:
                                        print >> nacthop, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                                        cstate = nstate
                                        index += 1
                                        continue
                                if dtypeq == 3:
                                    print >> peshop, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in pess)
                                    if line <= hsteps - 1:
                                        print >> nacthop, '%s%04d' % (NEXMD,dir), '%d' % (cstate), '%d' % (nstate), ' '.join(str('%.10f') % (x) for x in nacts)
                            cstate = nstate
                            index += 1
                        if tflag == 0:
                            ctraj += 1
                            if tsteps == tsmax:
                                etraj += 1
                    else:
                        print >> error, '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                        errflag = 1
                    print '%s%04d' % (NEXMD,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ttraj += 1
            ## Summary of results ##
            if ctraj == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(tsmax)), tcoll)
                if dtypeq == 1:
                    os.remove('%s/pes_hop_ensemble.out' % (cwd))
                if dtypeq == 2:
                    os.remove('%s/nact_hop_ensemble.out' % (cwd))
                if dtypeq == 3:
                    os.remove('%s/pes_hop_ensemble.out' % (cwd))
                    os.remove('%s/nact_hop_ensemble.out' % (cwd))
            else:
                print 'Total trajectories:', '%04d' % (ttraj)
                print 'Completed trajectories:', '%04d' % (ctraj)
                print 'Excellent trajectories:', '%04d' % (etraj)
            if errflag == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(tsmax)),tcoll)
            else:
                os.remove('%s/pesnact.err' % (cwd))
