#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares input files for an adiabatic simulation,   |
| where geometries are taken from a separate NEXMD simulation       |
| located in another directory.                                     |
|                                                                   |
| This function searches for geometries at a user-defined time in   |
| directory A and generates input files for an adiabatic simulation |
| in directory B.  Typically, directory A is a non-adiabatic        |
| simulation.  In directory B, there must be a 'header' file,       |
| with all inputs set except for 'rand_seed' and                    |
| 'nucl_coord_veloc'.  Since B is an adiabatic simulation,          |
| 'exc_state_init' must be specified in 'header'.  It is also       |
| important to have the 'bo_dynamics_flag' set to '1', indicating   |
| adiabatic dynamics.  The user may choose to use a pre-generated   |
| list of random seeds.  The 'rseedslist' list will be generated in |
| the directory of the adiabatic simulation.  Any problems that may |
| occur during generation of input files will be stated in the      |
| error file, 'newsim.err'.                                         |
|                                                                   |
| NOTE: All NEXMD folders and rseedslists, inside directory B, will |
| be deleted if this function is completely executed!               |
|___________________________________________________________________|

'''


import numpy as np
import random
import os
import sys
import glob
import shutil
import fileinput

cwd = os.getcwd()

def newsim(header):

    print 'Preparing input files for a new simulation with geometries taken from another simulation.'

    ## Directory names ##
    gsdir = raw_input('Ground-state dynamics directory: ')
    if not os.path.exists(gsdir):
        print 'Path %s does not exist.' % (gsdir)
        sys.exit()
    olddir = raw_input('NEXMD directory where geometries should be taken from: ')
    if not os.path.exists(olddir):
        print 'Path %s does not exist.' % (olddir)
        sys.exit()
    newdir = raw_input('NEXMD directory for new simulation: ')
    if not os.path.exists(newdir):
        print 'Path %s does not exist.' % (newdir)
        sys.exit()

    ## Delete previous NEXMD folders and rseedslist ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (newdir))
    NEXMDs.sort()
    if len(NEXMDs) != 0:
        contq = input('** WARNING ** All NEXMD folders and rseedslist inside %s will be deleted!\nContinue? Answer yes [1] or no [0]: ' % (newdir))
        if contq not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if contq == 0:
            sys.exit()
    for NEXMD in NEXMDs:
        print 'Deleting %s.' % (NEXMD)
        shutil.rmtree(NEXMD)
    rseedslist = glob.glob('%s/rseedslist*' % (newdir))
    for rseeds in rseedslist:
        os.remove(rseeds)
        print 'Deleting %s.' % (rseeds)

    ## Information from old header (completed NEXMD geometries) ##
    if not os.path.exists('%s/header' % (olddir)):
        print 'Path %s/header does not exist.' % (olddir)
        sys.exit()
    old_header = header('%s/header' % (oldir))

    ## Check output data ##
    if old_header.out_data_steps == 0:
        print 'No data have been printed to files because out_data_steps = 0 in %s/header.' % (olddir)
        sys.exit()
    if old_header.out_coords_steps == 0:
        print 'No coordinates have been printed to coords.xyz because out_coords_steps = 0 in %/header.' % (olddir)
        sys.exit()

    ## Tell user available information from header regarding coordinates ##
    cprint = old_header.time_step*old_header.out_data_steps*old_header.out_coords_steps
    print 'The coordinates in %s began at %.2f fs and were printed to coords.xyz every %.2f fs.' % (olddir,old_header.time_init,cprint)

    ## Information from new header (desired NEXMD simulations) ##
    if not os.path.exists('%s/header' % (newdir)):
        print 'Path %s/header does not exist.' % (newdir)
        sys.exit()
    new_header = header('%s/header' % (newdir))

    ## Check running dynamics ##
    if new_header.n_class_steps <= 0:
        print 'Must change n_class_steps in %s/header to greater than 0 to run dynamics.' % (newdir)
        sys.exit()

    ## Extract atomic numbers ##
    if not os.path.exists('%s/restart.out' % (gsdir)):
        print 'Path %s/restart.out does not exist.' % (gsdir)
        sys.exit()
    anum = open('%s/restart.out' % (gsdir),'r')
    anum = anum.readlines()
    top = None
    bottom = None
    index = 0
    for line in anum:
        if '$coord' in line:
            top = index
        if '$endcoord' in line:
            bottom = index
            break
        index += 1
    if isinstance(top, int) == True and isinstance(bottom, int) == True:
        anum = [ line.split()[0] for line in anum[top+1:bottom:1] ]
    else:
        print 'There is a problem with %s/restart.out.' % (gsdir)
        sys.exit()

    ## Choose random seeds ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (olddir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (olddir)
        sys.exit()
    with open('%s/totdirlist' % (newdir),'w') as data:
        for NEXMD in NEXMDs:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            inputdata = fileinput.input('%s/dirlist1' % (NEXMD))
            data.writelines(inputdata)
    dirlist1 = np.int_(np.genfromtxt('%s/totdirlist' % (newdir)))
    if isinstance(dirlist1,int) == True:
        dirlist1 = np.array([dirlist1])
    os.remove('%s/totdirlist' % (newdir))
    ntraj = len(dirlist1)
    randq = input('New random seeds? Answer yes [1] or no [0]: ')
    if randq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if randq == 1:
        rseeds = random.sample(np.arange(1,1000001), ntraj)
    else:
        rseeds = raw_input('Path to random-seeds list: ')
        if not os.path.exists(rseeds):
            print 'Path %s does not exist.' % (rseeds)
            sys.exit()
        rseeds = np.int_(np.genfromtxt('%s' % (rseeds)))
        if isinstance(rseeds,int) == True:
            rseeds = np.array([rseeds])
        lenrseeds = len(rseeds)
        if lenrseeds < ntraj:
            print 'Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories requested is %d.' % (lenrseeds,ntraj)
            sys.exit()
    for rseed in rseeds:
        if rseed < 0:
            print 'A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (rseed)
            sys.exit()
    rseeds = np.int_(rseeds)

    ## Choose time at which geometries should be taken from the old simulation ##
    startq = input('At what time, in femtoseconds, should geometries be taken from %s? ' % (olddir))
    if isinstance(startq, int) == False and isinstance(startq, float) == False:
        print 'Time must be integer or float'
        sys.exit()
    if startq < 0:
        print 'Time must be integer or float greater than zero.'
        sys.exit()
    startq = np.float(startq)
    nsteps = 0
    while nsteps*cprint + old_header.time_init <= startq:
        nsteps += 1
    if (nsteps - 1)*cprint + old_header.time_init != startq:
        print 'The time, %.2f, minus the initial time, %.2f, is not divisible by %.2f.' % (startq,old_header.time_init,cprint)
        sys.exit()

    ## Find geometries ##
    print 'Finding coordinates and velocities from %s.  please wait ...' % (olddir)
    tarrayc = np.array([])
    tarrayv = np.array([])
    error = open('%s/newsim.err' % (cwd),'w')
    errflag = 0
    terrflag = 0
    traj = 0
    NEXMDs = glob.glob('%s/NEXMD*/' % (olddir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (olddir)
        sys.exit()
    for NEXMD in NEXMDs:
        dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,dir)):
                print >> error, '%s/%04d/coords.xyz' % (NEXMD,dir), 'does not exist'
                errflag = 1
                terrflag = 1
            if not os.path.exists('%s/%04d/velocity.out' % (NEXMD,dir)):
                print >> error, '%s/%04d/velocity.out' % (NEXMD,dir), 'does not exist'
                errflag = 1
                terrflag = 1
            if errflag == 1:
                tarrayc = np.append(tarrayc, np.array([-9999,-9999]))
                tarrayv = np.append(tarrayv, np.array([-9999,-9999]))
                errflag = 0
                traj += 1
                continue
            datac = open('%s/%04d/coords.xyz' % (NEXMD,dir),'r')
            datac = datac.readlines()
            lenc = len(datac)
            datav = open('%s/%04d/velocity.out' % (NEXMD,dir),'r')
            datav = datav.readlines()
            lenv = len(datav)
            ncoords = 0
            index = 0
            arrayc = np.array([])
            for line in datac:
                if 'time' in line:
                    if ncoords == 0:
                        dtinit = np.around(np.float(line.split()[-1]), decimals = 3)
                        if dtinit != old_header.time_init:
                            print >> error, 'Initial time in %s%04d/coords.xyz does not match initial time in %s/header.' (NEXMD,dir,olddir)
                            errflag = 1
                            terrflag = 1
                            break
                    else:
                        time = np.around(np.float(line.split()[-1]), decimals = 3)
                        if time == startq:
                            tindex = index
                    ncoords += 1
                    arrayc = np.append(arrayc, index)
                index += 1
            if errflag == 1:
                tarrayc = np.append(tarrayc, np.array([-9999,-9999]))
                tarrayv = np.append(tarrayv, np.array([-9999,-9999]))
                errflag = 0
                traj += 1
                continue
            arrayc = np.append(arrayc, lenc + 1)
            cindex = np.where(arrayc == tindex)
            if len(cindex[0]) == 0:
                print >> error, 'Coordinates at %.2f fs in %s%04d/coords.xyz were not found.' % (startq,NEXMD,dir)
                tarrayc = np.append(tarrayc, np.array([-9999,-9999]))
                tarrayv = np.append(tarrayv, np.array([-9999,-9999]))
                terrflag = 1
                traj += 1
                continue
            arrayc = np.int_(arrayc[[cindex[0][0], cindex[0][0] + 1]])
            index = 0
            arrayv = np.array([0])
            for line in datav:
                if 'time' in line:
                    time = np.around(np.float(line.split()[-1]), decimals = 3)
                    if time == startq:
                        tindex = index
                    arrayv = np.append(arrayv, index)
                index += 1
            arrayv = np.append(arrayv, lenv)
            vindex = np.where(arrayv == tindex)
            if len(vindex[0]) == 0:
                print >> error, 'Velocities at %.2f fs in %s%04d/velocity.out were not found.' % (startq,NEXMD,dir)
                tarrayc = np.append(tarrayc, np.array([-9999,-9999]))
                tarrayv = np.append(tarrayv, np.array([-9999,-9999]))
                terrflag = 1
                traj += 1
                continue
            arrayv = np.int_(arrayv[[vindex[0][0], vindex[0][0] + 1]])
            tarrayc = np.append(tarrayc, arrayc)
            tarrayv = np.append(tarrayv, arrayv)
            traj += 1
    if terrflag == 1:
        print 'One or more trajectories of the previous simulation cannot be used for in the current simulation, check newsim.err.'
        contq = input('Continue preparing input files? Answer yes [1] or no [0]: ')
        if contq not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if contq == 0:
            sys.exit()
    tarrayc = np.int_(np.split(tarrayc, traj))
    tarrayv = np.int_(np.split(tarrayv, traj))

    ## Prepare NEXMD input files ##
    errflag = 0
    traj = 0
    for NEXMD in np.arange(1, len(NEXMDs) + 1):
        os.makedirs('%s/NEXMD%d' % (newdir,NEXMD))
        dirlist1 = np.int_(np.genfromtxt('%s/NEXMD%d/dirlist1' % (olddir,NEXMD)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        dirlist = open('%s/NEXMD%d/dirlist' % (newdir,NEXMD),'w')
        for dir in dirlist1:
            os.makedirs('%s/NEXMD%d/%04d' % (newdir,NEXMD,dir))
            if -9999 in tarrayc[traj]:
                print >> error, 'The input file for %s/NEXMD%d/%04d cannot be generated.' % (olddir,NEXMD,dir)
                errflag = 1
                traj += 1
                continue
            datac = open('%s/NEXMD%d/%04d/coords.xyz' % (olddir,NEXMD,dir),'r')
            datac = datac.readlines()
            datav = open('%s/NEXMD%d/%04d/velocity.out' % (olddir,NEXMD,dir),'r')
            datav = datav.readlines()
            coords = datac[tarrayc[traj][0] + 1:tarrayc[traj][1] - 1:1]
            velocs = datav[tarrayv[traj][0] + 2:tarrayv[traj][1] - 1:1]
            inputfile = open('%s/NEXMD%d/%04d/input.ceon' % (newdir,NEXMD,dir),'w')
            for line in new_header.file:
                if 'rnd_seed' in line:
                    inputfile.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseeds[traj]))
                else:
                    if 'nucl_coord_veloc' in line:
                        inputfile.write('&coord\n')
                        aindex = 0
                        for line in coords:
                            val = line.split()
                            inputfile.write('{:>6}  {:>12}  {:>12}  {:>12}'.format(anum[aindex],val[1],val[2],val[3]))
                            inputfile.write('\n')
                            aindex += 1
                        inputfile.write('&endcoord\n\n&veloc\n')
                        for line in velocs:
                            inputfile.write(line)
                        inputfile.write('&endveloc')
                    else:
                        inputfile.write(line)
            print >> dirlist, '%04d' % (dir)
            print '%s/NEXMD%d/%04d' % (newdir,NEXMD,dir)
            traj += 1
        dirlist.close()
        shutil.copyfile('%s/NEXMD%d/dirlist' % (newdir,NEXMD), '%s/NEXMD%d/dirlist1' % (newdir,NEXMD))
    np.savetxt('%s/rseedslist' % (newdir), np.transpose(rseeds[0:traj:1]))
    if errflag == 1:
        print 'One or more trajectories cannot be prepared, check newsim.err.'
        sys.exit()
    else:
        os.remove('%s/newsim.err' % (cwd))
