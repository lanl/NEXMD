#/usr/bin/python

'''

This function prepares input files for non-adiabatic
excited-state molecular dynamics (NEXMD).

A general header called 'header' must be in the NEXMD directory
(e.g. NEXMD) and must have all inputs set except for:

1) Random seed (rnd_seed_flag)
2) Initial excited state (exc_state_init_flag)
3) Initial nuclear coordinates and velocities (nucl_coord_veloc_flag)
4) Initial quantum amplitudes and phase (quant_amp_phase_flag)

In parentheses shown above, are flags in 'header' that label
these inputs. This function finds these labels and fills them in
accordingly. A 'ceo.err' file will be generated in case there
exists incomplete single-point calculations.

NOTE: All NEXMD folders and rseedslists, inside the NEXMD
directory, will be deleted if this function is completely
executed!

'''

import numpy as np
import random
import os
import sys
import shutil
import glob
import fileinput

cwd = os.getcwd()

def nexmd(header):
    
    print('Preparing input files for NEXMD.')
    
    ## Directory names ##
    gsdir = raw_input('Ground-state dynamics directory: ')
    if not os.path.exists(gsdir):
        print('Path %s does not exist.' % (gsdir))
        sys.exit()
    outdir = raw_input('Output directory [e.g. NEXMD]: ')
    if not os.path.exists(outdir):
        print('Path %s does not exist.' % (outdir))
        sys.exit()

    ## Delete previous NEXMD folders and rseedslists ##
    NEXMDs = glob.glob('%s/NEXMD*/' % (outdir))
    NEXMDs.sort()
    if len(NEXMDs) != 0:
        contq = input('** WARNING ** All NEXMD folders and rseedslists inside %s will be deleted!\ndo you want to continue? Answer yes [1] or no [0]: ' % (outdir))
        if contq not in [1,0]:
            print('answer must be 1 or 0.')
            sys.exit()
        if contq == 0:
            sys.exit()
    for NEXMD in NEXMDs:
        print('Deleting', '%s' % (NEXMD))
        shutil.rmtree(NEXMD)
    rseedslist = glob.glob('%s/rseedslist*' % (outdir))
    for rseeds in rseedslist:
        os.remove(rseeds)
        print('Deleting', '%s' % (rseeds))

    ## Information from header ##
    if not os.path.exists('%s/header' % (outdir)):
        print('Path %s/header does not exist.' % (outdir))
        sys.exit()
    header = header('%s/header' % (outdir))

    ## Check running dynamics ##
    if header.n_class_steps <= 0:
        print('Must change n_class_steps in %s/header to greater than 0 for dynamics.' % (outdir))
        sys.exit()

    ## Check type of initial excitation ##
    #waves = 1 # only for old code
    #wavec = 1 # only for old code
    ## possible combinations ##
    ## state set, coefficients set          (state_set = 1, coeff_set = 1) - single state
    ## state set, coefficients not set      (state_set = 1, coeff_set = 0) - single state
    ## state not set, coefficients set      (state_set = 0, coeff_set = 1) - ridiculous (exit system)
    ## state not set, coefficients not set  (state_set = 0, coeff_set = 0) - photoexcited wavepacket
    try:
        header.exc_state_init
        state_set = 1
    except AttributeError:
        state_set = 0
    try:
        header.quant_amp_phase_flag
        coeff_set = 0
    except AttributeError:
        coeff_set = 1
    if state_set == 1 and coeff_set == 1 or state_set == 1 and coeff_set == 0:
        print('All trajectories will begin on state %d.' % (header.exc_state_init))
    if state_set == 0 and coeff_set == 1:
        print('There is an inconsistency in header.\nInput exc_state_init is not set, while coefficients are set.')
        sys.exit()
    if state_set == 0 and coeff_set == 0:
        print('Initial excited states will model a photoexcited wavepacket according to the optical spectrum.')
        spdir = raw_input('Single-point calculations directory: ')
        if not os.path.exists(spdir):
            print('Path %s does not exist.' % (spdir))
            sys.exit()
    
    ## Find geometries ##
    if not os.path.exists('%s/coords.xyz' % (gsdir)):
        print('Path %s/coords.xyz does not exist.' % (gsdir))
        sys.exit()
    datac = open('%s/coords.xyz' % (gsdir),'r')
    datac = datac.readlines()
    lenc = len(datac)
    if not os.path.exists('%s/velocity.out' % (gsdir)):
        print('path %s/velocity.out does not exist.' % (gsdir))
        sys.exit()
    datav = open('%s/velocity.out' % (gsdir),'r')
    datav = datav.readlines()
    lenv = len(datav)
    ncoords = 0
    index = 0
    arrayc = np.array([])
    for line in datac:
        if 'time' in line:
            if ncoords == 0:
                tinit = np.float(line.split()[-1])
            else:
                time = np.float(line.split()[-1])
            ncoords += 1
            arrayc = np.append(arrayc,index)
        index += 1
    arrayc = np.append(arrayc,lenc + 1)
    arrayc = np.int_(arrayc)
    if ncoords == 0:
        print('No coordinates were found.')
        sys.exit()
    if ncoords == 1:
        print('Only initial coordinates, at %.2f fs, were found.' % (tinit))
        sys.exit()
    if ncoords > 1:
        index = 0
        arrayv = np.array([])
        for line in datav:
            if 'time' in line:
                arrayv = np.append(arrayv,index)
            index += 1
        arrayv = np.append(arrayv,lenv)
        arrayv = np.int_(arrayv)
    tinc = time/(ncoords - 1)
    if state_set == 1 and coeff_set == 1 or state_set == 1 and coeff_set == 0:
        print('A total of %d coordinate files, ranging from %.2f to %.2f fs in increments of %.2f fs, were found.' % (ncoords,tinit,time,tinc))
    if state_set == 0 and coeff_set == 0:
        print('A total of %d coordinate files, ranging from %.2f to %.2f fs in increments of %.2f fs, were found.\nNote: only coordinate files used for single-point calculations can be used for NEXMD.' % (ncoords,tinit,time,tinc))

    ## Choose geometries for a photoexcited wavepacket ##
    if state_set == 0 and coeff_set == 0:
        NEXMDs = glob.glob('%s/NEXMD*/' % (spdir))
        NEXMDs.sort()
        if len(NEXMDs) == 0:
            print('There are no NEXMD folders in %s.' % (spdir))
            sys.exit()
        with open('%s/totdirlist' % (outdir),'w') as data:
            for NEXMD in NEXMDs:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print('Path %sdirlist1 does not exist.' % (NEXMD))
                    sys.exit()
                inputdata = fileinput.input('%s/dirlist1' % (NEXMD))
                data.writelines(inputdata)
        dirlist1 = np.int_(np.genfromtxt('%s/totdirlist' % (outdir)))
        if isinstance(dirlist1,int) == True:
            dirlist1 = np.array([dirlist1])
        os.remove('%s/totdirlist' % (outdir))
        ntraj = len(dirlist1)
        ntrajq = input('How many trajectories for NEXMD? Enter a number no greater than %d: ' % (ntraj))
        if isinstance(ntrajq, int) == False:
            print('Number of trajectories must be integer.')
            sys.exit()
        if ntrajq == 0:
            print('Number of trajectories must be positive integer.')
            sys.exit()
        if np.abs(ntrajq) > ntraj:
            print('Number of trajectories must be less than or equal to %d.' % (ntraj))
            sys.exit()
        if ntrajq < 0:
            ntraj = np.abs(ntrajq)
            coordsq = input('You have requested %d randomly-selected coordinate files in the range %d to %d.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,dirlist[0],dirlist[-1]))
        else:
            interval = np.int(np.ceil(ntraj/np.float(ntrajq)))
            if interval*ntrajq > ntraj:
                interval = interval - 1
                ntraj = len(dirlist1[0::interval])
            else:
                ntraj = ntrajq
            coordsq = input('You have requested %d evenly-spaced coordinate files in the range %d to %d for NEXMD.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,dirlist1[0],dirlist1[0::interval][-1]))
        if coordsq not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if coordsq == 0:
            sys.exit()

    ## Choose geometries for a single excited state ##
    if state_set == 1 and coeff_set == 1 or state_set == 1 and coeff_set == 0:
        coords = input('Enter requested range of the ground-state sampling by coordinate files and the number of trajectories.\nInput an array of the form [start, end, number]: ')
        if not isinstance(coords,list):
            print('Input must be an array of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.')
            sys.exit()
        if len(coords) != 3:
            print('Input must be an array with three elements of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.')
            sys.exit()
        index = 0
        for i in coords:
            if isinstance(i, int) == False:
                print('Element number %d of input array must be integer.\nUser inputted [%s, %s, %s], which is not allowed.' % (index + 1,coords[0],coords[1],coords[2]))
                sys.exit()
            if index in [0,1]:
                if i not in np.arange(ncoords):
                    print('Element number %d of input array must be integer between 0 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (index + 1,ncoords - 1,coords[0],coords[1],coords[2]))
                    sys.exit()
            else:
                if i not in np.arange(1,ncoords + 1) and i not in -np.arange(1,ncoords + 1):
                    print('Element number %d of input array must be integer between 1 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (index + 1,ncoords,coords[0],coords[1],coords[2]))
                    sys.exit()
            index += 1
        if coords[0] > coords[1]:
            print('Second element of input array must be greater than first element.\nUser inputted [%d, %d, %d], which is not allowed.' % (coords[0],coords[1],coords[2]))
            sys.exit()
        if (coords[1] - coords[0]) + 1 < np.abs(coords[2]):
            print('Number of coordinate files requested must be less than or equal to the number of coordinate files in the sample.\nUser inputted [%d, %d, %d], which is not allowed.' % (coords[0],coords[1],coords[2]))
            sys.exit()
        if coords[2] < 0:
            ntraj = np.abs(coords[2])
            coordsq = input('You have requested %d randomly-selected coordinate files in the range %d to %d.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,coords[0],coords[1]))
        else:
            interval = np.int(np.ceil((coords[1] - coords[0] + 1)/np.float(coords[2])))
            if interval*coords[2] > (coords[1] - coords[0] + 1):
                interval = interval - 1
                ntraj = len(np.arange(coords[0],coords[1] + 1,interval))
            else:
                ntraj = coords[2]
            coordsq = input('You have requested %d evenly-spaced coordinate files in the range %d to %d.\nContinue? Answer yes [1] or no [0]: ' % (ntraj,coords[0],coords[1]))
        if coordsq not in [1,0]:
            print('Answer must be 1 or 0.')
            sys.exit()
        if coordsq == 0:
            sys.exit()

    ## Split geometries ##
    split = input('Number of trajectories per NEXMD folder: ')
    if isinstance(split, int) == False:
        print('Number of trajectories per NEXMD folder must be integer.')
        sys.exit()
    if split <= 0:
        print('Number of trajectories per NEXMD folder must be integer greater than zero.')
        sys.exit()
    dirsplit = split*np.arange(1,np.ceil(np.float(ntraj)/split)+1)
    if state_set == 1 and coeff_set == 1 or state_set == 1 and coeff_set == 0:
        dirsplit[-1] = coords[1] + 1
        if coords[2] < 0:
            dirsplit = np.split(np.sort(random.sample(np.arange(coords[0],coords[1]+1),ntraj)),dirsplit)
        else:
            dirsplit = np.split(np.arange(coords[0],coords[1]+1,interval),dirsplit)
    if state_set == 0 and coeff_set == 0:
        dirsplit[-1] = dirlist1[-1]
        if ntrajq < 0:
            dirsplit = np.split(np.sort(random.sample(dirlist1,ntraj)),dirsplit)
        else:
            dirsplit = np.split(dirlist1[0::interval],dirsplit)

    ## Extract atomic numbers ##
    if not os.path.exists('%s/restart.out' % (gsdir)):
        print('Path %s/restart.out does not exist.' % (gsdir))
        sys.exit()
    anum = open('%s/restart.out' % (gsdir),'r')
    anum = anum.readlines()
    top = None
    bottom = None
    index = 0
    for line in anum:
        if '$COORD' in line:
            top = index
        if '$ENDCOORD' in line:
            bottom = index
            break
        index += 1
    if isinstance(top, int) == True and isinstance(bottom, int) == True:
        anum = [ line.split()[0] for line in anum[top+1:bottom:1] ]
    else:
        print('There is a problem with %s/restart.out.' % (gsdir))
        sys.exit()

    ## Choose random seeds ##
    randq = input('New random seeds? answer yes [1] or no [0]: ')
    if randq not in [1,0]:
        print('Answer must be 1 or 0.')
        sys.exit()
    if randq == 1:
        rseeds = random.sample(np.arange(1,1000001), ntraj)
    else:
        rseeds = raw_input('path to random-seeds list: ')
        if not os.path.exists(rseeds):
            print('Path %s does not exist.' % (rseeds))
            sys.exit()
        rseeds = np.int_(np.genfromtxt('%s' % (rseeds)))
        if isinstance(rseeds,int) == True:
            rseeds = np.array([rseeds])
        lenrseeds = len(rseeds)
        if lenrseeds < ntraj:
            print('Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories requested is %d.' % (lenrseeds,ntraj))
            sys.exit()
    for rseed in rseeds:
        if rseed < 0:
            print('A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (rseed))
            sys.exit()
    rseeds = np.int_(rseeds)

    ## Prepare NEXMD input files with a photoexcited wavepacket ##
    if state_set == 0 and coeff_set == 0:
        spNEXMDs = glob.glob('%s/NEXMD*/' % (spdir))
        spNEXMDs.sort()
        error = open('%s/ceo.err' % (cwd),'w')
        stype = input('Spectral lineshape? Answer Gaussian [0] or Lorentzian [1]: ')
        if stype not in [0,1]:
            print('Answer must be 0 or 1.')
            sys.exit()
        excen = input('Laser excitation energy in eV: ')
        if isinstance(excen, int) == False and isinstance(excen, float) == False:
            print('Excitation energy must be integer or float.')
            sys.exit()
        if stype == 0:
            specb = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
        else:
            specb = input('Spectral broadening (i.e. Lorentzian fwhm) in eV [e.g. 0.36]: ')
        if isinstance(specb, int) == False and isinstance(specb, float) == False:
            print('Spectral broadening must be integer or float.')
            sys.exit()
        traj = 0
        index = 0
        ceoflag = 0
        for NEXMD in np.arange(1,np.ceil(np.float(ntraj)/split)+1):
            os.makedirs('%s/NEXMD%d' % (outdir,NEXMD))
            dirlist = open('%s/NEXMD%d/dirlist' % (outdir,NEXMD),'w')
            for dir in dirsplit[index]:
                os.makedirs('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
                for spNEXMD in spNEXMDs:
                    if os.path.exists('%s/%04d/ceo.out' % (spNEXMD,dir)):
                        data = np.genfromtxt('%s/%04d/ceo.out' % (spNEXMD,dir))
                        iceoflag = 0
                        break
                    else:
                        iceoflag = 1
                if iceoflag == 1:
                    error.write( '%s%04d/ceo.out does not exist' % (spNEXMD,dir))
                    rseeds[traj] = -123456789
                    ceoflag = 1
                    iceoflag = 0
                    traj += 1
                    continue
                lendata = len(data)
                qpop = np.zeros(lendata)
                oindex = 0
                if stype == 0:
                    for state, energy, osx, osy, osz, oscstren in data:
                        qpop[oindex] = oscstren*np.exp(-(energy - excen)**(2.0)/(2.0*specb**(2.0)))/np.sqrt(2.0*np.pi*specb**(2.0))
                        oindex += 1
                else:
                    for state, energy, osx, osy, osz, oscstren in data:
                        qpop[oindex] = oscstren/(1.0 + ((energy - excen)**2.0)/(specb/2.0)**(2.0))/(specb*np.pi/2.0)
                        oindex += 1
                qpop = qpop/np.sum(qpop)
                state = np.searchsorted(np.cumsum(qpop),np.random.uniform()) + 1
                qpop = np.zeros(lendata)
                qpop[state-1] = 1.0
                coords = datac[arrayc[dir]+1:arrayc[dir+1]-1:1]
                velocs = datav[arrayv[dir]+2:arrayv[dir+1]-1:1]
                inputfile = open('%s/NEXMD%d/%04d/input.ceon' % (outdir,NEXMD,dir),'w')
                for line in header.file:
                    if 'rnd_seed_flag' in line:
                        inputfile.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseeds[traj]))
                    else:
                        if 'exc_state_init_flag' in line:
                            inputfile.write('   exc_state_init=%d, ! initial excited state (0 - ground state) [0]\n' % (state))
                        else:
                            if 'nucl_coord_veloc_flag' in line:
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
                                inputfile.write('&endveloc\n')
                            else:
                                if 'quant_amp_phase_flag' in line:
                                    inputfile.write('&coeff\n')
                                    for line in qpop:
                                        inputfile.write('  %.3f  %.3f\n' % (line,0.0))
                                    inputfile.write('&endcoeff\n')
                                else:
                                    inputfile.write(line)
                dirlist.write( '%04d' % (dir))
                print('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
                traj += 1
            dirlist.close()
            shutil.copyfile('%s/NEXMD%d/dirlist' % (outdir,NEXMD), '%s/NEXMD%d/dirlist1' % (outdir,NEXMD))
            index += 1
        np.savetxt('%s/rseedslist' % (outdir), np.transpose(rseeds[0:traj:1]))
        if ceoflag == 1:
            print('One or more NEXMD trajectories cannot be prepared, check ceo.err.')
            sys.exit()
        else:
            os.remove('%s/ceo.err' % (cwd))

    ## Prepare NEXMD input files with a single excited state ##
    if state_set == 1 and coeff_set == 1 or state_set == 1 and coeff_set == 0:
        qpop = np.zeros(header.exc_state_init)
        qpop[header.exc_state_init - 1] = 1.0
        traj = 0
        index = 0
        for NEXMD in np.arange(1,np.ceil(np.float(ntraj)/split)+1):
            os.makedirs('%s/NEXMD%d' % (outdir,NEXMD))
            dirlist = open('%s/NEXMD%d/dirlist' % (outdir,NEXMD),'w')
            for dir in dirsplit[index]:
                os.makedirs('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
                coords = datac[arrayc[dir]+1:arrayc[dir+1]-1:1]
                velocs = datav[arrayv[dir]+2:arrayv[dir+1]-1:1]
                inputfile = open('%s/NEXMD%d/%04d/input.ceon' % (outdir,NEXMD,dir),'w')
                for line in header.file:
                    if 'rnd_seed_flag' in line:
                        inputfile.write('   rnd_seed=%d, ! seed for the random number generator\n' % (rseeds[traj])) ## for new code (NEXMD)
                        #inputfile.write('%d ! seed for the random number generator\n' % (rseeds[traj])) ## for old code (NAESMD)
                    else:
                        if 'nucl_coord_veloc_flag' in line:
                            inputfile.write('&coord\n')  ## for new code (NEXMD)
                            #inputfile.write('$COORD\n')  ## for old code (NAESMD)
                            aindex = 0
                            for line in coords:
                                val = line.split()
                                inputfile.write('{:>6}  {:>12}  {:>12}  {:>12}'.format(anum[aindex],val[1],val[2],val[3]))
                                inputfile.write('\n')
                                aindex += 1
                            inputfile.write('&endcoord\n\n&veloc\n')  ## for new code (NEXMD)
                            #inputfile.write('$ENDCOORD\n\n$VELOC\n') ## for old code (NAESMD)
                            for line in velocs:
                                inputfile.write(line)
                            inputfile.write('&endveloc\n')
                            #inputfile.write('$ENDVELOC\n') ## for old code (NAESMD)
                        else:
                            if 'quant_amp_phase_flag' in line:
                                inputfile.write('&coeff\n')
                                for line in qpop:
                                    inputfile.write('  %.3f  %.3f\n' % (line,0.0))
                                inputfile.write('&endcoeff\n')
                            else:
                                inputfile.write(line)
                dirlist.write( '%04d' % (dir))
                print('%s/NEXMD%d/%04d' % (outdir,NEXMD,dir))
                traj += 1
            dirlist.close()
            shutil.copyfile('%s/NEXMD%d/dirlist' % (outdir,NEXMD), '%s/NEXMD%d/dirlist1' % (outdir,NEXMD))
            index += 1
        np.savetxt('%s/rseedslist' % (outdir), np.transpose(rseeds[0:traj:1]))
