#/usr/bin/python

'''
 _____________________________________________________________________
|                                                                     |
| This function calculates a user-specified dihedral angle.  The user |
| must supply the line number of four atoms, where the first atom     |
| between '&coord' and '&endcoord' in 'input.ceon' is always labeled  |
| with line number = 0.  The user specifies the three atoms that make |
| the first plane and another three atoms that make the second plane. | 
| Since there are four unique atoms, two atoms will be shared between |
| both planes.  These two atoms should be in input[1] and input[2] of |
| the input array, where input[1] is the second atom defining the     |
| first plane (first atom defining the second plane) and input[2] is  |
| the second atom defining the second plane (last atom defining the   |
| first plane).  The output is 'dihedral.out', where first column     |
| is time in femtoseconds and second column is dihedral angle in      |
| degrees.  The user must specify whether the calculation is for a    |
| single trajectory or an ensemble of trajectories.  For the latter,  |
| the second column in 'dihedral.out' is an average of the dihedral   |
| angle as a function of time.  There will also be a third column,    |
| showing standard deviation of the dihedral angle as a function of   |
| time. An error file called 'dihedral.err' is generated if an error  |
| occurs, such as nonexistent or incomplete files, etc.               |
|_____________________________________________________________________|

'''

import numpy as np
import os
import sys
import glob
import fileinput
import math

CWD = os.getcwd()

def DIHEDRAL():

    print 'Calculating a dihedral angle as a function of time.'

## DIRECTORY NAMES ##
    DYNQ = input('Calculate dihedral along one trajectory or an ensemble of trajectories? Answer ONE [1] or ENSEMBLE [0]: ')
    if DYNQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if DYNQ == 0:
        NEXMDIR = raw_input('Ensemble directory [e.g. nexmd]: ')
    else:
        NEXMDIR = raw_input('Single trajectory directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()

## USER-DEFINED LENGTH OF ANALYSIS AND INITIALIZE ARRAYS ##
    if DYNQ == 1:
        if not os.path.exists('%s/input.ceon' % (NEXMDIR)):
            print 'Path %s/input.ceon does not exist.' % (NEXMDIR)
            sys.exit()
        HEADER = open('%s/input.ceon' % (NEXMDIR),'r')
        HEADER = HEADER.readlines()
    else:
        if not os.path.exists('%s/header' % (NEXMDIR)):
            print 'Path %s/header does not exist.' % (NEXMDIR)
            sys.exit()
        HEADER = open('%s/header' % (NEXMDIR),'r')
        HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'time_init' in LINE:
            TINITH = np.float(LINE.split()[0][len('time_init='):-1])
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1]) + 1
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
        if 'out_coords_steps' in LINE:
            CDATA = np.int(LINE.split()[0][len('out_coords_steps='):-1])
        if 'natoms' in LINE:
            NATOMS = np.int(LINE.split()[0][len('natoms='):-1])
    TCOLL = input('Calculate dihedral angle up to what time in femtoseconds: ')
    if isinstance(TCOLL, int) == False and isinstance(TCOLL, float) == False:
        print 'Time must be integer or float.'
        sys.exit()
    if TCOLL < 0:
        print 'Time must be integer or float greater than zero.'
        sys.exit()
    TCOLL = np.float(TCOLL)
    if TCOLL > (TSMAX - 1)*DT:
        TCOLL = (TSMAX - 1)*DT
    TSCOL = 0
    while TSCOL*DT*ODATA <= TCOLL:
        TSCOL += 1
    CCOLL = 0
    NUM = 0
    while CCOLL <= TCOLL:
        CCOLL += DT*ODATA*CDATA
        NUM += 1
    TIMES = np.linspace(TINITH, CCOLL - DT*ODATA*CDATA, NUM)
    if DYNQ == 0:
        FDIHEDRAL = np.zeros(len(TIMES) - 1)

## FOUR UNIQUE ATOMS DEFINED BY USER ##
    LINES = input('Input the line numbers labeling the coordinates of the four atoms.\nInput an array of the form [ .., .., .., .. ]: ')
    if isinstance(LINES, list) == False:
        print 'Input must be an array of the form [atom 1, atom2, atom3, atom4], where atom# = line number of atom#.'
        sys.exit()
    if len(LINES) != 4:
        print 'Input must be an array with four elements labeling the line numbers of four atoms.'
        sys.exit()
    INDEX = 0
    for i in LINES:
        if isinstance(i, int) == False:
            print 'Element number %d of input array must be integer.\nUser inputted [%s, %s, %s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1], LINES[2], LINES[3])
            sys.exit()
        if i < 0:
            print 'Element number %d of input array must be a positive integer.\nUser inputted [%s, %s, %s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1], LINES[2], LINES[3])
            sys.exit()
        if i > NATOMS - 1:
            print 'Element number %d of input array must be less than the max number of atoms (-1).\nUser inputted [%s, %s, %s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1], LINES[2], LINES[3])
            sys.exit()
        INDEX += 1
    if len(np.unique(LINES)) != 4:
        print 'All elements of input array must be unique.\nUser inputted [%s, %s, %s, %s], which is not allowed.' % (LINES[0], LINES[1], LINES[2], LINES[3])
        sys.exit()

## CALCULATE DIHEDRAL ANGLE ALONG A SINGLE TRAJECTORY ##
    if DYNQ == 1:
        if not os.path.exists('%s/energy-ev.out' % (NEXMDIR)):
            print 'Path %s/energy-ev.out does not exist.' % (NEXMDIR)
            sys.exit()
        DATA = open('%s/energy-ev.out' % (NEXMDIR),'r')
        DATA = DATA.readlines()
        TSTEPS = len(DATA) - 1
        OUTPUT = open('%s/dihedral.out' % (CWD),'w')
        if TSTEPS >= TSCOL:
            if not os.path.exists('%s/coords.xyz' % (NEXMDIR)):
                print 'Path %s%coords.xyz does not exist.' % (NEXMDIR)
                sys.exit()
            DATA = open('%s/coords.xyz' % (NEXMDIR),'r')
            DATA = DATA.readlines()
            LENC = len(DATA)
            NCOORDS = 0
            CINDEX = 0
            TFLAG1 = 0
            TFLAG2 = 0
            ARRAY = np.array([])
            for LINE in DATA:
                if 'time' in LINE:
                    if NCOORDS == 0:
                        TINIT = np.float(LINE.split()[-1])
                        if TINIT != TINITH:
                            TFLAG1 = 1
                            continue
                    else:
                        TIME = np.around(np.float(LINE.split()[-1]), decimals = 3)
                        if TIME > TCOLL:
                            continue
                        if TIME != TIMES[NCOORDS]:
                            TFLAG2 = 1
                            continue
                    NCOORDS += 1
                    ARRAY = np.append(ARRAY,CINDEX)
                CINDEX += 1
            if TFLAG1 == 1:
                print 'Initial time in %s/coords.xyz does not match time_init in %s/input.ceon.' % (NEXMDIR,NEXMDIR)
                sys.exit()
            if TFLAG2 == 1:
                print 'There is an inconsistency in time-step in %s/coords.xyz.' % (NEXMDIR)
                sys.exit()
            if TIME < TCOLL:
                ARRAY = np.append(ARRAY, LENC + 1)
            ARRAY = np.int_(ARRAY)
        ## CHECKS TO ENSURE DIHEDRAL CALCULATION ##
            if NCOORDS == 0:
                print 'No coordinates were found in %s/coords.xyz' % (NEXMDIR)
                sys.exit()
            if NCOORDS == 1:
                print 'Only initial coordinates, at %.2f fs, were found in %s/coords.xyz.' % (TINIT,NEXMDIR)
                sys.exit()
        ## CALCULATE DIHEDRAL ANGLE ALONG A SINGLE TRAJECTORY ##
            SDIHEDRAL = np.zeros(NCOORDS - 1)
            for NCOORD in np.arange(NCOORDS - 1):
                COORDS = DATA[ARRAY[NCOORD]+1:ARRAY[NCOORD+1]-1:1]
                VEC0 = np.float_(COORDS[LINES[0]].split()[1:])
                VEC1 = np.float_(COORDS[LINES[1]].split()[1:])
                VEC2 = np.float_(COORDS[LINES[2]].split()[1:])
                VEC3 = np.float_(COORDS[LINES[3]].split()[1:])
                A = np.subtract(VEC1, VEC0)
                B = np.subtract(VEC2, VEC1)
                C = np.subtract(VEC3, VEC2)
                N1 = np.cross(A, B)
                N2 = np.cross(B, C)
                UN1 = N1/np.linalg.norm(N1)
                UN2 = N2/np.linalg.norm(N2)
                UB = B/np.linalg.norm(B)
                M1 = np.cross(UN1, UB)
                X = np.dot(UN1, UN2)
                Y = np.dot(M1, UN2)
                SDIHEDRAL[NCOORD] = np.abs(np.degrees(-math.atan2(Y,X)))
            print '%s' % (NEXMDIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
            CTRAJ = 1
            if TSTEPS == TSMAX:
                ETRAJ = 1
            else:
                print '%s' % (NEXMDIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
        TTRAJ = 1
        if CTRAJ == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(TSMAX)),TCOLL)
        else:
            print 'Total Trajectories:', '%04d' % (TTRAJ)
            print 'Completed Trajectories:', '%04d' % (CTRAJ)
            print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            print >> OUTPUT, 'Total Trajectories: ', '%04d' % (TTRAJ)
            print >> OUTPUT, 'Completed Trajectories: ', '%04d' % (CTRAJ)
            print >> OUTPUT, 'Excellent Trajectories: ', '%04d' % (ETRAJ)
            for NCOORD in np.arange(NCOORDS - 1):
                print >> OUTPUT, '%0*.2f' % (len(str((TSMAX))) + 2,DT*ODATA*CDATA*(NCOORD + 1)), '%08.3f' % (SDIHEDRAL[NCOORD])
                    
## CALCULATE DIHEDRAL ANGLE ALONG AN ENSEMBLE OF TRAJECTORIES ##
    else:
        NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
        NEXMDS.sort()
        if len(NEXMDS) == 0:
            print 'There are no NEXMD folders in %s.' % (NEXMDIR)
            sys.exit()
        with open('%s/totdirlist' % (NEXMDIR),'w') as DATA:
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %NEXMDIRlist1 does not exist.' % (NEXMD)
                    sys.exit()
                INPUT = fileinput.input('%s/dirlist1' % (NEXMD))
                DATA.writelines(INPUT)
        DIRLIST1 = np.int_(np.genfromtxt('%s/totdirlist' % (NEXMDIR)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        os.remove('%s/totdirlist' % (NEXMDIR))
        OUTPUT = open('%s/dihedral.out' % (CWD),'w')
        ERROR = open('%s/dihedral.err' % (CWD),'w')
        TTRAJ = 0
        CTRAJ = 0
        ETRAJ = 0
        EDIHEDRAL = np.zeros((len(TIMES) - 1, len(DIRLIST1)))
        ERRFLAG = 0
        for NEXMD in NEXMDS:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %NEXMDIRlist1 does not exist.' % (NEXMD)
                sys.exit()
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1, int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                TSTEPS = len(DATA) - 1
                if TSTEPS >= TSCOL:
                    if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/coords.xyz' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                ## GENERATE ARRAY WITH INDICES OF THE COORDINATE BLOCKS ALONG A SINGLE TRAJECTORY ##
                    DATA = open('%s/%04d/coords.xyz' % (NEXMD,DIR),'r')
                    DATA = DATA.readlines()
                    LENC = len(DATA)
                    NCOORDS = 0
                    CINDEX = 0
                    TFLAG1 = 0
                    TFLAG2 = 0
                    ARRAY = np.array([])
                    for LINE in DATA:
                        if 'time' in LINE:
                            if NCOORDS == 0:
                                TINIT = np.float(LINE.split()[-1])
                                if TINIT != TINITH:
                                    TFLAG1 = 1
                                    continue
                            else:
                                TIME = np.around(np.float(LINE.split()[-1]), decimals = 3)
                                if TIME > TCOLL:
                                    continue
                                if TIME != TIMES[NCOORDS]:
                                    TFLAG2 = 1
                                    continue
                            NCOORDS += 1
                            ARRAY = np.append(ARRAY,CINDEX)
                        CINDEX += 1
                    if TFLAG1 == 1:
                        print >> ERROR, 'Initial time in %s%04d/coords.xyz does not match time_init in %s/header.' % (NEXMD,DIR,NEXMDIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    if TFLAG2 == 1:
                        print >> ERROR, 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (NEXMD,DIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    if TIME < TCOLL:
                        ARRAY = np.append(ARRAY, LENC + 1)
                    ARRAY = np.int_(ARRAY)
                ## CHECKS TO ENSURE DIHEDRAL CALCULATION ##
                    if NCOORDS == 0:
                        print >> ERROR, 'No coordinates were found in %s%04d/coords.xyz' % (NEXMD,DIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    if NCOORDS == 1:
                        print >> ERROR, 'Only initial coordinates, at %.2f fs, were found in %s%04d/coords.xyz.' % (TINIT,NEXMD,DIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                ## CALCULATE DIHEDRAL ANGLE ALONG A SINGLE TRAJECTORY ##
                    SDIHEDRAL = np.zeros(NCOORDS - 1)
                    for NCOORD in np.arange(NCOORDS - 1):
                        COORDS = DATA[ARRAY[NCOORD]+1:ARRAY[NCOORD+1]-1:1]
                        VEC0 = np.float_(COORDS[LINES[0]].split()[1:])
                        VEC1 = np.float_(COORDS[LINES[1]].split()[1:])
                        VEC2 = np.float_(COORDS[LINES[2]].split()[1:])
                        VEC3 = np.float_(COORDS[LINES[3]].split()[1:])
                        A = np.subtract(VEC1, VEC0)
                        B = np.subtract(VEC2, VEC1)
                        C = np.subtract(VEC3, VEC2)
                        N1 = np.cross(A, B)
                        N2 = np.cross(B, C)
                        UN1 = N1/np.linalg.norm(N1)
                        UN2 = N2/np.linalg.norm(N2)
                        UB = B/np.linalg.norm(B)
                        M1 = np.cross(UN1, UB)
                        X = np.dot(UN1, UN2)
                        Y = np.dot(M1, UN2)
                        SDIHEDRAL[NCOORD] = np.abs(np.degrees(-math.atan2(Y,X)))
                        EDIHEDRAL[NCOORD,CTRAJ] = SDIHEDRAL[NCOORD]
                    FDIHEDRAL += SDIHEDRAL
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    CTRAJ += 1
                    if TSTEPS == TSMAX:
                        ETRAJ += 1
                else:
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    ERRFLAG = 1
                TTRAJ += 1
        if CTRAJ == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(TSMAX)),TCOLL)
        else:
            EDIHEDRAL = np.delete(EDIHEDRAL, np.arange(CTRAJ, TTRAJ), axis = 1)
            EDIHEDRAL = np.std(EDIHEDRAL, axis = 1)
            FDIHEDRAL = FDIHEDRAL/CTRAJ
            print 'Total Trajectories:', '%04d' % (TTRAJ)
            print 'Completed Trajectories:', '%04d' % (CTRAJ)
            print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            print >> OUTPUT, 'Total Trajectories: ', '%04d' % (TTRAJ)
            print >> OUTPUT, 'Completed Trajectories: ', '%04d' % (CTRAJ)
            print >> OUTPUT, 'Excellent Trajectories: ', '%04d' % (ETRAJ)
            for NCOORD in np.arange(NCOORDS - 1):
                print >> OUTPUT, '%0*.2f' % (len(str((TSMAX))) + 2,DT*ODATA*CDATA*(NCOORD + 1)), '%08.3f' % (FDIHEDRAL[NCOORD]), '%07.3f' % (EDIHEDRAL[NCOORD])
        if ERRFLAG == 1:
            print 'One or more trajectories did not finish within %0*.2f femtoseconds, check dihedral.err.' % (len(str(TSMAX)),TCOLL)
        else:
            os.remove('%s/dihedral.err' % (CWD))
