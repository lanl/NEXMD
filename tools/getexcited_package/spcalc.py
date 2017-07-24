#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares input files for single-point calculations. |
|                                                                   |
| A general header called 'header' must be in the single-point      |
| directory (e.g. singlepoint) and must have all inputs set except  |
| for:                                                              |
|                                                                   |
| 1) Initial nuclear coordinates and velocities (nucl_coord_veloc)  |
|                                                                   |
| In parentheses shown above, is the flag in 'header' that labels   |
| this input.  This function finds this label and fills it in       |
| accordingly.                                                      |
|                                                                   |
| NOTE: All NEXMD folders, inside the single-point directory, will  |
| be deleted if this function is completely executed!               |
|___________________________________________________________________|

'''

import numpy as np
import random
import os
import sys
import shutil
import glob

def SPCALC():

    print 'Preparing input files for single-point calculations.'

## DIRECTORY NAMES ##
    GSDIR = raw_input('Ground-state dynamics directory: ')
    if not os.path.exists(GSDIR):
        print 'Path %s does not exist.' % (GSDIR)
        sys.exit()
    OUTDIR = raw_input('Output directory [e.g. singlepoint]: ')
    if not os.path.exists(OUTDIR):
        print 'Path %s does not exist.' % (OUTDIR)
        sys.exit()

## DELETE PREVIOUS NEXMD FOLDERS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (OUTDIR))
    NEXMDS.sort()
    if len(NEXMDS) != 0:
        CONTQ = input('** WARNING ** All NEXMD folders inside %s will be deleted!\nDo you want to continue? Answer YES [1] or NO [0]: ' % (OUTDIR))
        if CONTQ not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if CONTQ == 0:
            sys.exit()
    for NEXMD in NEXMDS:
        print 'Deleting', '%s' % (NEXMD)
        shutil.rmtree(NEXMD)

## CHECK RUNNING SINGLE-POINT ##
    if not os.path.exists('%s/header' % (OUTDIR)):
        print 'Path %s/header does not exist.' % (OUTDIR)
        sys.exit()
    HEADER = open('%s/header' % (OUTDIR),'r')
    HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1])
            break
    if TSMAX != 0:
        print 'User must change n_class_steps in %s/header to 0 for single-point calculations.' % (OUTDIR)
        sys.exit()

## FIND GEOMETRIES ##
    if not os.path.exists('%s/coords.xyz' % (GSDIR)):
        print 'Path %s/coords.xyz does not exist.' % (GSDIR)
        sys.exit()
    DATAC = open('%s/coords.xyz' % (GSDIR),'r')
    DATAC = DATAC.readlines()
    LENC = len(DATAC)
    if not os.path.exists('%s/velocity.out' % (GSDIR)):
        print 'Path %s/velocity.out does not exist.' % (GSDIR)
        sys.exit()
    DATAV = open('%s/velocity.out' % (GSDIR),'r')
    DATAV = DATAV.readlines()
    LENV = len(DATAV)
    NCOORDS = 0
    INDEX = 0
    ARRAYC = np.array([])
    for LINE in DATAC:
        if 'time' in LINE:
            if NCOORDS == 0:
                TINIT = np.float(LINE.split()[-1])
            else:
                TIME = np.float(LINE.split()[-1])
            NCOORDS += 1
            ARRAYC = np.append(ARRAYC,INDEX)
        INDEX += 1
    ARRAYC = np.append(ARRAYC,LENC + 1)
    ARRAYC = np.int_(ARRAYC)
    if NCOORDS == 0:
        print 'No coordinates were found.'
        sys.exit()
    if NCOORDS == 1:
        print 'Only initial coordinates, at %.2f fs, were found.' % (TINIT)
        sys.exit()
    if NCOORDS > 1:
        INDEX = 0
        ARRAYV = np.array([0])
        for LINE in DATAV:
            if 'time' in LINE:
                ARRAYV = np.append(ARRAYV,INDEX)
            INDEX += 1
        ARRAYV = np.append(ARRAYV,LENV)
        ARRAYV = np.int_(ARRAYV)
    TINC = TIME/(NCOORDS - 1)
    print 'A total of %d coordinates files, ranging from %.2f to %.2f fs in increments of %.2f fs, were found.' % (NCOORDS,TINIT,TIME,TINC)

## CHOOSE GEOMETRIES ##
    COORDS = input('Enter requested range of the ground-state sampling by coordinate files and the number of single-point calculations.\nInput an array of the form [start, end, number]: ')
    if not isinstance(COORDS,list):
        print 'Input must be an array of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.'
        sys.exit()
    if len(COORDS) != 3:
        print 'Input must be an array with three elements of the form [start, end, number].\nFor example, [1, 1000, 500] requests 500 coordinate files sampled from 1 to 1000.'
        sys.exit()
    INDEX = 0
    for i in COORDS:
        if isinstance(i, int) == False:
            print 'Element number %d of input array must be integer.\nUser inputted [%s, %s, %s], which is not allowed.' % (INDEX + 1,COORDS[0],COORDS[1],COORDS[2])
            sys.exit()
        if INDEX in [0,1]:
            if i not in np.arange(NCOORDS):
                print 'Element number %d of input array must be integer between 0 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (INDEX + 1,NCOORDS - 1,COORDS[0],COORDS[1],COORDS[2])
                sys.exit()
        else:
            if i not in np.arange(1,NCOORDS + 1) and i not in -np.arange(1,NCOORDS + 1):
                print 'Element number %d of input array must be integer between 1 and %d.\nUser inputted [%d, %d, %d], which is not allowed.' % (INDEX + 1,NCOORDS,COORDS[0],COORDS[1],COORDS[2])
                sys.exit()
        INDEX += 1
    if COORDS[0] > COORDS[1]:
        print 'Second element of input array must be greater than first element.\nUser inputted [%d, %d, %d], which is not allowed.' % (COORDS[0],COORDS[1],COORDS[2])
        sys.exit()
    if (COORDS[1] - COORDS[0]) + 1 < np.abs(COORDS[2]):
        print 'Number of coordinate files requested must be less than or equal to the number of coordinate files in the sample.\nUser inputted [%d, %d, %d], which is not allowed.' % (COORDS[0],COORDS[1],COORDS[2])
        sys.exit()
    if COORDS[2] < 0:
        NTRAJ = np.abs(COORDS[2])
        COORDSQ = input('You have requested %d randomly-selected coordinate files in the range %d to %d.\nContinue? Answer YES [1] or NO [0]: ' % (NTRAJ,COORDS[0],COORDS[1]))
    else:
        INTERVAL = np.int(np.ceil((COORDS[1] - COORDS[0] + 1)/np.float(COORDS[2])))
        if INTERVAL*COORDS[2] > (COORDS[1] - COORDS[0] + 1):
            INTERVAL = INTERVAL - 1
            NTRAJ = len(np.arange(COORDS[0],COORDS[1] + 1,INTERVAL))
        else:
            NTRAJ = COORDS[2]
        COORDSQ = input('You have requested %d evenly-spaced coordinate files in the range %d to %d.\nContinue? Answer YES [1] or NO [0]: ' % (NTRAJ,COORDS[0],COORDS[1]))
    if COORDSQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if COORDSQ == 0:
        sys.exit()

## SPLIT GEOMETRIES ##
    SPLIT = input('Number of single-point calculations per NEXMD folder [e.g. 100]: ')
    if isinstance(SPLIT, int) == False:
        print 'Number of single-point calculations per NEXMD folder must be integer.'
        sys.exit()
    if SPLIT < 0:
        print 'Number of single-point calculations per NEXMD folder must be integer greater than zero.'
        sys.exit()
    DIRSPLIT = SPLIT*np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT) + 1)
    DIRSPLIT[-1] = COORDS[1] + 1
    if COORDS[2] < 0:
        DIRSPLIT = np.split(np.sort(random.sample(np.arange(COORDS[0],COORDS[1]+1),NTRAJ)),DIRSPLIT)
    else:
        DIRSPLIT = np.split(np.arange(COORDS[0],COORDS[1]+1,INTERVAL),DIRSPLIT)

## EXTRACT ATOMIC NUMBERS ##
    if not os.path.exists('%s/restart.out' % (GSDIR)):
        print 'Path %s/restart.out does not exist.' % (GSDIR)
        sys.exit()
    ANUM = open('%s/restart.out' % (GSDIR),'r')
    ANUM = ANUM.readlines()
    TOP = None
    BOTTOM = None
    INDEX = 0
    for LINE in ANUM:
        if '$COORD' in LINE:
            TOP = INDEX
        if '$ENDCOORD' in LINE:
            BOTTOM = INDEX
            break
        INDEX += 1
    if isinstance(TOP, int) == True and isinstance(BOTTOM, int) == True:
        ANUM = [ LINE.split()[0] for LINE in ANUM[TOP+1:BOTTOM:1] ]
    else:
        print 'There is a problem with %s/restart.out.' % (GSDIR)
        sys.exit()

## PREPARE SINGLE-POINT INPUT FILES ##
    INDEX = 0
    for NEXMD in np.arange(1,np.ceil(np.float(NTRAJ)/SPLIT)+1):
        os.makedirs('%s/NEXMD%d' % (OUTDIR,NEXMD))
        DIRLIST = open('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD),'w')
        for DIR in DIRSPLIT[INDEX]:
            COORDS = DATAC[ARRAYC[DIR]+1:ARRAYC[DIR+1]-1:1]
            VELOCS = DATAV[ARRAYV[DIR]+2:ARRAYV[DIR+1]-1:1]
            os.makedirs('%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR))
            INPUT = open('%s/NEXMD%d/%04d/input.ceon' % (OUTDIR,NEXMD,DIR),'w')
            for LINE in HEADER:
                if 'nucl_coord_veloc' in LINE:
                    INPUT.write('&coord\n')
                    AINDEX = 0
                    for LINE in COORDS:
                        VAL = LINE.split()
                        INPUT.write('{:>6}  {:>12}  {:>12}  {:>12}'.format(ANUM[AINDEX],VAL[1],VAL[2],VAL[3]))
                        INPUT.write('\n')
                        AINDEX += 1
                    INPUT.write('&endcoord\n\n&veloc\n')
                    for LINE in VELOCS:
                        INPUT.write(LINE)
                    INPUT.write('&endveloc\n')
                else:
                    INPUT.write(LINE)
            print >> DIRLIST, '%04d' % (DIR)
            print '%s/NEXMD%d/%04d' % (OUTDIR,NEXMD,DIR)
        DIRLIST.close()
        shutil.copyfile('%s/NEXMD%d/dirlist' % (OUTDIR,NEXMD), '%s/NEXMD%d/dirlist1' % (OUTDIR,NEXMD))
        INDEX += 1
