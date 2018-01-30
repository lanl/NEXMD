#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares input files for single-point calculations  |
| before the push pulse.  This function is used when simulating     |
| pump-push-probe spectroscopy.                                     |
|                                                                   |
| A general header called 'header' must be in the single-point      |
| directory (e.g. singlepoint_push) and must have all inputs set    |
| except for:                                                       |
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
import os
import sys
import shutil
import glob

def SPCALC_PUSH():

    print 'Preparing input files for single-point calculations before push pulse.'

## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory after pump pulse: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()
    OUTDIR = raw_input('Output directory [e.g. singlepoint_push]: ')
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

## LENGTH OF NEXMD BEFORE PUSH PULSE ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header' % (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    for LINE in HEADER:
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1])
        if 'out_data_steps' in LINE:
            ODATS = np.int(LINE.split()[0][len('out_data_steps='):-1])
        if 'out_coords_steps' in LINE:
            OCORS = np.int(LINE.split()[0][len('out_coords_steps='):-1])
    print 'NEXMD trajectories before push pulse were set to run for %.2f fs.\nThe last restart.out file was generated at %.2f fs.\nTherefore, the pump-push delay time is %.2f fs.' % (TSMAX*DT, TSMAX*DT - np.round(TSMAX*DT % DT*ODATS*OCORS, 2), TSMAX*DT - np.round(TSMAX*DT % DT*ODATS*OCORS, 2))
    CONTQ = input('Continue with single-point calculations before push pulse? Answer YES [1] or NO [0]: ')
    if CONTQ not in [1,0]:
    	print 'Answer must be 1 or 0.'
        sys.exit()
    if CONTQ == 0:
    	sys.exit()

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

## PREPARE SINGLE-POINT INPUT FILES BEFORE PUSH PULSE ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()
    ERROR = open('%s/spcalc_push.err' % (CWD),'w')
    PSHFLAG = 0
    for NEXMD in NEXMDS:
        SPNEXMD = os.path.basename(os.path.normpath('%s' % (NEXMD)))
        os.makedirs('%s/%s' % (OUTDIR,SPNEXMD))
        DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        DIRLIST = open('%s/%s/dirlist' % (OUTDIR,SPNEXMD),'w')
        for DIR in DIRLIST1:
            if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'does not exist'
                PSHFLAG = 1
                continue
            DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
            DATA = DATA.readlines()
            TSTEPS = len(DATA)- 2
            if TSTEPS == TSMAX:
                if not os.path.exists('%s/%04d/restart.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/restart.out' % (NEXMD,DIR), 'does not exist'
                    PSHFLAG = 1
                    continue
                DATA = open('%s/%04d/restart.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                INDEX = 0
                ARRAY = np.array([])
                for LINE in DATA:
                    if '$COORD' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$ENDCOORD' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$VELOC' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$ENDVELOC' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    INDEX += 1
                ARRAY = np.int_(ARRAY)
                if len(ARRAY) != 4:
                    print >> ERROR, '%s%04d/restart.out' % (NEXMD,DIR), 'is incomplete'
                    PSHFLAG = 1
                    continue
                COORDS = DATA[ARRAY[0]:ARRAY[1]+1:1]
                VELOCS = DATA[ARRAY[2]:ARRAY[3]+1:1]
                INPUT = open('%s/%04d/input.ceon' % (NEXMD,DIR),'w')
                for LINE in HEADER:
                    if 'nucl_coord_veloc' in LINE:
                        for LINE in COORDS:
                            INPUT.write(LINE)
                        INPUT.write('\n')
                        for LINE in VELOCS:
                            INPUT.write(LINE)
                    else:
                        INPUT.write(LINE)
                print >> DIRLIST, '%04d' % (DIR)
                print '%s/%s/%04d' % (OUTDIR,SPNEXMD,DIR)
            else:
                print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'is incomplete'
                PSHFLAG = 1
        DIRLIST.close()
        shutil.copyfile('%s/%s/dirlist' % (OUTDIR,SPNEXMD), '%s/%s/dirlist1' % (OUTDIR,SPNEXMD))
    if PSHFLAG == 1:
        print 'One or more single-point calculations cannot be prepared, check spcalc_push.err.'
    if os.stat('%s/spcalc_push.err' % (CWD)).st_size == 0:
        os.remove('%s/spcalc_push.err' % (CWD))

