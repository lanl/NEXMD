#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function prepares restart input files for NEXMD.             |
|                                                                   |
| If restart input files are requested, the function searches for   |
| trajectories that did not complete up to the user-defined         |
| number of classical steps and prepares input files with the name  |
| 'input.ceon' in their respective directories.  The input file is  |
| prepared with coordinates, velocities, quantum coefficients, and  |
| last-residing surface taken from the 'restart.out' file.  These   |
| are data from the last-generated time-step.  Within each 'NEXMD#' |
| folder, with directories that contain restart input files, a new  |
| 'dirlist' will be generated with a list of restart directories.   |
| An error file, 'restart.err', lists the directories that either   |
| (1) do not contain 'coeff-n.out' files, from which populations    |
| are collected, or (2) have incomplete 'restart.out' files.  The   |
| former generally means the trajectory did not start when NEXMD    |
| was first attempted.  The 'restart.err' file is not generated if  |
| there are no such trajectories.  During every iteration of        |
| requesting restart input files, a file containing the random      |
| seeds is generated with the name 'rseedslist#'.  Part of this     |
| function deletes the data at all time-steps after the time-step   |
| of the 'restart.out' file.  The purpose of this is to have        |
| continuous set of data along the trajectory with no repeated      |
| information.  An error file called 'delextra.out' is generated if |
| one or more output files do not exist.                            |
|___________________________________________________________________|

'''


import numpy as np
import os
import sys
import re
import glob
import filecmp
import subprocess
import shlex

CWD = os.getcwd()

def EXTRACT(FILE):
    NUM = re.findall('\d+$', FILE)
    return (np.int(NUM[0]) if NUM else 0, FILE)

def RESTART(PATHTODEL):

    print 'Preparing restart input files for NEXMD.'

## DIRECTORY NAMES ##
    NEXMDIR = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()

## CHOOSE CLASSICAL TIME-STEPS, GET NUMBER OF QUANTUM STEPS AND VERBOSITY ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header' % (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    NUM = 0
    TLINE = len(HEADER)
    VERB = None
    for LINE in HEADER:
        if 'bo_dynamics_flag' in LINE:
            BOFLAG = np.int(LINE.split()[0][len('bo_dynamics_flag='):-1])
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1])
        if 'n_quant_steps' in LINE:
            NQSTEP = np.int(LINE.split()[0][len('n_quant_steps='):-1])
            if NQSTEP == 0:
                NQSTEP = 1
        if '&moldyn' in LINE:
            TLINE = NUM
        if 'verbosity' in LINE and NUM > TLINE and VERB is None:
            VERB = np.int(LINE.split()[0][len('verbosity='):-1])
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
            if ODATA == 0:
                print 'No data has been printed to files because out_data_steps = 0 in header.'
                sys.exit()
        NUM += 1
    print 'Currently, trajectories are set to run for %d classical steps with a time-step of %.2f fs.\nThis is a total of %.2f fs.' % (TSMAX,DT,TSMAX*DT)
    TSMAXQ = input('Keep this trajectory length? Answer YES [1] or NO [0]: ')
    if TSMAXQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if TSMAXQ == 0:
        NTSMAX = input('Enter new number of classical time-steps: ')
        if isinstance(NTSMAX, int) == False:
            print 'Answer must be integer.'
            sys.exit()
        if NTSMAX <= TSMAX:
            print 'Answer must be greater than or equal to the previous number of classical steps used, which was %d.' % (TSMAX)
            sys.exit()
        NHEADER = open('%s/nheader' % (NEXMDIR),'w')
        for LINE in HEADER:
            if 'n_class_steps' in LINE:
                NHEADER.write('   n_class_steps=%d, ! Number of classical steps [1]\n' % (NTSMAX))
            else:
                NHEADER.write(LINE)
        NHEADER.close()
        TSMAX = NTSMAX
        os.rename('%s/nheader' % (NEXMDIR), '%s/header' % (NEXMDIR))
        
## CHOOSE RANDOM SEEDS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()
    NTRAJ = 0
    for NEXMD in NEXMDS:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print 'Path %sdirlist1 does not exist.' % (NEXMD)
            sys.exit()
        DATA = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DATA,int) == True:
            DATA = np.array([DATA])
        NTRAJ += len(DATA)
    RANDQ = input('New random seeds? Answer YES [1] or NO [0]: ')
    if RANDQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    DATA = glob.glob('%s/rseedslist*' % (NEXMDIR))
    MAXDIR = np.int(EXTRACT(max(DATA,key = EXTRACT))[0])
    if RANDQ == 0:
        RSEEDS = raw_input('Path to random-seeds list (** must be different from past random seeds **): ')
        if not os.path.exists(RSEEDS):
            print 'Path %s does not exist.' % (RSEEDS)
            sys.exit()
        RSEEDS = np.int_(np.genfromtxt('%s' % (RSEEDS)))
        if isinstance(RSEEDS,int) == True:
            RSEEDS = np.array([RSEEDS])
        LEN = len(RSEEDS)
        if LEN < NTRAJ:
            print 'Length of random-seeds list must be equal to or greater than the number of trajectories.\nUser inputted a random-seeds list of length %d, while the number of trajectories is %d.' % (LEN,NTRAJ)
            sys.exit()
        for RSEED in RSEEDS:
            if RSEED < 0:
                print 'A negative random seed was detected, %d.\nWithin the getexcited_package, a negative seed is assigned to a trajectory that could not be prepared due to some problem.' % (RSEED)
                sys.exit()

## PREPARE NEXMD RESTART INPUT FILES ##
    HEADER = open('%s/header' % (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    RSEEDSLIST = open('%s/rseedslist%d' % (NEXMDIR,MAXDIR + 1),'w')
    ERROR = open('%s/restart.err' % (CWD),'w')
    RTIMES = np.array([])
    RSTFLAG = 0
    TRAJ = 0
    for NEXMD in NEXMDS:
        DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        DIRLIST = open('%s/dirlist' % (NEXMD),'w')
        for DIR in DIRLIST1:
            if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'does not exist'
                print >> RSEEDSLIST, '%d' % (-123456789)
                RSTFLAG = 1
                TRAJ += 1
                continue
            DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
            DATA = DATA.readlines()
            TSTEPS = len(DATA) - 2
            if TSTEPS != TSMAX:
                if not os.path.exists('%s/%04d/restart.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s/%04d/restart.out' % (NEXMD,DIR), 'does not exist'
                    print >> RSEEDSLIST, '%d' % (-123456789)
                    RSTFLAG = 1
                    TRAJ += 1
                    continue
                DATA = open('%s/%04d/restart.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                INDEX = 0
                ARRAY = np.array([])
                for LINE in DATA:
                    if 'time' in LINE:
                        TIME = np.float(LINE.split()[-1])
                    if 'State' in LINE:
                        STATE = np.int(LINE.split()[-1])
                    if 'Seed' in LINE:
                        RSEED = np.int(LINE.split()[-1])
                    if '$COORD' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$ENDCOORD' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$VELOC' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$ENDVELOC' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$COEFF' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    if '$ENDCOEFF' in LINE:
                        ARRAY = np.append(ARRAY,INDEX)
                    INDEX += 1
                ARRAY = np.int_(ARRAY)
                if len(ARRAY) != 6:
                    print >> ERROR, '%s%04d/restart.out' % (NEXMD,DIR), 'is incomplete'
                    print >> RSEEDSLIST, '%d' % (-123456789)
                    RSTFLAG = 1
                    TRAJ += 1
                    continue
                COORDS = DATA[ARRAY[0]:ARRAY[1]+1:1]
                VELOCS = DATA[ARRAY[2]:ARRAY[3]+1:1]
                COEFFS = DATA[ARRAY[4]:ARRAY[5]+1:1]
                DATA = glob.glob('%s/view*' % (NEXMDIR))
                DATA = [ x[:-10] for x in DATA ]
                if len(DATA) != 0:
                    MAX = np.int(EXTRACT(max(DATA,key = EXTRACT))[0])
                else:
                    MAX = 0
                INPUT = open('%s/%04d/input.ceon' % (NEXMD,DIR),'w')
                for LINE in HEADER:
                    if 'rnd_seed' in LINE:
                        INPUT.write('   rnd_seed=%d, ! Seed for the random number generator\n' % (RSEED if RANDQ == 1 else RSEEDS[TRAJ]))
                    else:
                        if 'exc_state_init_flag' in LINE:
                            INPUT.write('   exc_state_init=%d, ! Initial excited state (0 - ground state) [0]\n' % (STATE))
                        else:
                            if 'time_init' in LINE:
                                INPUT.write('   time_init=%.1f, ! Initial time, fs [0.00]\n' % (TIME))
                            else:
                                if 'n_class_steps' in LINE:
                                    INPUT.write('   n_class_steps=%d, ! Number of classical steps [1]\n' % (np.int(TSMAX-TIME/DT)))
                                else:
                                    if 'out_count_init' in LINE:
                                        INPUT.write('   out_count_init=%d, ! Initial count for output files [0]\n' % (MAX))
                                    else:
                                        if 'nucl_coord_veloc' in LINE:
                                            for LINE in COORDS:
                                                INPUT.write(LINE)
                                            INPUT.write('\n')
                                            for LINE in VELOCS:
                                                INPUT.write(LINE)
                                        else:
                                            if 'quant_amp_phase' in LINE:
                                                for LINE in COEFFS:
                                                    INPUT.write(LINE)
                                            else:
                                                INPUT.write(LINE)
                RTIMES = np.append(RTIMES,TIME)
                print >> DIRLIST, '%04d' % (DIR)
                print '%s%04d' % (NEXMD,DIR)
                print >> RSEEDSLIST, '%d' % (RSEED if RANDQ == 1 else RSEEDS[TRAJ])
            else:
                print >> RSEEDSLIST, '%d' % (-123456789)
            TRAJ += 1
        DIRLIST.close()
    if filecmp.cmp('%s/rseedslist%d' % (NEXMDIR,MAXDIR + 1), '%s/rseedslist%d' % (NEXMDIR,MAXDIR + 1)):
        os.remove('%s/rseedslist%d' % (NEXMDIR,MAXDIR + 1))
    if RSTFLAG == 1:
        print 'One or more trajectories cannot be restarted, check restart.err.'
    else:
        os.remove('%s/restart.err' % (CWD))

## DELETE EXTRANEOUS DATA IN OUTPUT FILES ##
    if RSTFLAG == 1:
        CONTQ = input('Continue to delete extraneous data in output files? Answer YES [1] or NO [0]: ')
        if CONTQ not in [1,0]:
            print 'Answer must be 1 or 0.'
            sys.exit()
        if CONTQ == 0:
            sys.exit()
    print 'Deleting extraneous data in output files. Please wait ...'
    if BOFLAG == 0:
        FILES = ['coeff-n.out', 'energy-ev.out', 'nact.out', 'order.out', 'pes.out', 'temperature.out', 'transition-densities.out']
        OFILES = ['hops.out', 'hops-trial.out']
    else:
        FILES = ['energy-ev.out', 'temperature.out']
        OFILES = []
    ERROR = open('%s/delextra.err' % (CWD),'w')
    RSTFLAG = 0
    TRAJ = 0
    for NEXMD in NEXMDS:
        DIRLIST = np.int_(np.genfromtxt('%s/%s/dirlist' % (CWD,NEXMD)))
        if isinstance(DIRLIST,int) == True:
            DIRLIST = np.array([DIRLIST])
        for DIR in DIRLIST:
            if not os.path.exists('%s/%s/%04d' % (CWD,NEXMD,DIR)):
                print >> ERROR, '%s/%s%04d' % (CWD,NEXMD,DIR), 'does not exist'
                RSTFLAG = 1
                TRAJ += 1
                continue
            os.chdir('%s/%s/%04d' % (CWD,NEXMD,DIR))
            for INDEX in np.arange(len(FILES)):
                if not os.path.exists('%s/%s/%04d/%s' % (CWD,NEXMD,DIR,FILES[INDEX])):
                    print >> ERROR, '%s/%s%04d/%s' % (CWD,NEXMD,DIR,FILES[INDEX]), 'does not exist'
                    RSTFLAG = 1
                    continue
                ## DERIVATION OF THE FOLLOWING ALGORITHM IS PROVIDED AT THE END OF THIS SCRIPT ##
                DATA = subprocess.check_output(['tail','-1','%s' % (FILES[INDEX])])
                LTIME = np.float(np.fromstring(DATA,dtype=float,sep=' ')[1 if INDEX == 0 else 0])
                if RTIMES[TRAJ] > LTIME:
                    print >> ERROR, 'Last time-step in', '%s/%s%04d/%s' % (CWD,NEXMD,DIR,FILES[INDEX]), 'exceeds time-step in', '%s/%s%04d/%s' % (CWD,NEXMD,DIR,'restart.out')
                    RSTFLAG = 1
                    continue
                if VERB == 3 and INDEX in [2,4]:
                    NCSTEPS = 0
                    while RTIMES[TRAJ] + NCSTEPS*(ODATA*DT) <= LTIME:
                        NCSTEPS += 1
                    LCTIME = RTIMES[TRAJ] + (NCSTEPS - 1)*(ODATA*DT)
                    if LTIME == LCTIME:
                        NLINES = np.ceil((LCTIME - RTIMES[TRAJ])*(ODATA*(NQSTEP - 1) + 1)/(ODATA*DT) + NQSTEP*(LTIME - LCTIME)/DT + 1)
                    else:
                        NCSTEPS = 0
                        while RTIMES[TRAJ] + NCSTEPS*DT < LTIME:
                            NCSTEPS += 1
                        NLINES = np.ceil((LCTIME - RTIMES[TRAJ])*(ODATA*(NQSTEP - 1) + 1)/(ODATA*DT) + NQSTEP*(LTIME - LCTIME)/DT - (NCSTEPS - 1) + 1)
                else:
                    NLINES = np.ceil((LTIME - RTIMES[TRAJ])/(ODATA*DT) + 1)
                if BOFLAG == 0:
                    subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (PATHTODEL,(NLINES - 1 if INDEX in [2,3,5] else NLINES),FILES[INDEX])))
                else:
                    subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (PATHTODEL,(NLINES - 1 if INDEX == 1 else NLINES),FILES[INDEX])))
                os.rename('%s.restart' % (FILES[INDEX]), '%s' % (FILES[INDEX]))
            for INDEX in np.arange(len(OFILES)):
                if not os.path.exists('%s/%s/%04d/%s' % (CWD,NEXMD,DIR,OFILES[INDEX])):
                    print >> ERROR, '%s/%s%04d/%s' % (CWD,NEXMD,DIR,OFILES[INDEX]), 'does not exist'
                    RSTFLAG = 1
                    continue
                LTIME = 1000000
                NLINES = 0
                while LTIME >= RTIMES[TRAJ]:
                    DATA = subprocess.check_output(['tail','%d' % (-(NLINES + 1)),'%s' % (OFILES[INDEX])])
                    LTIME = np.float(np.fromstring(DATA,dtype=float,sep=' ')[0])
                    NLINES += 1
                subprocess.call(shlex.split('sh %s/getexcited_package/cutdata.sh %d %s' % (PATHTODEL, NLINES - 1, OFILES[INDEX])))
                os.rename('%s.restart' % (OFILES[INDEX]), '%s' % (OFILES[INDEX]))
            print '%s%04d' % (NEXMD,DIR)
            TRAJ += 1
    if RSTFLAG == 1:
        print 'One or more trajectories cannot be restarted properly, check delextra.err.'
    else:
        os.remove('%s/delextra.err' % (CWD))

'''
    
    DT      = Time-step
    VERB    = Verbosity of molecular dynamics data
    ODATA   = Data are printed every ODATA classical time-steps
    LTIME   = Last printed time-step (may be classical or quantum)
    LCTIME  = Last printed classical time-step
    RTIME   = Restart time-step (always classical)
    NQSTEPS = Number of quantum steps per classical step
    
    A = (LCTIME - RTIME)/DT                                       = Total # of classical steps computed
    B = (LCTIME - RTIME)/(ODATA x DT)                             = Total # of classical steps printed
    C = (LCTIME - RTIME)/DT x ((ODATA - 1)/ODATA)                 = Total # of classical steps not printed
    D = (NQSTEPS - 1) x (LCTIME - RTIME)/DT x ((ODATA - 1)/ODATA) = Total # of quantum steps from unprinted classical steps
    E = NQSTEPS x (LCTIME - RTIME)/(ODATA x DT)                   = Total # of classical and quantum steps from printed classical steps
    
    If VERB != 3:
        Total # of lines to delete = B
    If VERB == 3 and FILE not in [nact.out, pes.out]:
        Total # of lines to delete = B
    If VERB == 3 and FILE in [nact.out, pes.out]:
        Total # of lines to delete = D + E + Resdiual quantum steps after LCTIME
        If LTIME == LCTIME:
            Residual quantum steps after LCTIME = NQSTEP x (LTIME - LCTIME)/DT
            Total # of lines to delete = D + E + NQSTEP x (LTIME - LCTIME)/DT + 1 ( + 1 to delete the restart time-step )
        If LTIME != LCTIME:
            Residual quantum steps after LCTIME = NQSTEP x (LTIME - LCTIME)/DT
            Number of unprinted classical steps after LCTIME = NCSTEPS - 1
            Total # of lines to delete = D + E + NQSTEP x (LTIME - LCTIME)/DT - (NCSTEPS - 1) + 1 ( + 1 to delete the restart time-step)
'''
