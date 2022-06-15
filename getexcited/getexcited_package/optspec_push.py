#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function generates an optical spectrum from single-point     |
| calculations after the pump-push delay time.                      |
|                                                                   |
| One of two files will be generated in the current working         |
| directory if an optical spectrum is requested, 'muab_<type>.out'  |
| or 'muab.err'.  In 'muab_<type>.out', the spectrum is given       |
| either a Gaussian or Lorentzian lineshape, where the type and     |
| width of the lineshape are defined by the user. The outputted     |
| spectrum is a sum of all spectra determined from initial          |
| geometries, divided the number of geometries.  In                 |
| 'muab_<type>.out', first column is energy in eV and second column |
| is relative absorbance in arbitrary units.  In 'muab.err',        |
| incomplete single-point-calculation files that are required for   |
| generating an optical spectrum are listed.  The 'muab_<type>.out' |
| file cannot be generated unless all single-point calculations are |
| complete.                                                         |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import glob

CWD = os.getcwd()

def OPTSPEC(PATHTOCEO):

    print('Generating optical spectrum from single-point calculations after pump-push delay time.')

## DIRECTORY NAMES ##
    SPDIR = raw_input('Single-point calculations directory: ')
    if not os.path.exists(SPDIR):
        print('Path %s does not exist.' % (SPDIR))
        sys.exit()
    
## CHECK EXCITATION ENERGIES AND EXCITED-TO-EXCITED OSCILLATOR STRENGTHS ##
    print('Checking energies and oscillator strengths. Please wait ...')
    if not os.path.exists('%s/header' % (SPDIR)):
        print('Path %s/header does not exist.' % (SPDIR))
        sys.exit()
    HEADER = open('%s/header'% (SPDIR),'r')
    for LINE in HEADER:
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
            break
    NEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print('There are no NEXMD folders in %s.' % (SPDIR))
        sys.exit()
    ERROR = open('%s/muab.err' % (CWD),'w')
    MUABFLAG = 0
    for NEXMD in NEXMDS:
        if not os.path.exists('%s/dirlist1' % (NEXMD)):
            print('Path %sdirlist1 does not exist.' % (NEXMD))
            sys.exit()
        DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        for DIR in DIRLIST1:
            if not os.path.exist('%s/%04d/restart.out' % (NEXMD,DIR)):
                ERROR.wirte( '%s%04d/restart.out' % (NEXMD,DIR), 'does not exist')
                MUABFLAG = 1
            if not os.path.exists('%s/%04d/muab.out' % (NEXMD,DIR)):
                ERROR.wirte( '%s%04d/muab.out' % (NEXMD,DIR), 'does not exist')
                MUABFLAG = 1
            else:
                with open('%s/%04d/muab.out' % (NEXMD,DIR)) as DATA:
                    if len(DATA.readlines()) != np.sum(np.arange(NSTATES + 2)):
                        ERROR.wirte( '%s%04d/muab.out' % (NEXMD,DIR), 'is incomplete')
                        MUABFLAG = 1
    if MUABFLAG == 1:
        print('One or more single-point calculations did not finish, check muab.err.')
        sys.exit()
    else:
        os.remove('%s/muab.err' % (CWD))

## INDICES FOR EXCITED-TO-EXCITED OSCILLATOR STRENGTHS ##
    INDICES = np.cumsum(np.insert(np.arange(1,NSTATES + 2)[::-1],[0],0,axis=0))

## GENERATE OPTICAL SPECTRUM ##
    STYPE = input('Spectral lineshape? Answer GAUSSIAN [0] or LORENTZIAN [1]: ')
    if STYPE not in [0,1]:
        print('Answer must be 0 or 1.')
        sys.exit()
    if STYPE == 0:
        SPECB = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
    else:
        SPECB = input('Spectral broadening (i.e. Lorentzian FWHM) in eV [e.g. 0.36]: ')
    if isinstance(SPECB, int) == False and isinstance(SPECB, float) == False:
        print('Spectral broadening must be integer or float.')
        sys.exit()
    if SPECB < 0:
        print('Spectral broadening must be integer or float greater than zero.')
        sys.exit()
    NPOINTS = 100000
    DPOINTS = np.zeros((2,NPOINTS))
    DPOINTS[0] = np.linspace(0.0,10.0,NPOINTS)
    TRAJ = 0
    if STYPE == 0:
        for NEXMD in NEXMDS:
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                DATA = open('%s/%04d/restart.out' % (NEXMD,DIR),'r')
                for LINE in DATA:
                    if 'State' in LINE:
                        STATE = np.int(LINE.split()[-1])
                        break
                DATA = np.genfromtxt('%s/%04d/muab.out' % (NEXMD,DIR))[INDICES[STATE-1]:INDICES[STATE]-1:1]
                for LOW, HIGH, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    DPOINTS[1] = DPOINTS[1] + OSCSTREN*np.exp(-(ENERGY-DPOINTS[0])**(2.0)/(2.0*SPECB**(2.0)))/np.sqrt(2.0*np.pi*SPECB**(2.0))
                print('%s%04d' % (NEXMD,DIR))
                TRAJ += 1
    else:
        for NEXMD in NEXMDS:
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                DATA = open('%s/%04d/restart.out' % (NEXMD,DIR),'r')
                for LINE in DATA:
                    if 'State' in LINE:
                        STATE = np.int(LINE.split()[-1])
                        break
                DATA = np.genfromtxt('%s/%04d/muab.out' % (NEXMD,DIR))[INDICES[STATE-1]:INDICES[STATE]-1:1]
                for LOW, HIGH, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    DPOINTS[1] = DPOINTS[1] + OSCSTREN/(1.0+((ENERGY-DPOINTS[0])**2.0)/(SPECB/2.0)**(2.0))/(SPECB*np.pi/2.0)
                print('%s%04d' % (NEXMD,DIR))
                TRAJ += 1
    DPOINTS[0] = DPOINTS[0]*TRAJ
    np.savetxt('%s/muab.out' % (CWD),np.transpose(DPOINTS/TRAJ),fmt='%10.5e')
    
