#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function generates an optical spectrum from single-point     |
| calculations.                                                     |
|                                                                   |
| One of two files will be generated in the current working         |
| directory if an optical spectrum is requested, 'ceo_<type>.out'   |
| or 'ceo.err'.  Here, <type> is 'gauss' or 'lorentz'. These refer  |
| to Gaussian or Lorentzian lineshapes, respectively.  In           |
| 'ceo.out', the spectrum is generated with one of these            |
| lineshapes.  Both the lineshape and the width are defined by the  |
| user.  The outputted spectrum is a sum of all spectra determined  |
| from initial geometries, divided the number of geometries.  The   |
| 'collectceo.sh' script, located in the getexcited_package is      |
| called to grep excitation energies and oscillator strengths from  |
| single-point calculations.  In 'ceo.out', first column is energy  |
| in eV, followed by the spectrum associated to each state,         |
| followed by the total spectrum.  Spectra should be interpreted as |
| relative absorbance in arbitrary units.  In 'ceo.err',            |
| directories of single-point calculations that are incomplete are  |
| listed. The 'ceo.out' file cannot be generated unless all         |
| single-point calculations are complete.                           |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import subprocess
import shlex
import glob

CWD = os.getcwd()

def OPTSPEC(PATHTOCEO):

    print 'Generating optical spectrum.'

## DIRECTORY NAMES ##
    SPDIR = raw_input('Single-point calculations directory: ')
    if not os.path.exists(SPDIR):
        print 'Path %s does not exist.' % (SPDIR)
        sys.exit()
    
## CHECK EXCITATION ENERGIES ##
    print 'Checking energies and oscillator strengths. Please wait ...'
    if not os.path.exists('%s/getexcited_package/collectceo.sh' % (PATHTOCEO)):
        print 'The script, collectceo.sh, must be in the getexcited_package.'
        sys.exit()
    if not os.path.exists('%s/header' % (SPDIR)):
        print 'Path %s/header does not exist.' % (SPDIR)
        sys.exit()
    HEADER = open('%s/header'% (SPDIR),'r')
    for LINE in HEADER:
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
            break
    NEXMDS = glob.glob('%s/NEXMD*/' % (SPDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (SPDIR)
        sys.exit()
    ERROR = open('%s/ceo.err' % (CWD),'w')
    CEOFLAG = 0
    for NEXMD in NEXMDS:
        if not os.path.exists('%s/%s/dirlist1' % (CWD,NEXMD)):
            print 'Path %s/%sdirlist1 does not exist.' % (CWD,NEXMD)
            sys.exit()
        DIRLIST1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (CWD,NEXMD)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        for DIR in DIRLIST1:
            if not os.path.exists('%s/%s/%04d' % (CWD,NEXMD,DIR)):
                print >> ERROR, '%s%04d' % (NEXMD,DIR), 'does not exist'
                CEOFLAG = 1
                continue
            os.chdir('%s/%s/%04d' % (CWD,NEXMD,DIR))
            if not os.path.exists('%s/%s/%04d/md.out' % (CWD,NEXMD,DIR)):
                print >> ERROR, '%s%04d/md.out' % (NEXMD,DIR), 'does not exist'
                CEOFLAG = 1
                continue
            subprocess.call(shlex.split('sh %s/getexcited_package/collectceo.sh %d' % (PATHTOCEO,NSTATES+1)))
            if not os.path.exists('%s/%s/%04d/ceo.out' % (CWD,NEXMD,DIR)):
                print >> ERROR, '%s/%04d/ceo.out' % (NEXMD,DIR), 'does not exist'
                CEOFLAG = 1
                continue
            with open('%s/%s/%04d/ceo.out' % (CWD,NEXMD,DIR),'r') as DATA:
                if len(DATA.readlines()) != NSTATES:
                    print >> ERROR, '%s%04d/ceo.out' % (NEXMD,DIR), 'is incomplete'
                    CEOFLAG = 1
    if CEOFLAG == 1:
        print 'One or more single-point calculations did not finish, check ceo.err.'
        sys.exit()
    else:
        os.remove('%s/ceo.err' % (CWD))

## GENERATE OPTICAL SPECTRUM ##
    STYPE = input('Spectral lineshape? Answer GAUSSIAN [0] or LORENTZIAN [1]: ')
    if STYPE not in [0,1]:
        print 'Answer must be 0 or 1.'
        sys.exit()
    if STYPE == 0:
        SPECB = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
    else:
        SPECB = input('Spectral broadening (i.e. Lorentzian FWHM) in eV [e.g. 0.36]: ')
    if isinstance(SPECB, int) == False and isinstance(SPECB, float) == False:
        print 'Spectral broadening must be integer or float.'
        sys.exit()
    if SPECB < 0:
        print 'Spectral broadening must be integer or float greater than zero.'
        sys.exit()
    NPOINTS = 100000
    DPOINTS = np.zeros((NSTATES+2,NPOINTS))
    DPOINTS[0] = np.linspace(0.0,10.0,NPOINTS)
    TRAJ = 0
    if STYPE == 0:
        for NEXMD in NEXMDS:
            DIRLIST1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (CWD,NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                DATA = np.genfromtxt('%s/%s/%04d/ceo.out' % (CWD,NEXMD,DIR))
                for STATE, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    DPOINTS[STATE] = DPOINTS[STATE] + OSCSTREN*np.exp(-(ENERGY-DPOINTS[0])**(2.0)/(2.0*SPECB**(2.0)))/np.sqrt(2.0*np.pi*SPECB**(2.0))
                print '%s%04d' % (NEXMD,DIR)
                TRAJ += 1
    else:
        for NEXMD in NEXMDS:
            DIRLIST1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (CWD,NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                DATA = np.genfromtxt('%s/%s/%04d/ceo.out' % (CWD,NEXMD,DIR))
                for STATE, ENERGY, OSX, OSY, OSZ, OSCSTREN in DATA:
                    DPOINTS[STATE] = DPOINTS[STATE] + OSCSTREN/(1.0+((ENERGY-DPOINTS[0])**2.0)/(SPECB/2.0)**(2.0))/(SPECB*np.pi/2.0)
                print '%s%04d' % (NEXMD,DIR)
                TRAJ += 1
    DPOINTS[0] = DPOINTS[0]*TRAJ
    DPOINTS[-1] = np.sum(DPOINTS[1:-1:1], 0)
    np.savetxt('%s/ceo_%s.out' % (CWD,'gauss' if STYPE == 0 else 'lorentz'), np.transpose(DPOINTS/TRAJ), fmt='%10.5e')
