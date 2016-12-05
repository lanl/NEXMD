QMMM: SCF Convergence Information
QMMM: (*) = Pseudo Diagonalisation
QMMM: Cycle      Energy          dE             dP              |FP-PF|
QMMM:     1     -30940.72      -30940.72      0.1797693+309  -1.000000    
QMMM: KJ/mol    -129456.0      -129456.0           Infinity
QMMM:     2     -30940.72     -0.5719432E-04  0.7096641E-04  -1.000000    
QMMM: KJ/mol    -129456.0     -0.2393011E-03  0.2969235E-03
QMMM:     3     -30940.72     -0.3422720E-05  0.1876907E-04  -1.000000    
QMMM: KJ/mol    -129456.0     -0.1432066E-04  0.7852980E-04
QMMM:     4 *   -30940.72     -0.4033536E-06  0.5015599E-05  -1.000000    
QMMM: KJ/mol    -129456.0     -0.1687632E-05  0.2098527E-04
QMMM:     5 *   -30940.72      0.1205553E-06  0.5677023E-06  -1.000000    
QMMM: KJ/mol    -129456.0      0.5044035E-06  0.2375266E-05
QMMM:     6 *   -30940.72     -0.8367351E-10  0.1632498E-06  -1.000000    
QMMM: KJ/mol    -129456.0     -0.3500900E-09  0.6830372E-06
QMMM:     7     -30940.72     -0.2910383E-10  0.6850644E-07  -1.000000    
QMMM: KJ/mol    -129456.0     -0.1217704E-09  0.2866310E-06
QMMM: SCF Converged to 0.1000E-07 in:     7 Cycles 
  Restoring quirality in HF orbitals
 
 Start Equilibrium Vertical Excitation Solvent Calculation
 
 SCF Step,  Excitation Energy,    DeltaE_sol, abs(error), error,  COSMO SCF Tole
 rance 
 --------------------------------------------------
 Davidson parameters
 mdflag=                     2
 irflag=                     0
 M4=                   225
 Mx=                     1
 Mj=                     0
 nd=                   500
 nd1=                     3
 j1=                     1
 istore=                     1
 istore_M=                     1
 
 Davidson batch                      1
 So far found                      0  states
 out of requested                      1  states
 This batch will seek                     1  vectors
  Entering davidson0
 Restart Davidson with guesses from the previous step
 Restart loop =                      0
 Currently have                      1  states
 With                      2  initial guesses
                     2                     1  9.212178804568923E-015
 COUNT=                     1 Exp=                     2
 nd1,nd1_old,j0,j1                     2                     0
                     0                     1
 info  -4624372248013176832
    3.82    
 eigenvalues and residual norm
    1    1    0    3.816858870   3.816858877 0.548E-07 0.220E-07   Converged!
 All vectors found after loop                     0 , Expansion 
                     2
 @@@@ Davidson subroutine Found vectors                     1
  i, e0(i), ferr(i), ftol0
  1 +++    3.816858877298256     0.77E-07 0.10E-06
 -------------------------------------------------
 
 Overlaps=   1.01916430722759        1.01916430722759     
 Frequencies (eV) and Oscillator strengths (unitless)
        Omega            fx              fy              fz          ftotal
   1     3.81685887729826        0.121207752238523E-08    0.741527977343448E-08    0.552121995712992E-12    0.862790941781542E-08
 
 Frequencies (eV) and Transition Dipole Moments (AU)
        Omega            fx              fy              fz          ftotal
   1     3.81685887729826        0.113850040565024E-03    0.281599423665135E-03    0.242988385554859E-05    0.922659714807453E-07
 
 Total energy of the ground state (eV,AU)
                     0  -850.381316945535       -31.2509255569009     
 Total energies of excited states (eV,AU)
                     1  -846.562583892832       -31.1105897569795     
 Davidson parameters
 mdflag=                     0
 irflag=                     0
 M4=                   225
 Mx=                     1
 Mj=                     0
 nd=                   500
 nd1=                     3
 j1=                     1
 istore=                     1
 istore_M=                     1
 
 Davidson batch                      1
 So far found                      0  states
 out of requested                      1  states
 This batch will seek                     1  vectors
  Entering davidson0
 Restart Davidson with guesses from the previous step
 Restart loop =                      0
 Currently have                      1  states
 With                      2  initial guesses
 COUNT=                     1 Exp=                     2
 nd1,nd1_old,j0,j1                     2                     0
                     0                     1
 info  -4624372248013176832
    3.82    
 eigenvalues and residual norm
    1    1    0    3.816858870   3.816858877 0.548E-07 0.220E-07   Converged!
 All vectors found after loop                     0 , Expansion 
                     2
 @@@@ Davidson subroutine Found vectors                     1
  i, e0(i), ferr(i), ftol0
  1 +++    3.816858877317500     0.77E-07 0.10E-06
 -------------------------------------------------
 
 Frequencies (eV) and Oscillator strengths (unitless)
        Omega            fx              fy              fz          ftotal
   1     3.81685887731750        0.121207751384096E-08    0.741527970506745E-08    0.552121994253071E-12    0.862790934090267E-08
 
 Frequencies (eV) and Transition Dipole Moments (AU)
        Omega            fx              fy              fz          ftotal
   1     3.81685887731750        0.113850040163457E-03    0.281599422366287E-03    0.242988385232992E-05    0.922659706577830E-07
 
 Total energy of the ground state (eV,AU)
                     0  -850.381316945535       -31.2509255569009     
 Total energies of excited states (eV,AU)
                     1  -846.562583892832       -31.1105897569795     
 Overlaps=   1.01916538712993        1.01916538712993     
  1    3.816858877317500             0.192E-10         0.192E-10        -0.192E-10         0.100E-04        
 
 Final Results of Equilibrium Vertical Excitation Solvent Calculation 
  ES nonlinear term energy=   -0.6459629677E-08 eV
 -------------------------------------------------
 Frequencies (eV) and Oscillator strengths (unitless)
        Omega            fx              fy              fz          ftotal
   1     3.81685887731750        0.121207751384096E-08    0.741527970506745E-08    0.552121994253071E-12    0.862790934090267E-08
 
 Frequencies (eV) and Transition Dipole Moments (AU)
        Omega            fx              fy              fz          ftotal
   1     3.81685887731750        0.113850040163457E-03    0.281599422366287E-03    0.242988385232992E-05    0.922659706577830E-07
 
 Total energy of the ground state (eV,AU)
                     0  -850.381316945535       -31.2509255569009     
 Total energies of excited states (eV,AU)
                     1  -846.564458068218       -31.1106586316294     
 Ground State Molecular Dipole Moment (A.U.)
                         dx              dy              dz          ftotal
                    -0.1801926E-03  0.1147294E-02 -0.5127450E-05  0.1161369E-02
 
 Frequencies (eV) and Total Molecular Dipole Moments (Debye)
        Omega            dx              dy              dz          ftotal
   1   3.816859     -0.4765415E-03  0.2942320E-02  0.2390300E-04  0.2980757E-02
 
 Frequencies (eV) and Total Molecular Dipole Moments (AU)
        Omega            dx              dy              dz          ftotal
   1   3.816859     -0.1874859E-03  0.1157598E-02  0.9404165E-05  0.1172720E-02
 
 Frequencies (eV) Unrelaxed Difference Dipole Moments (AU)
        Omega            dx              dy              dz          ftotal
   1   3.816859     -0.1052468E-03  0.2374845E-03  0.1608941E-04  0.2602588E-03
 
 Frequencies (eV) Relaxed Difference Dipole Moments (AU)
        Omega            dx              dy              dz          ftotal
   1   3.816859     -0.7293285E-05  0.1030408E-04  0.1453161E-04  0.1924926E-04
QMMM: Occupied MO Energies (eV):
      -39.19165319      -31.49449372      -31.49416093      -23.19496285      -23.19477296
      -17.94608424      -16.28071301      -15.48496214      -14.25397939      -14.25376219
      -13.48404852      -11.98039197      -11.98023190       -9.78474386       -9.78438305
QMMM: Virtual MO Energies (eV):
        0.38567901        0.38603439        2.80069134        4.05741708        4.05760894
        4.06620981        4.30286436        4.53514231        4.53532860        5.10090558
        5.10114345        5.54323169        5.54351454        5.60403663        6.03981608
QMMM:
QMMM: SCF Energy =       20.98462080 KCal/mol,        87.79965345 KJ/mol
QMMM: SCF Energy = Heat of formation
QMMM:
QMMM:        Electronic energy =     -1341.69037013 eV (   -30940.72162562 KCal/mol)
QMMM: Dielectric(COSMO) energy =        -0.07189945 eV (       -1.65807331 KCal/mol)
QMMM: QM core - QM core energy =       491.30905319 eV (    11330.07807553 KCal/mol)
QMMM: QM core - MM atom energy =         0.00000000 eV (        0.00000000 KCal/mol)
QMMM: Total core - core energy =       491.30905319 eV (    11330.07807553 KCal/mol)
QMMM:             Total energy =      -850.38131695 eV (   -19610.64355008 KCal/mol)
  Restoring quirality in cmdqt
QMMM: Forces on QM atoms ground state calculation (eV/A)
QMMM: state=  1  0
QMMM: Atm      1:    -0.75187917647493    0.12585183011913   -0.00034961557326
QMMM: Atm      2:    -0.26601826599842    0.71461180561451    0.00018851293693
QMMM: Atm      3:     0.48451664725155    0.58758840092086    0.00017663619060
QMMM: Atm      4:     0.75170861271743   -0.12606315993657   -0.00004953377380
QMMM: Atm      5:     0.26654900124373   -0.71369765529192   -0.00011026149860
QMMM: Atm      6:    -0.48406489596264   -0.58812632792064    0.00000925298092
QMMM: Atm      7:     0.53057386016895   -0.08855649193092   -0.00000679027828
QMMM: Atm      8:     0.18796987729201   -0.50427700266600   -0.00000352482137
QMMM: Atm      9:    -0.34233149704544   -0.41473182080101    0.00004128726932
QMMM: Atm     10:    -0.53070036914164    0.08864453365880    0.00004484506754
QMMM: Atm     11:    -0.18845670263035    0.50426420927877    0.00006502685347
QMMM: Atm     12:     0.34213290857975    0.41449167895499   -0.00000583535346
QMMM: Forces on QM atoms excited state calculation (eV/A)
QMMM: state=  1  0
QMMM: Atm      1:     1.05194857541998   -0.17620800459108   -0.00004723017429
QMMM: Atm      2:     0.37329520593377   -0.99920095399967    0.00002959300212
QMMM: Atm      3:    -0.67846764705180   -0.82312974886978    0.00005493847762
QMMM: Atm      4:    -1.05176051064686    0.17617283619343    0.00002530744441
QMMM: Atm      5:    -0.37325837152349    0.99927267431346   -0.00000168644572
QMMM: Atm      6:     0.67824588394300    0.82309305703942   -0.00002828885180
QMMM: Atm      7:    -0.00389232955257    0.00065357449017    0.00001114094125
QMMM: Atm      8:    -0.00139535639744    0.00371673823027   -0.00001025815952
QMMM: Atm      9:     0.00253987228047    0.00307866949664   -0.00003282642847
QMMM: Atm     10:     0.00385830624365   -0.00064414069973   -0.00001937974981
QMMM: Atm     11:     0.00140586216222   -0.00375438311211   -0.00000248018007
QMMM: Atm     12:    -0.00251949081094   -0.00305031849104    0.00002117012427
QMMM: Excited State Forces on QM atoms from SCF calculation (KJ/mol)
QMMM: Atm      1:   101.49959783302855  -17.00182120957259   -0.00455710840622
QMMM: Atm      2:    36.01821815305740  -96.41012626957945    0.00285534662442
QMMM: Atm      3:   -65.46345983784019  -79.42150446026196    0.00530086119622
QMMM: Atm      4:  -101.48145198513811   16.99842790851114    0.00244184505800
QMMM: Atm      5:   -36.01466410306184   96.41704636356225   -0.00016272046561
QMMM: Atm      6:    65.44206253108662   79.41796416740766   -0.00272951278065
QMMM: Atm      7:    -0.37556007342043    0.06306158823418    0.00107495849353
QMMM: Atm      8:    -0.13463406527965    0.35861775417363   -0.00098978133505
QMMM: Atm      9:     0.24506522565730    0.29705227334428   -0.00316733095693
QMMM: Atm     10:     0.37227725879174   -0.06215134798279   -0.00186989826075
QMMM: Atm     11:     0.13564773735960   -0.36225000431987   -0.00023930568971
QMMM: Atm     12:    -0.24309867424098   -0.29431676351644    0.00204264652275

  QMMM: QM Region Cartesian Coordinates (*=link atom) 
  QMMM: QM_NO. MM_NO.  ATOM         X         Y         Z
  QMMM:     1        1       1      -7.9895    0.6830    0.0000
  QMMM:     2        2       2      -7.0994    1.7625   -0.0000
  QMMM:     3        3       3      -5.7195    1.5314   -0.0000
  QMMM:     4        4       4      -5.2296    0.2208    0.0000
  QMMM:     5        5       5      -6.1197   -0.8587    0.0000
  QMMM:     6        6       6      -7.4996   -0.6276   -0.0000
  QMMM:     7        7       7      -9.0598    0.8622   -0.0000
  QMMM:     8        8       8      -7.4793    2.7791    0.0000
  QMMM:     9        9       9      -5.0291    2.3688    0.0000
  QMMM:    10       10      10      -4.1592    0.0416   -0.0000
  QMMM:    11       11      11      -5.7397   -1.8753   -0.0000
  QMMM:    12       12      12      -8.1900   -1.4649   -0.0000
 
 |^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|
 | MD normal termination at   Tue Nov  8 10:01:26 MST 2016   |
 | MD total CPU time	         0.30000     seconds |
 |      0 days  0 hours  0 minutes  0 seconds     |
 
  SQM (ground state) took overall [s]:
       0.110E-01
  Davidson (excited states) took overall [s]:
       0.106    
  deriv (adiabatic forces) took overall [s]:
       0.115    
  nacT (NA derivatives) took overall [s]:
        0.00    
 |________________________________________________|
 
