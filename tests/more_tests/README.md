# NEXMD Tests

This folder contains the test input files and corresponding outputs for NEXMD.

- [Outputs of Concern](#outputs-of-concern)
  - [Energies](#energies)
  - [Dipole Moments](#dipole-moments)
  - [MO Energies](#mo-energies)
  - [Coordinates and Velocities](#coordinates-and-velocities)
- [Tests](#tests)
  - [HOW-TO](#how-to)
  - [Geometry Optimization](#geometry-optimization)
  - [Single Point](#single-point)
  - [Born-Oppenheimer Dynamics](#born-oppenheimer-dynamics)
  - [Non-Born-Oppenheimer Dynamics](#non-born-oppenheimer-dynamics)

## Outputs of Concern

- ### Energies

  The energies are printed out in output file `energy-ev.out`, containing the
  kinetic, potential, and total energies for dynamics simulations.

  ```text
  ##     time(fs)           T                  T-T0               U                  U-U0               E                  E-E0
       0.0000000000       0.0000000000       0.0000000000   -1011.1054312585       0.0000000000   -1011.1054312585       0.0000000000
       0.1000000000       0.0004073341       0.0004073341   -1011.1058287368      -0.0003974784   -1011.1054214028       0.0000098557
  ```

  They are also printed to standard output for non-dynamics simulations.

  ```text
  Total energy of the ground state (eV,AU)
            0  -850.30588667453048       -31.248153546583932 
  Total energies of excited states (eV,AU)
            1  -846.48778618470760       -31.107840992911186
  ```

- ### Dipole Moments
  
  The dipole moments are printed in the standard output like the following
  example.

  ```text
  The following dipole moments and excitation energies:
  
  Ground State Molecular Dipole Moment (A.U.)
                           dx              dy              dz          ftotal
                      -0.4815844E-05  0.2294108E-06  0.6140443E-14  0.4821305E-05
   Frequencies (eV) and Total Molecular Dipole Moments (AU)
          Omega            dx              dy              dz          ftotal
     1   3.818106     -0.5477786E-05  0.4489792E-06  0.6351349E-14  0.5496155E-05
  
   Frequencies (eV) Unrelaxed Difference Dipole Moments (AU)
          Omega            dx              dy              dz          ftotal
     1   3.818106     -0.2277329E-05  0.5989293E-06  0.2530236E-15  0.2354770E-05
  
   Frequencies (eV) Relaxed Difference Dipole Moments (AU)
          Omega            dx              dy              dz          ftotal
     1   3.818106     -0.6619418E-06  0.2195684E-06  0.2109053E-15  0.6974075E-06
  ```

  Additionally, transition dipole moments are also printed in a separate file
  called `tdipole.out`.

- ### MO Energies
  
  MO energies are printed in the standard output when the QM/MM verbosity is
  set to 5.

  ```text
  QMMM: Occupied MO Energies (eV):
        -40.43827978      -38.07712635      -35.99393829      -31.99824377      -31.40750968
        -29.05213325      -26.57380282      -23.22100905      -22.73455229      -19.71959530
        -18.37467763      -17.09878667      -16.47470737      -15.68633201      -15.19327751
        -14.57920084      -14.55916978      -14.42180229      -13.72857311      -13.44095597
        -13.13856412      -12.93912014      -12.92217623      -12.89745824      -11.82797833
        -10.85364751      -10.04187937       -8.42406695
  QMMM: Virtual MO Energies (eV):
         -0.04495286        0.33624900        1.84764612        1.97425778        2.08384772
          3.03830386        3.29528755        3.41937862        3.49413172        3.69265996
          3.77592779        3.82572108        4.07215740        4.12328027        4.12754261
          4.19720897        4.34858962        4.43540956        4.76418508        5.00351116
          5.11055190        5.41269531        5.55346819        5.91220581        6.24615063
          6.95468106
  ```

- ### Coordinates and Velocities

  The coordinates are printed in the standard xyz format in file `coords.xyz`.
  For AIMC, the file is named `coords_xxxx.xyz`, where xxxx is the clone index
  like 0001 and 0002, etc.

  ```text
               8
  FINAL HEAT OF FORMATION =      -918.0963448452  time =       0.0000000000
    N   -1.2135490000   -0.5974060000    0.0680630000
    C    0.0000010000    0.1026030000   -0.0006850000
    N    1.2133780000   -0.5978260000   -0.0679130000
    O    0.0002810000    1.3583540000    0.0001900000
    H   -1.2494380000   -1.5212790000   -0.2888710000
    H   -2.0156140000   -0.0362140000   -0.1068760000
    H    1.2485800000   -1.5215190000    0.2895600000
    H    2.0154140000   -0.0368140000    0.1077320000
  ```

  The velocities are printed to `velocity.out` in the unit of Å/fs. Similarly,
  for AIMC, the file is named `velocity_xxxx.out`.

  ```text
  FINAL HEAT OF FORMATION =      -917.1136389812  time =    3463.9999999979
  $VELOC
       9.3412357125     0.0059061964    -3.1297071337
     -11.0244530268     2.6955204955     2.4710651192
      -2.6865336667     3.2497630698    -3.8813984985
       2.1260653419    -6.6755217877     0.2298677958
       3.2769484192     8.9516993088    12.4387036795
     -10.6344153234    33.9182031838     2.8130787320
      23.3431001833   -13.1554552393    55.5928589199
     -10.7991093614    -0.8881600461    -5.9977946989
   $ENDVELOC
   ```

  For geometry optimization, they are also printed at the end of the standard
  output.

## Tests

The details of each test is listed in this section with a hyperlink to the
corresponding folder.

- ### HOW-TO

  Run the input files with your own executable. The output files should be
  identical to the examples provided in the repo.  

- ### Geometry Optimization

  - [Test 1](geomopt/test1): ground state geometry optimization
  - [Test 2](geomopt/test2): first excited state geometry
    optimization
  - [Test 3](geomopt/test3): third excited state geometry
    optimization

- ### Single Point

  - [Test 4](singlepoint/test4): single point ground state and excited state
    SCF properties
  - [Test 12](singlepoint/test4): single point ground state and excited state
    scf properties with state specific solvent

- ### Born-Oppenheimer Dynamics

  - [Test 5](bodynamics/test5): ground state BOMD in NVT ensemble
  - [Test 6](bodynamics/test6): ground state BOMD in NVE ensemble
  - [Test 7](bodynamics/test7): excited state BO Newtonian dynamics (RPA)
  - [Test 8](bodynamics/test8): excited state BO Newtonian dynamics (CIS)
  - [Test 9](bodynamics/test9): excited state BO Langevin dynamics (RPA)
  - [Test 10](bodynamics/test10): excited state BO Newtonian dynamics (RPA) with
    linear response solvent
  - [Test 11](bodynamics/test11): excited state BO Newtonian dynamics (CIS) with
    vertical excitation solvent

- ### Non-Born-Oppenheimer Dynamics

  - [Test 13](nonbodynamics/test13): nonadiabatic molecular dynamics (NAMD)
    using CIS, and Langevin dynamics with decoherence. The initial excited
    state is $S_3$ and $S_5$ for benzene and nitromethane, respectively.
  - [Test 14](nonbodynamics/test14): same as Test 13, but with different random
    seed.
  - [Test 15](nonbodynamics/test15): same as Test 14, but with different
    initial states. Now it is $S_4$ and $S_6$ for benzene and nitromethane,
    respectively.
  - [Test 16](nonbodynamics/test16): NAMD using RPA, and Langevin dynamics with
    decoherence. The initial excited state is $S_3$ and $S_5$ for benzene and
    nitromethane, respectively.
  - [Test 17](nonbodynamics/test17): NAMD using CIS, and Newtonian dynamics
    with decoherence. The initial excited state is $S_3$ and $S_5$ for benzene
    and nitromethane, respectively.
  - [Test 18](nonbodynamics/test18): NAMD using CIS, and Langevin dynamics with
    decoherence only after a successful hop. The initial excited state is $S_3$
    and $S_5$ for benzene and nitromethane, respectively.
  - [Test 19](nonbodynamics/test19): NAMD using CIS, and Langevin dynamics with
    decoherence in linear-response solvent model. The initial excited state is
    $S_4$ and $S_5$ for benzene and nitromethane, respectively.
  - [Test 20](nonbodynamics/test20): NAMD using CIS, and Langevin dynamics with
    decoherence in linear-response solvent model with the Onsager potential.
    The initial excited state is $S_4$ and $S_5$ for benzene and nitromethane,
    respectively.
  - [Test 21](nonbodynamics/test22): NAMD using CIS, and Langevin dynamics with
    decoherence. The trivial crossing routine is tested. The initial excited
    state is $S_5$ and $S_8$ for benzene and nitromethane, respectively.
  - [Test 22](nonbodynamics/test22): excited state non-adiabatic dynamics with
    [ab-initio multiple cloning (AIMC)](nonbodynamics/test22/AIMC) vs
    [Ehrenfest (EHR)](nonbodynamics/test22/EHR)
  - [Test 23](nonbodynamics/test23): NAMD with certain bond distances being
    frozen
  - [Test 24](nonbodynamics/test24): NAMD with certain normal modes being
    frozen
  - [Test 25](nonbodynamics/test25): NAMD with states reduction
