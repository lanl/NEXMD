Bugs
1.) ground state dipole prints during geometry optimization when printdipole=1, this should only print before and after geometry optimization
2.) Possible error: profiler logs calc_rhotz instead of packing/unpacking subroutines possibly because of compilation error from 'Entry'.

Optimization
1.) Use 2-electron integrals from ground state for excited state
	Currently this is qm2_fock in ES and qm2_fock_d in GS. Redundant and expensive calculation because these could just be stored from the GS and reused. Also, this adds d-orbitals to ES calculations.
2.) Compare speed of diagonalization subroutines from GS with alternatives from LAPACK or others
3.) Optimize packing/unpacking subroutines and the number of calls to these subroutines
4.) Combine modules (qm2_davidson,qmmm_struct,etc.) which are pointed to by Communism module for clarity and possible memory management/optimization
	This requires a revision of the entire code and room should be left for integrated ambertools using qmmm modules.

Features
1.) State Specific Solvent Z-Vector equation
2.) State Specific Solvent Gradients
3.) QM/MM
4.) Add davidson alternatives from recent Tretiak paper
5.) Write documentation
