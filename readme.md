Bugs
1.) Possible error: profiler logs calc_rhotz instead of packing/unpacking subroutines possibly because of compilation error from 'Entry'.
2.) Excited state energies do not print during iterative solvent calculations instead is junk

Optimization
1.) Use 2-electron integrals from ground state for excited state
	Currently this is qm2_fock in ES and qm2_fock_d in GS. Redundant and expensive calculation because these could just be stored from the GS and reused. Also, this adds d-orbitals to ES calculations. Potential operator in linearized Liouville equation does not use subroutines for d-orbitals, but d-orbitals are allowed in gradients.  
2.) Compare speed of diagonalization subroutines from GS with alternatives from LAPACK or others
3.) Optimize packing/unpacking subroutines and the number of calls to these subroutines
4.) Combine modules (qm2_davidson,qmmm_struct,etc.) which are pointed to by Communism module for clarity and possible memory management/optimization
	This requires a revision of the entire code and room should be left for integrated ambertools using qmmm modules.
5.) Comparison of old and new with unoptimized BLAS and B-O dynamics gives new as faster, nonB-O gives old as faster -- NOT SURE HOW TO USE OPTIMIZED BLAS IN THE NEW CODE
6.) For NACVs specifically, there is a huge speed-up when going to optimized BLAS in the old code

Features
1.) QM/MM
2.) Write documentation
