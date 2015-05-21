Notes
1.) The nonadiabatic coupling vector NACR from nacT_analytic_module is calculated using a derivative subroutine, i.e. dcart1_xpm, while in the previous NAESMD version it was calculated using two Fockian constructions and finite difference. This was much faster before, since the gradient calculates the energy on an atom by atom basis. The gradient in the NACR is only necessary along the nuclear trajectory and thus we should go back to the original form of calculating NACR.

Bugs
1.) Possible error: profiler logs calc_rhotz instead of packing/unpacking subroutines possibly because of compilation error from 'Entry'.
2.) Excited state energies do not print during iterative solvent calculations instead is junk
3.) Crash when printing excited state dipoles and geometry optimization of the ground state 

Optimization
1.) Use 2-electron integrals from ground state for excited state
	Currently this is qm2_fock in ES and qm2_fock_d in GS. Redundant and expensive calculation because these could just be stored from the GS and reused. Also, this adds d-orbitals to ES calculations. Potential operator in linearized Liouville equation does not use subroutines for d-orbitals, but d-orbitals are allowed in gradients.  
2.) Compare speed of diagonalization subroutines from GS with alternatives from LAPACK or others
3.) Optimize packing/unpacking subroutines and the number of calls to these subroutines
4.) Combine modules (qm2_davidson,qmmm_struct,etc.) which are pointed to by Communism module for clarity and possible memory management/optimization
	This requires a revision of the entire code and room should be left for integrated ambertools using qmmm modules.
6.) For NACVs specifically, there is a huge speed-up when going to optimized BLAS in the old code

Features
1.) QM/MM
2.) Write documentation
