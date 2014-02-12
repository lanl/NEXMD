#define STATEINF_FLD_C 5
#define TITR_RES_C 50
#define TITR_STATES_C 200
#define ATOM_CHRG_C 1000
#define MAX_H_COUNT 4

type :: const_ph_info
   sequence
   integer :: num_states, first_atom, num_atoms, first_state, first_charge
end type const_ph_info

type (const_ph_info), parameter :: NULL_CPH_INFO = const_ph_info(0,0,0,0,0)
