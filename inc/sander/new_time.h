!+ Data declarations for new_time routines.

! All but BC_TIME_PAR and Tpar_P (used as a pointer to the whole
! time_ptr common) are strictly internal to new_time.f.

integer maxtime
parameter (maxtime = 200)
integer t_maxtree
parameter (t_maxtree = 10)
integer bc_time_par  ! size in integers of common time_ptr
parameter (bc_time_par = 9 * maxtime)

integer tpar_p,tchild_p,tsib_p,t_level,t_state, &
      t_added,tnext_p,tprev_p,t_print
common/TIME_ptr/tpar_p(maxtime),tchild_p(maxtime), &
      tsib_p(maxtime),t_level(maxtime),t_state(maxtime), &
      t_added(maxtime),tnext_p(maxtime),tprev_p(maxtime), &
      t_print(maxtime)

_REAL_ t_accum,tch_acc,t_curr,t_other
common/TIME_val/t_accum(maxtime),tch_acc(maxtime), &
      t_curr(maxtime),t_other(maxtime)

character(len=20) t_string(maxtime)
common/TIME_lab/t_string

integer t_numtree,t_first,t_last
common/timeglob/t_numtree,t_first(t_maxtree), &
      t_last(t_maxtree)
