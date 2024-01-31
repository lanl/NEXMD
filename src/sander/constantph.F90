! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

! This module handles all of the constant pH MD capabilities in Amber.
module constantph

#  include "dynph.h"
#  include "random.h"

! Public variable declaration
    _REAL_, public  :: chrgdat(0:ATOM_CHRG_C - 1)
    integer, public :: cph_igb, cphfirst_sol, relaxations, tot_relax

! Public for pH REM
    _REAL_, public  :: target_ph

! Private variable declaration
    type(const_ph_info), private  :: stateinf(0:TITR_RES_C - 1)
    type(rand_gen_state), private :: cph_gen

    character(len=40), private :: resname(0:TITR_RES_C)
    _REAL_, private           :: statene(0:TITR_STATES_C - 1)
    integer, private           :: hindex(0:TITR_RES_C*(MAX_H_COUNT + 1) - 1), &
        tpair(0:TITR_RES_C*4), resstate(0:TITR_RES_C - 1), &
        iselres(0:2), iselstat(0:2), &
        protcnt(0:TITR_STATES_C - 1)
    integer, private           :: trescnt, nres_titrate
    logical, private           :: cph_success, first

!  Variable descriptions
!
!  cph_gen        : random number generator for constant pH MC
!  stateinf       : information about titratable residues. See dynph.h
!  chrgdat        : partial atomic charges for every state of every titr. res.
!  statene        : relative energy of every state for each titr. res.
!  hindex         : location of each titratable hydrogen; see below
!  tpair          : neighborlist for interacting titr. residues; see below
!  resstate       : list of current protonation states for each titr. res.
!  iselres        : array of selected, titrating residues
!  iselstat       : array of selected, proposed, states for selected titr. res.
!  protcnt        : number of protons for each state of each titr. res.
!  trescnt        : number of titrating residues
!  cphfirst_sol   : atom number of the first solvent residue
!  cph_igb        : igb model to use for MC evals. in explicit solvent pH MD
!  cph_success    : did the proton exchange attempt succeed or not
!  first          : is this the first pass through? (for CPOUT printing)
!  nres_titrate   : how many residues we are going to titrate
!  relaxations    : how many times we have run water relaxations (for timing)
!  tot_relax      : how many relaxations we've done in total (for timing)

! tpair and hindex arrays:
!  These arrays contain index information for each titr. res., and have similar
!  structures. The head of each array is a pointer section, from indices 0 to
!  trescnt.  The residue number (0 to trescnt-1) is an index that shows the
!  location of the FIRST element belonging to that residue (be it neighbors or
!  titrating protons). You can get the *number* of elements belonging to residue
!  "x" by subtracting array(x+1) - array(x). The trescnt-th element of each
!  array is a tail index so that trick also works for the last titr. res.
!  Let me show you an example using hindex: Consider 2 titr. res.  The 0th res
!  has titrating proton numbers 12, 15, 16, and 17. The 1st res has titrating
!  proton numbers 22 and 23. The hindex array would then be arranged as follows:
!  hindex = 3, 7, 9, 12, 15, 16, 17, 22, 23

! Make necessary subroutines public
    public cnstphinit, cnstphupdatepairs, cnstphbeginstep, cnstphendstep, &
        cnstphwriterestart, cnstphwrite, cnstph_explicitmd, cnstphread, &
#ifdef MPI
        cnstph_bcast, &
#endif
        cnstph_zero

! Make internal subroutines private
    private cnstphupdatechrg

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Zero all of the constant pH arrays
    subroutine cnstph_zero()
        use constants, only : ZERO
        implicit none

        ! I zero out all of the arrays at first because when the cprestrt is written,
        ! it simply writes the whole namelist, but unused values would not otherwise
        ! be initialized.

        ! NULL_CPH_INFO is a zero-filled const_ph_info type. See dynph.h

        stateinf(0:TITR_RES_C - 1) = NULL_CPH_INFO
        chrgdat(0:ATOM_CHRG_C - 1) = 0
        hindex(0:TITR_RES_C*(MAX_H_COUNT + 1) - 1) = 0
        tpair(0:TITR_RES_C*4 - 1) = 0
        resstate(0:TITR_RES_C - 1) = 0
        iselres(0:2) = 0
        iselstat(0:2) = 0
        statene(0:TITR_STATES_C - 1) = ZERO
        protcnt(0:TITR_STATES_C - 1) = 0
        resname(0:TITR_RES_C) = ' '
        relaxations = 0
        tot_relax = 0

    end subroutine cnstph_zero

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the CPIN file and initialize arrays
    subroutine cnstphread(charge)
        use constants, only : AMBER_ELECTROSTATIC, ZERO, ONE
        use file_io_dat, only : CNSTPH_UNIT, cpin
        implicit none
#  include "md.h"

        ! Variable Descriptions:
        !
        ! Passed:
        !  charge         : partial charge array
        !
        ! Internal:
        !  itres          : residue loop pointer
        !  res            : residue pointer
        !  istat          : state pointer
        !  iatom          : atom pointer
        !  icumstat       : sum of states
        !  icumchrg       : sum of charge
        !
        ! Common memory:
        !  gbgamma        : (md.h) OBC GB parameter
        !  gbbeta         : (md.h) OBC GB parameter
        !  gbalpha        : (md.h) OBC GB parameter

        _REAL_, intent(inout) :: charge(*)

        integer  :: itres, res, istat, iatom, icumstat, icumchrg

        ! This subroutine reads in and initializes the constant pH data. It
        ! multiplies in the AMBER_ELECTROSTATIC factor and adjusts the initial charge
        ! array for the initial protonation states.

        namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
            trescnt, resname, cphfirst_sol, cph_igb

        icumstat = -1
        icumchrg = 0

        write (6, '(a,a)') 'reading charge increments from file: ', cpin
        call amopen(CNSTPH_UNIT, cpin, 'O', 'F', 'R')
        read (CNSTPH_UNIT, nml=cnstph)
        do iatom = 0, ATOM_CHRG_C - 1
            chrgdat(iatom) = chrgdat(iatom)*AMBER_ELECTROSTATIC
        end do
        close (CNSTPH_UNIT)

        ! Set initial charges
        do itres = 0, trescnt - 1
            do iatom = 0, stateinf(itres)%num_atoms - 1
                charge(iatom + stateinf(itres)%first_atom) &
                    = chrgdat(stateinf(itres)%first_charge + iatom + resstate(itres)* &
                    stateinf(itres)%num_atoms)
            end do
        end do

        ! Set proper GB stuff
        if (cph_igb .eq. 2) then
            gbgamma = 2.90912499999d0
            gbbeta = ZERO
            gbalpha = 0.8d0
        else if (cph_igb .eq. 5) then
            gbgamma = 4.851d0
            gbbeta = 0.8d0
            gbalpha = ONE
        end if

        ! Overflow checking
        if (trescnt > TITR_RES_C) then
            write (6, *) 'Too many titrating residues; alter dynph.h and recompile'
            call mexit(6, 1)
        end if
        do itres = 0, trescnt - 1
            if (stateinf(itres)%first_state > icumstat) then
                icumstat = stateinf(itres)%first_state
                res = itres
            end if
        end do

        icumstat = stateinf(res)%first_state + stateinf(res)%num_states
        icumchrg = stateinf(res)%first_charge + stateinf(res)%num_atoms* &
            stateinf(res)%num_states

        if (icumstat > TITR_STATES_C) then
            write (6, *) 'Too many titrating states; alter dynph, recompile'
            call mexit(6, 1)
        end if
        if (icumchrg > ATOM_CHRG_C) then
            write (6, *) 'Too much charge data; alter dynph.h, recompile'
            call mexit(6, 1)
        end if

        if (icnstph > 1) then
            nres_titrate = trescnt
            write (6, '(a,i4,a,i4,a)') '| Attempting ', nres_titrate, &
                ' MC protonation changes every ', ntcnstph, ' steps.'
        end if

        target_ph = solvph

    end subroutine cnstphread

#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Broadcasts all of the constant pH information to each thread
    subroutine cnstph_bcast(ierr)
        implicit none
#include "parallel.h"
        include 'mpif.h'
        integer, intent(out) :: ierr

        ! Broadcast everything
        call mpi_bcast(stateinf, TITR_RES_C*STATEINF_FLD_C, MPI_INTEGER, 0, &
            commsander, ierr)
        call mpi_bcast(trescnt, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(chrgdat, ATOM_CHRG_C, MPI_DOUBLE_PRECISION, 0, &
            commsander, ierr)
        call mpi_bcast(resstate, TITR_RES_C, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(statene, TITR_STATES_C, MPI_DOUBLE_PRECISION, 0, &
            commsander, ierr)
        call mpi_bcast(protcnt, TITR_STATES_C, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(resname, TITR_RES_C + 1, MPI_CHARACTER, 0, commsander, ierr)
        call mpi_bcast(cphfirst_sol, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(cph_igb, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(nres_titrate, 1, MPI_INTEGER, 0, commsander, ierr)

    end subroutine cnstph_bcast
#endif /* MPI */

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find all titr. protons in each titr. res. for cnstphupdatepairs
    subroutine cnstphinit(x, ig)
        implicit none

#ifdef MPI
        include 'mpif.h'
#  include "parallel.h"
#endif

        ! Variable Descriptions:
        !
        ! Passed:
        !  x              : global position array
        !  ig             : random seed -- used for initializing cph_gen
        !
        ! Local:
        !  res, state     : descriptive pointers in hindex initialization loops
        !  atom, i        : (descriptive) pointers in hindex initialization loops
        !  index_pointer  : index for initializing portions of hindex array
        !  used_atoms     : memory for which titr. protons have been found
        !  is_used        : whether the current proton has been found yet or not

        _REAL_, intent(in)  :: x(*)
        integer, intent(in) :: ig

        integer :: res, state, atom, i, index_pointer
        logical :: is_used
        integer, dimension(MAX_H_COUNT) :: used_atoms
        integer :: seed, ierr

        ! This subroutine will initialize the hindex array according to the format
        ! described above. We must find each titr. proton in each res. exactly once.
        ! We will use the fact that a titr. proton must have a charge of 0 in at
        ! least 1 state to find them. However, it can be 0 in multiple states, so
        ! take care not to double-count. At the end, initialize the tpair array for
        ! the first time. It also initializes the cph_gen random generator

        ! Every thread in parallel MUST have an identical random number stream.
        ! However, setting ig=-1 in the input file forces desynchronization of the
        ! random numbers. For simplicity's sake, we'll just broadcast the master's
        ! "ig" value (but we don't want to clobber it, so copy it here)

        seed = ig
#ifdef MPI
        call mpi_bcast(seed, 1, mpi_integer, 0, commsander, ierr)
#endif
        call amrset_gen(cph_gen, seed)

        ! The first element points past the header section
        index_pointer = trescnt + 1

        do res = 0, trescnt - 1
            hindex(res) = index_pointer
            used_atoms(1:MAX_H_COUNT) = 0

            do atom = 0, stateinf(res)%num_atoms - 1
                do state = 0, stateinf(res)%num_states - 1
                    if (chrgdat(stateinf(res)%first_charge + &
                        state*stateinf(res)%num_atoms + atom) .eq. 0.d0) then

                        is_used = .false.
                        do i = 1, MAX_H_COUNT
                            if (used_atoms(i) .eq. (stateinf(res)%first_atom + atom)) &
                                is_used = .true.
                        end do ! i = 1, MAX_H_COUNT

                        if (.not. is_used) then
                            hindex(index_pointer) = stateinf(res)%first_atom + atom
                            index_pointer = index_pointer + 1

                            do i = 1, MAX_H_COUNT
                                if (used_atoms(i) .eq. 0) then
                                    used_atoms(i) = stateinf(res)%first_atom + atom
                                    exit
                                end if ! used_atoms(i) .eq. 0
                            end do ! i = 1, MAX_H_COUNT

                        end if ! .not. is_used
                    end if ! chrgdat(stateinf(res)%first_charge ...) .eq. 0.d0
                end do ! state = 0, ...
            end do ! atom = 0, ...
        end do ! res = 0, ...

        ! set tail index
        hindex(trescnt) = index_pointer

        ! this is our first pass through
        first = .true.

        call cnstphupdatepairs(x)

    end subroutine cnstphinit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Sets the neighborlist (tpair) for "coupled" titrating residues
    subroutine cnstphupdatepairs(x)
        implicit none

        ! Variable descriptions:
        !
        ! Passed:
        !  x              : global coordinate array
        !
        ! Local:
        !  resi, resj     : residue loop counters
        !  hi, hj         : proton loop counters for residues i, j
        !  tpair_index    : index for tpair array initialization
        !  ti, tj         : pointer for atom i, j in global position arrays
        !  xij, yij, zij  : xj - xi, yj - yi, zj - zi
        !  CUTOFF         : parameter that determines "neighbor" cutoff criteria

        _REAL_, intent(in)   :: x(*)

        integer  :: resi, resj, hi, hj, ti, tj, tpair_index
        _REAL_   :: xij, yij, zij

        _REAL_, parameter :: CUTOFF = 4.0d0

        ! This subroutine sets the tpair array according to the description above.
        ! Two residues are considered "neighbors" if the square of their distance is
        ! less than or equal to CUTOFF defined above

        ! First element points beyond the header
        tpair_index = trescnt + 1

        do resi = 0, trescnt - 1
            tpair(resi) = tpair_index
            do resj = 0, trescnt - 1
                if (resi .eq. resj) cycle

                hloop: do hi = 0, hindex(resi + 1) - hindex(resi) - 1
                    ti = 3*hindex(hi + hindex(resi))

                    do hj = 0, hindex(resj + 1) - hindex(resj) - 1
                        tj = 3*hindex(hj + hindex(resj))
                        xij = x(ti - 2) - x(tj - 2)
                        yij = x(ti - 1) - x(tj - 1)
                        zij = x(ti) - x(tj)
                        if (4.0d0 .gt. xij*xij + yij*yij + zij*zij) then
                            if (tpair_index < TITR_STATES_C*4) then
                                tpair(tpair_index) = resj
                                tpair_index = tpair_index + 1
                                exit hloop
                            else
                                write (6, *) "Constant pH pair list overflow. &
                                &Increase TITR_RES_C in dynph.h; recompile"
                                call mexit(6, 1)
                            end if
                        end if
                    end do
                end do hloop
            end do ! resj
        end do ! resi

        tpair(trescnt) = tpair_index

    end subroutine cnstphupdatepairs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
    subroutine cnstphbeginstep(dcharge)
        implicit none

        ! Variable descriptions:
        !
        ! Passed:
        !  dcharge        : charge array for the proposed state
        !
        ! Internal:
        !  i, j           : loop counters
        !  neighborcount  : how many found neighbors we have
        !  randval        : random value from random number generator

        _REAL_, intent(inout) :: dcharge(1:*)

        integer  :: i, j, neighborcount
        _REAL_   :: randval

        ! This subroutine will randomly select what to titrate. There is a 25% chance
        ! that we will attempt a multi-site move among neighboring residues. This is
        ! done so we can at times observe a proton transfer when that would otherwise
        ! be a very very rare event. The code only currently supports a 2-site move,
        ! but iselres/iselstat and associated loops have been generalized to an
        ! arbitrary number of sites.
        ! First determine if we will do a multisite move. Then randomly select a res
        ! to titrate. If we tried a multisite move, see if our selected site has any
        ! neighbors, then select from among those to titrate also. Update the dcharge
        ! array for the proposed states. If we have no neighbors, don't do multisite.

        call amrand_gen(cph_gen, randval)
        if (randval .lt. 0.25d0) then
            iselres(0) = 2
        else
            iselres(0) = 1
        end if

        ! select a random residue
        call amrand_gen(cph_gen, randval)
        iselres(1) = int((randval*0.9999999d0)*trescnt)

        ! select a different, random state
        call amrand_gen(cph_gen, randval)
        iselstat(1) = int((randval*0.9999999d0)* &
            (stateinf(iselres(1))%num_states - 1))
        if (iselstat(1) .ge. resstate(iselres(1))) &
            iselstat(1) = iselstat(1) + 1
        call cnstphupdatechrg(dcharge, iselres(1), iselstat(1))

        if (iselres(0) .eq. 2) then
            neighborcount = tpair(iselres(1) + 1) - tpair(iselres(1))
            if (neighborcount .gt. 0) then
                call amrand_gen(cph_gen, randval)
                iselres(2) = tpair(tpair(iselres(1)) + int(randval*0.9999999d0* &
                    neighborcount))
                call amrand_gen(cph_gen, randval)
                iselstat(2) = int((randval*0.9999999d0)* &
                    (stateinf(iselres(2))%num_states - 1))
                if (iselstat(2) .ge. resstate(iselres(2))) iselstat(2) = iselstat(2) + 1
                call cnstphupdatechrg(dcharge, iselres(2), iselstat(2))
            else
                iselres(0) = 1
            end if
        end if

    end subroutine cnstphbeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Same as cnstphbeginstep, but for explicit solvent, so we can do multi-MC
    subroutine excnstphbeginstep(dcharge, done_res, num_done)
        implicit none

        ! Variable descriptions:
        !
        ! Passed:
        !  dcharge        : charge array for the proposed state
        !  done_res       : array for which residues have been done
        !  num_done       : how many residues have been finished
        !
        ! Internal:
        !  i, j           : loop counters
        !  neighborcount  : how many found neighbors we have
        !  randval        : random value from random number generator
        !  holder_res     : holder variable so we insert selected residue in order

        _REAL_, intent(inout)  :: dcharge(1:*)
        integer, intent(inout) :: done_res(1:TITR_RES_C)
        integer, intent(in)    :: num_done

        integer  :: i, j, neighborcount, holder_res(1:TITR_RES_C)
        _REAL_   :: randval

        ! This subroutine will randomly select what to titrate, as long as it has.
        ! not yet been titrated. That is the only difference from cnstphbeginstep

        call amrand_gen(cph_gen, randval)
        if (randval .lt. 0.25d0) then
            iselres(0) = 2
        else
            iselres(0) = 1
        end if

        ! select a random residue
        call amrand_gen(cph_gen, randval)
        iselres(1) = int((randval*0.9999999d0)*(trescnt - num_done))

        if (num_done .eq. 0) then
            done_res(1) = iselres(1)
        else
            ! now move the selected residue off of those that have been selected already
            ! THIS ARRAY SHOULD ALREADY BE SORTED (see below)
            do i = 1, num_done
                if (iselres(1) .ge. done_res(i)) iselres(1) = iselres(1) + 1
            end do

            ! If the residue we chose was the larger than our last number, add it to
            ! the end
            if (iselres(1) .gt. done_res(num_done)) &
                done_res(num_done + 1) = iselres(1)
        end if ! (num_done .eq. 0)

        ! now add this to the list of residues that we've titrated already, but
        ! keep this list in order!
        do i = 1, num_done
            if (iselres(1) .lt. done_res(i)) then
                holder_res(i:num_done) = done_res(i:num_done)
                done_res(i) = iselres(1)
                done_res(i + 1:num_done + 1) = holder_res(i:num_done)
                exit
            end if
        end do

        ! select a different, random state
        call amrand_gen(cph_gen, randval)
        iselstat(1) = int((randval*0.9999999d0)* &
            (stateinf(iselres(1))%num_states - 1))
        if (iselstat(1) .ge. resstate(iselres(1))) &
            iselstat(1) = iselstat(1) + 1
        call cnstphupdatechrg(dcharge, iselres(1), iselstat(1))

        if (iselres(0) .eq. 2) then
            neighborcount = tpair(iselres(1) + 1) - tpair(iselres(1))
            if (neighborcount .gt. 0) then
                call amrand_gen(cph_gen, randval)
                iselres(2) = tpair(tpair(iselres(1)) + int(randval*0.9999999d0* &
                    neighborcount))
                call amrand_gen(cph_gen, randval)
                iselstat(2) = int((randval*0.9999999d0)* &
                    (stateinf(iselres(2))%num_states - 1))
                if (iselstat(2) .ge. resstate(iselres(2))) iselstat(2) = iselstat(2) + 1
                call cnstphupdatechrg(dcharge, iselres(2), iselstat(2))
            else
                iselres(0) = 1
            end if
        end if

    end subroutine excnstphbeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
    subroutine cnstphendstep(dcharge, charge, dvdl, temp, solvph)
        use constants, only : KB, LN_TO_LOG
        implicit none
#include "extra.h"

        ! Variable descriptions
        !
        ! Passed:
        !  dcharge        : charge array for the proposed state
        !  charge         : charge array for the current state
        !  dvdl           : change in energy between the two states
        !  temp           : system temperature
        !  solvph         : solvent pH
        !
        ! Internal:
        !  randval        : random value from random number generator
        !  deltae         : Non-electrostatic factors adjusting dvdl
        !  i              : loop counter
        !  statebase      : the first state of the selected res. in state arrays

        _REAL_, intent(in)    :: temp, solvph
        _REAL_, intent(inout) :: dcharge(1:*), charge(1:*)
        _REAL_, intent(inout) :: dvdl

        _REAL_   :: randval, deltae
        integer  :: i, statebase

        ! This subroutine adjusts the delta E passed from electrostatic calculations
        ! for non-EEL factors, and uses the final result to evaluate a MC transition.
        ! If successful, it will call cnstphupdatechrg charge to modify the charge
        ! array to be the same as the dcharge array and update the resstate to be
        ! the newly updated residue state. If it fails, it will reset the
        ! dcharge array to be the same as the charge array. Then it will modify the
        ! cph_success variable based on whether or not we succeeded.

        deltae = 0.d0
        do i = 1, iselres(0)
            statebase = stateinf(iselres(i))%first_state
            !  deltae = E(proposed state) - E(current state)
            deltae = deltae + statene(iselstat(i) + statebase) - &
                statene(resstate(iselres(i)) + statebase)
            !  correct for pH (delta protons * pH * LN_TO_LOG * KB * T)
            deltae = deltae - (protcnt(iselstat(i) + statebase) - &
                protcnt(resstate(iselres(i)) + statebase))*solvph* &
                LN_TO_LOG*KB*temp
        end do
        call amrand_gen(cph_gen, randval)

        dvdl = dvdl - deltae

        if ((dvdl .lt. 0.d0) .or. (randval .le. exp(-dvdl/(KB*temp)))) then
            cph_success = .true.
            do i = 1, iselres(0)
                call cnstphupdatechrg(charge, iselres(i), iselstat(i))
                resstate(iselres(i)) = iselstat(i)
            end do
        else
            do i = 1, iselres(0)
                call cnstphupdatechrg(dcharge, iselres(i), resstate(iselres(i)))
            end do
        end if

    end subroutine cnstphendstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Update charges to reflect new state
    subroutine cnstphupdatechrg(charge, res, inewstate)
        implicit none

        ! Variables
        !
        ! Passed:
        !  charge         : charge array to be updated
        !  res            : the titr. res. to change the charges of
        !  inewstate      : the charge state to change to
        !
        ! Internal:
        !  i              : loop counter

        _REAL_, intent(inout) :: charge(*)
        integer, intent(in)   :: res, inewstate

        integer :: i

        do i = 0, stateinf(res)%num_atoms - 1
            charge(i + stateinf(res)%first_atom) &
                = chrgdat(stateinf(res)%first_charge &
                + inewstate*stateinf(res)%num_atoms + i)
        end do

    end subroutine cnstphupdatechrg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write the cprestrt file
    subroutine cnstphwriterestart(inchrgdat)
        use constants, only : AMBER_ELECTROSTATIC
        use file_io_dat, only : CNSTPH_UNIT, cprestrt, owrite, CPOUT_UNIT
        implicit none

        ! Variable descriptions
        !
        ! Passed:
        !  inchrgdat      : chrgdat array - sent as a copy to multiply EEL factor
        !
        ! Internal:
        !  chrgdat        : chgrdat array adjusted by AMBER_ELECTROSTATIC
        !  iatom          : atom counter
        !  stat           : file status
        !  first          : is this the first pass through?

        _REAL_, intent(in)   :: inchrgdat(0:ATOM_CHRG_C - 1)

        _REAL_            :: chrgdat(0:ATOM_CHRG_C - 1)
        integer           :: iatom
        character(len=7)  :: stat
        logical, save     :: first = .true.

        namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
            trescnt, resname, cphfirst_sol, cph_igb

        do iatom = 0, ATOM_CHRG_C - 1
            chrgdat(iatom) = inchrgdat(iatom)/AMBER_ELECTROSTATIC
        end do

        if (first) then
            if (owrite == 'N') then
                stat = 'NEW'
            else if (owrite == 'O') then
                stat = 'OLD'
            else if (owrite == 'R') then
                stat = 'REPLACE'
            else if (owrite == 'U') then
                stat = 'UNKNOWN'
            end if
            open (unit=CNSTPH_UNIT, file=cprestrt, status=stat, form='FORMATTED', &
                delim='APOSTROPHE')
            first = .false.
        else
            open (unit=CNSTPH_UNIT, file=cprestrt, status='OLD', form='FORMATTED', &
                delim='APOSTROPHE')
        end if

        write (CNSTPH_UNIT, nml=cnstph)
        close (CNSTPH_UNIT)

        ! flush all cpout data
        call amflsh(CPOUT_UNIT)

    end subroutine cnstphwriterestart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write data to the cpout file
    subroutine cnstphwrite(rem)
        use file_io_dat, only : CPOUT_UNIT, ntwx
        implicit none
#  include "md.h"

        ! Variables
        !
        ! Passed:
        !  rem            : remd method used
        !
        ! Internal:
        !  full           : do we write a full record or not
        !  i              : loop counter
        !
        ! Common memory:
        !  irespa         : (md.h) step counter for updating slow-varying terms
        !  nstlim         : (md.h) limit of how many steps we run
        !  solvph         : (md.h) solvent pH
        !  ntcnstph       : (md.h) how often MC protonation jumps are attempted
        !  t              : (md.h) time

        integer, intent(in) :: rem

        logical  :: full
        integer  :: i

        if (ntwx .gt. 0) then
            full = (first .or. irespa .eq. nstlim .or. mod(irespa, ntwx) .eq. 0)
        else
            full = (first .or. irespa .eq. nstlim)
        end if

        if (full) then
            write (CPOUT_UNIT, '(a,f8.5)') 'Solvent pH: ', solvph
            write (CPOUT_UNIT, '(a,i8)') 'Monte Carlo step size: ', ntcnstph
            write (CPOUT_UNIT, '(a,i8)') 'Time step: ', irespa
            write (CPOUT_UNIT, '(a,f10.3)') 'Time: ', t
            if (rem == 4) then
                do i = 0, trescnt - 1
                    write (CPOUT_UNIT, '(a,i4,a,i2,a,f7.3)') &
                        'Residue ', i, ' State: ', resstate(i), ' pH: ', solvph
                end do
            else
                do i = 0, trescnt - 1
                    write (CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ', i, ' State: ', resstate(i)
                end do
            end if
        else
            if (rem == 4) then
                do i = 1, iselres(0)
                    write (CPOUT_UNIT, '(a,i4,a,i2,a,f7.3)') 'Residue ', iselres(i), &
                        ' State: ', resstate(iselres(i)), ' pH: ', solvph
                end do
            else
                do i = 1, iselres(0)
                    write (CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ', iselres(i), ' State: ', &
                        resstate(iselres(i))
                end do
            end if
        end if
        write (CPOUT_UNIT, '()')
        first = .false.

    end subroutine cnstphwrite

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculates the total number of "active" protons
    subroutine total_protonation(tot_prot)

        implicit none

! Passed Variables
        integer, intent(out) :: tot_prot

! Local Variables
        integer              :: i ! counter

        tot_prot = 0

        do i = 0, trescnt - 1
            tot_prot = tot_prot + protcnt(stateinf(i)%first_state + resstate(i))
        end do

    end subroutine total_protonation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Sets up explicit constant pH transition and water relaxation
    subroutine cnstph_explicitmd(xx, ix, ih, ipairs, x, winv, amass, f, v, vold, &
        xr, xc, conp, skip, nsp, tma, erstop, qsetup, &
        do_list_update, rem, mdloop)
        use file_io_dat, only : ntpr, ntwr, ntwx
        use sgld, only : tlangv
        use state
        implicit none

#include "box.h"
#include "extra.h"
#include "md.h"
#include "memory.h"
#include "nmr.h"

!  Passed Variables
!
!  XX    : Global real array
!  X     : Global coordinate array
!  WINV  : Inverted masses
!  AMASS : Atomic masses
!  F     : Force array
!  V     : Velocity array
!  VOLD  : Old velocity array (for velocity verlet integrator)
!  XR    : Coords relative to COM of molecule
!  XC    : Restraint reference coords
!  CONP  : Bond parameter for SHAKE
!  SKIP  : Skip array for shake with logical value for each atom
!  NSP   : Submolecule index array (?)
!  TMA   : Submolecular weight array (?)
!  ERSTOP: Should we stop in error (?)
!  QSETUP: Not quite sure...
!  DO_LIST_SETUP : logical that indicates if we should rebuild pairlist (?)

!  Local Variables
!
!  holder_cut     : temporarily store the &cntrl CUT value while we do GB
!  holder_ntb     : temporarily store the &cntrl NTB value while we do GB
!  holder_natom   : temporarily store the true number of atoms
!  holder_ntwx    : temporarily store the ntwx value
!  last_solute    : last atom of the solute (first atom of solvent - 1)
!  cph_ener       : record for storing energies
!  iselres        : constant pH residue selection (change prot state)
!  iselstat       : constant pH state selection for chosen residues
!  vtemp          : temporary array holder for the velocities
!  holder_ntt     : temporary holder for the thermostat being used
!  holder_tcouple : temporary holder for the temperature coupling constant.
!                   gamma_ln needs to be 0 to turn off langevin dynamics
!  holder_tlangv  : holder for value of tlangv
!  vtemp          : holder array for velocities
!  holder_nscm    : holder variable for nscm (see below)
!  holder_ndfmin  : holder variable for ndfmin (see below)
!  holder_irespa  : holder variable for irespa (see below)
!  natom3         : natom * 3
!  selres_holder  : residues that have been chosen to titrate
!  rem_holder     : holder value for REM
!  mdloop_holder  : holder value for mdloop

!  Other used variables from common memory
!
!  natom       : (memory.h) -- number of atoms in the system
!  igb         : (md.h) -- GB solvent model to use
!  ntb         : (box.h) -- PBC (1 -> V, 2 -> P, 0 -> none)
!  ntp         : (md.h) -- constant pressure control
!  cut         : (box.h) -- value of cutoff
!  relaxing    : (md.h) -- tag to determine if we are relaxing waters or not
!  master      : (extra.h) -- True if this is PID 0, False otherwise
!  ntt         : (md.h) -- thermostat
!  gamma_ln    : (md.h) -- langevin coupling constant
!  tlangv      : (sgld.h) -- logical for if this run is langevin dynamics
!  nscm        : (md.h) -- how often COM terms are removed
!  ndfmin      : (md.h) -- which COM terms are removed
!  irespa      : (md.h) -- counter to determine when to evaluate terms
!
!  NOTE: All pointer variables into IX and XX are defined in memory.h and
!        described in locmem.f, where they're used during the allocation
!        of the global integer, real, and hollerith arrays.

        ! Passed arguments
        _REAL_  :: xx(*), x(*), winv(*), amass(*), f(*), v(*), vold(*), &
            xr(*), xc(*), conp(*), tma(*), factt
        logical :: skip(*), erstop, qsetup, do_list_update
        integer :: ix(*), ipairs(*), nsp(*), rem, mdloop
        character(len=4) :: ih(*)

        ! local variable declarations
        integer          :: holder_cut, holder_ntb, holder_natom, holder_nstlim, &
            holder_ntcnstph, holder_ntwr, holder_ntwx, &
            holder_ntpr, holder_ntt, holder_nscm, holder_ndfmin, &
            holder_irespa, natom3, i, selres_holder(1:TITR_RES_C), &
            holder_ntp, holder_mdloop, holder_rem
        logical          :: holder_tlangv
        type(state_rec) :: cph_ener
        _REAL_           :: holder_tcouple
        _REAL_, dimension(natom*3) :: vtemp

        ! This subroutine is the main driver for running constant pH MD in explicit
        ! solvent. The first thing it does is call cnstphbeginstep to set up the MC.
        ! Then it removes the periodic boundary conditions (PBC), selects a GB model,
        ! then calls force with nstep = 0 which will force oncpstep to be .true., so
        ! we'll get dvdl back. This is then passed to cnstphendstep to evaluate the
        ! transition. If it fails, then we restore the PBC, turn off the GB model,
        ! restore the CUToff, and return to the calling routine. If it's successful,
        ! we still return the above variables, but then we also turn on belly and fix
        ! the protein while we call runmd to relax the solvent for ntrelax steps.
        ! Some important things we do are:
        !  o  turn off T-coupling, since ntt=3 doesn't work with belly
        !  o  turn off trajectory/energy printing
        !  o  zero-out the v and vold arrays for the solute
        !  o  turn on belly and assign belly arrays

        ! Local variable assignments
        natom3 = natom*3
        holder_cut = cut
        holder_natom = natom
        holder_ntb = ntb
        holder_nstlim = nstlim
        holder_ntcnstph = ntcnstph
        holder_ntwr = ntwr
        holder_ntwx = ntwx
        holder_ntpr = ntpr
        holder_ntt = ntt
        holder_tlangv = tlangv
        holder_tcouple = gamma_ln
        holder_ndfmin = ndfmin
        holder_nscm = nscm
        holder_irespa = irespa
        holder_ntp = ntp
        holder_rem = rem
        holder_mdloop = mdloop

        cut = 9801.0d0    ! 99 * 99
        natom = cphfirst_sol - 1
        ntb = 0
        igb = cph_igb

        do i = 1, nres_titrate
            selres_holder(i) = 0
        end do

        ! Initiate constant pH stuff
        cph_success = .false.
        do i = 1, nres_titrate
            call excnstphbeginstep(xx(l190), selres_holder, i - 1)

            ! call force to get energies
            call force(xx, ix, ih, ipairs, x, f, cph_ener, cph_ener%vir, &
                xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
                do_list_update, 0)

            ! end constant pH stuff
            call cnstphendstep(xx(l190), xx(l15), cph_ener%pot%dvdl, temp0, solvph)

            ! We do _not_ want this force call to count toward our nmropt counter,
            ! so we decrement it here if nmropt is on
            if (nmropt /= 0) call nmrdcp()
        end do

        ! restore PBC variables
        cut = holder_cut
        natom = holder_natom
        ntb = holder_ntb
        igb = 0
        ntp = holder_ntp

        ! set iselres so that every residue is printed out
        iselres(0) = trescnt
        do i = 1, trescnt
            iselres(i) = i - 1
        end do

        ! write out the results to the cpout file
        if (master) call cnstphwrite(rem)

        ! if our exchange didn't succeed, just return and continue dynamics
        if (.not. cph_success) return

        ! increase the relaxation count
        relaxations = relaxations + 1
        tot_relax = tot_relax + 1

        ! if we succeeded, turn on relaxing and adjust other variables. Don't
        ! print during relaxation steps. Also turn off T-coupling and zero
        ! velocity array. Also, don't remove COM motion.
        relaxing = 1
        nstlim = ntrelax
        ntcnstph = ntrelax + 1
        ntwr = ntrelax + 1
        ntwx = ntrelax + 1
        ntpr = ntrelax + 1
        ntt = 0
        tlangv = .false.
        gamma_ln = 0.0
        do i = 1, cphfirst_sol*3 - 3
            vtemp(i) = v(i)
            v(i) = 0
        end do
        ndfmin = 0
        nscm = 0
        irespa = 1
        rem = 0
        mdloop = 0

        ! call runmd to relax the waters
        call runmd(xx, ix, ih, ipairs, x, winv, amass, f, &
            v, vold, xr, xc, conp, skip, nsp, tma, erstop, qsetup)

        ! restore the original values
        nstlim = holder_nstlim
        ntcnstph = holder_ntcnstph
        ntwr = holder_ntwr
        ntwx = holder_ntwx
        ntpr = holder_ntpr
        relaxing = 0
        tlangv = holder_tlangv
        gamma_ln = holder_tcouple
        ntt = holder_ntt
        ndfmin = holder_ndfmin
        nscm = holder_nscm
        irespa = holder_irespa
        rem = holder_rem
        mdloop = holder_mdloop
        do i = 1, cphfirst_sol*3 - 3
            v(i) = vtemp(i)
        end do

    end subroutine cnstph_explicitmd

end module constantph
