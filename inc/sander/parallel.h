!+ Specification and control of Amber's parallel implementation.

! The parallel programs use MPI in a replicated data paradigm.
! The usual algorithms to distribute coordinates and forces,
! which are activated when the number of processors are a power of two,
! employ a recursive doubling scheme based on ideas from the book
! by Fox, et al., and Robert van de Geijn's papers.

! If the number of cpus is not a power of 2 then a regular all gather / scatter
! etc is used.

! If you change MPI_MAX_PROCESSORS make sure you update logtwo
! below as well as ew_parallel.h
#undef  MPI_MAX_PROCESSORS
#define MPI_MAX_PROCESSORS 256
integer, dimension(1:MPI_MAX_PROCESSORS), parameter :: logtwo = &
   (/ 0,1,0,2,0,0,0,3,0,0,0,0,0,0,0,4, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8 &
      /)

integer commworld, commsander, commmaster, comm_lpimd
integer worldrank, sanderrank, masterrank, lpimd_rank
integer worldsize, sandersize, mastersize, lpimd_size
integer groupmaster, worldmaster
logical ng_sequential
integer numtasks,mytaskid
integer iparpt  ! the atom partition among the processors
! processor i owns atoms from iparpt(i) + 1 to iparpt(i+1)
integer iparpt3,rcvcnt,rcvcnt3
logical mpi_orig
dimension iparpt(0:MPI_MAX_PROCESSORS)
dimension iparpt3(0:MPI_MAX_PROCESSORS)
dimension rcvcnt(0:MPI_MAX_PROCESSORS)
dimension rcvcnt3(0:MPI_MAX_PROCESSORS)

integer notdone
common/parallel/numtasks,mytaskid,notdone, &
      iparpt,iparpt3,rcvcnt,rcvcnt3,mpi_orig
common/parallel_multi/commworld, commsander, commmaster, comm_lpimd, &
      worldrank, sanderrank, masterrank, lpimd_rank, &
      worldsize, sandersize, mastersize, lpimd_size, &
      groupmaster, worldmaster, ng_sequential

