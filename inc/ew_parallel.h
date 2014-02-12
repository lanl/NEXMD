! If you change the value of MPI_MAX_PROCESSORS
! Make sure you update parallel.h as well.
#undef MPI_MAX_PROCESSORS
#define MPI_MAX_PROCESSORS 256

integer BC_SLABS
parameter ( BC_SLABS = 9+4*(MPI_MAX_PROCESSORS+1) )

integer indz,ind_tr_tmp,ind_tr_tmp1,ntxyslab,ntxzslab
integer mxyslabs,mxzslabs
integer nxyslab(0:MPI_MAX_PROCESSORS)
integer nxzslab(0:MPI_MAX_PROCESSORS)
integer mxystart(0:MPI_MAX_PROCESSORS)
integer mxzstart(0:MPI_MAX_PROCESSORS)
integer num_recip,num_direct

common /ew_slabs/ indz,ind_tr_tmp,ind_tr_tmp1, &
      ntxyslab,ntxzslab,mxyslabs,mxzslabs, &
      nxyslab, nxzslab, mxystart,mxzstart, num_recip, num_direct

integer ranks(MPI_MAX_PROCESSORS)
integer direct_group,recip_group,world_group
logical i_do_recip,i_do_direct
common /pme_groups/direct_group,recip_group,world_group, &
                   i_do_recip,i_do_direct

integer direct_comm, recip_comm, world_comm
common /pme_comms/direct_comm, recip_comm, world_comm

