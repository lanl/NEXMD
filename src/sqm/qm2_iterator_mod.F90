MODULE qm2_iterator_mod

  IMPLICIT NONE

  ! iterators / counters

  PUBLIC :: scf_iterator_value
  PUBLIC :: diis_iterator_value
  PUBLIC :: diis_iterator_prev_value
  PUBLIC :: remaining_diis_tokens

  PRIVATE

CONTAINS
  
  
  FUNCTION scf_iterator_value(ITER) RESULT(val)
    ! this is an unbounded linear iterator
    IMPLICIT NONE

    integer,intent(in),optional :: ITER
    integer :: val

    integer,save :: iter_val

    IF ( PRESENT(ITER) ) iter_val = ITER

    val = iter_val

  END FUNCTION scf_iterator_value


  FUNCTION diis_iterator_value() RESULT(val)
    ! this is a circular iterator

    use qmmm_module, only : qmmm_nml
    IMPLICIT NONE

    integer :: val
!    integer :: scf_iterator_value

    val = MOD( scf_iterator_value()-1 , qmmm_nml%ndiis_matrices ) + 1

  END FUNCTION diis_iterator_value



  FUNCTION diis_iterator_prev_value() RESULT(val)
    ! this is a circular iterator

    use qmmm_module, only : qmmm_nml
    IMPLICIT NONE

    integer :: val
!    integer :: scf_iterator_value

    val = MOD( scf_iterator_value()-2 , qmmm_nml%ndiis_matrices ) + 1
    IF ( val < 1 ) val = 1

  END FUNCTION diis_iterator_prev_value


  FUNCTION remaining_diis_tokens(SET_TO) RESULT(val)

    IMPLICIT NONE

    integer,intent(in),optional :: SET_TO
    integer :: val

    integer,save :: saved_val

    IF ( PRESENT(SET_TO) ) saved_val = SET_TO

    val = saved_val

  END FUNCTION remaining_diis_tokens


END MODULE qm2_iterator_mod
