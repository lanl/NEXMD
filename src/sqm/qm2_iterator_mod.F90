MODULE qm2_iterator_mod

  IMPLICIT NONE

  ! iterators / counters

  PUBLIC :: scf_iterator_value
  PUBLIC :: diis_iterator_value
  PUBLIC :: diis_iterator_prev_value
  PUBLIC :: remaining_diis_tokens

  PRIVATE

CONTAINS
  
  
  FUNCTION scf_iterator_value(qmmm_nml,ITER) RESULT(val)
    ! this is an unbounded linear iterator
    use qmmm_nml_module   , only : qmmm_nml_type
    IMPLICIT NONE

    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    integer,intent(in),optional :: ITER
    integer :: val

    IF ( PRESENT(ITER) ) qmmm_nml%iter_val_saved = ITER

    val = qmmm_nml%iter_val_saved

  END FUNCTION scf_iterator_value


  FUNCTION diis_iterator_value(qmmm_nml) RESULT(val)
    ! this is a circular iterator

  use qmmm_nml_module   , only : qmmm_nml_type
    IMPLICIT NONE
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    integer :: val

    val = MOD( scf_iterator_value(qmmm_nml)-1 , qmmm_nml%ndiis_matrices ) + 1

  END FUNCTION diis_iterator_value



  FUNCTION diis_iterator_prev_value(qmmm_nml) RESULT(val)
    ! this is a circular iterator

  use qmmm_nml_module   , only : qmmm_nml_type
    IMPLICIT NONE
    
    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    integer :: val

    val = MOD( scf_iterator_value(qmmm_nml)-2 , qmmm_nml%ndiis_matrices ) + 1
    IF ( val < 1 ) val = 1

  END FUNCTION diis_iterator_prev_value


  FUNCTION remaining_diis_tokens(qmmm_nml,SET_TO) RESULT(val)

  use qmmm_nml_module   , only : qmmm_nml_type
    IMPLICIT NONE

    type(qmmm_nml_type), intent(inout) :: qmmm_nml
    integer,intent(in),optional :: SET_TO
    integer :: val

    IF ( PRESENT(SET_TO) ) qmmm_nml%saved_val = SET_TO

    val = qmmm_nml%saved_val

  END FUNCTION remaining_diis_tokens


END MODULE qm2_iterator_mod
