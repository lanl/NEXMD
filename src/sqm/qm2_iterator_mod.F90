module qm2_iterator_mod

    implicit none

    ! iterators / counters

    public :: scf_iterator_value
    public :: diis_iterator_value
    public :: diis_iterator_prev_value
    public :: remaining_diis_tokens

    private

contains

    function scf_iterator_value(qmmm_nml, ITER) result(val)
        ! this is an unbounded linear iterator
        use qmmm_nml_module, only : qmmm_nml_type
        implicit none

        type(qmmm_nml_type), intent(inout) :: qmmm_nml
        integer, intent(in), optional :: ITER
        integer :: val

        if (PRESENT(ITER)) qmmm_nml%iter_val_saved = ITER

        val = qmmm_nml%iter_val_saved

    end function scf_iterator_value

    function diis_iterator_value(qmmm_nml) result(val)
        ! this is a circular iterator

        use qmmm_nml_module, only : qmmm_nml_type
        implicit none
        type(qmmm_nml_type), intent(inout) :: qmmm_nml
        integer :: val

        val = MOD(scf_iterator_value(qmmm_nml) - 1, qmmm_nml%ndiis_matrices) + 1

    end function diis_iterator_value

    function diis_iterator_prev_value(qmmm_nml) result(val)
        ! this is a circular iterator

        use qmmm_nml_module, only : qmmm_nml_type
        implicit none

        type(qmmm_nml_type), intent(inout) :: qmmm_nml
        integer :: val

        val = MOD(scf_iterator_value(qmmm_nml) - 2, qmmm_nml%ndiis_matrices) + 1
        if (val < 1) val = 1

    end function diis_iterator_prev_value

    function remaining_diis_tokens(qmmm_nml, SET_TO) result(val)

        use qmmm_nml_module, only : qmmm_nml_type
        implicit none

        type(qmmm_nml_type), intent(inout) :: qmmm_nml
        integer, intent(in), optional :: SET_TO
        integer :: val

        if (PRESENT(SET_TO)) qmmm_nml%saved_val = SET_TO

        val = qmmm_nml%saved_val

    end function remaining_diis_tokens

end module qm2_iterator_mod
