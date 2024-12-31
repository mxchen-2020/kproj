Module other_unfolding
        use mod_comp
        use mod_input
        contains

        subroutine other
        use mod_input
        use mod_lcao
        use mod_phonon
        use mod_metric
        call metric

        if (spin_texture .or. lcao_unfold == .true.) then
           call lcao
        endif

        if (phonon_unfold == .true.) then
           call unfold_phonon
        endif
        end subroutine


end Module

