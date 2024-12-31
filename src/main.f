       PROGRAM KPROJ
         ! Mingxing Chen University of Wisconsin, Milwaukee, USA, 2012
         use mod_comp
         use mod_input
         use mod_wave
         use mod_lattice
         use mod_fft
         use mod_kproj
         use mod_step_funciton
         use other_unfolding
         time_start = job_time()
         
         !other_unfold
         call input 
         IF(spin_texture.or.phonon_unfold.or.lcao_unfold==.true.) THEN
             call other
         ELSE !vasp-qe-abinit
         call init
         call init_wav
         call set_grid
         call strgn3

           ! do k-projection
           select case(task)
           case ("kproj")
             if(.not.LZLAYER.and..not.LXLAYER.and..not.LYLAYER)then
               call bands_wt_all
             else
               call bands_wt_layer
             endif
             call out_bands

           ! output wavefunction for a given band at a k-point
           case("psi2")
             call out_psi2
 
           ! output charge density for bands in [E1,E2]
           case("pchg")
             call out_pchg

           case default
             stop '   No task has been assigned'
         end select
         ENDIF
         time_end = job_time()
         consumetime = time_end - time_start
         call job_completed(consumetime)
       END
