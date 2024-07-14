      subroutine fft_cutoffs(nfft3)

!     generates the cutoffs for g-space on the fft mesh for the third
!     dimension for both stars and wave functions

      implicit none

      integer, intent(in) :: nfft3

      integer :: i3,k,k3

      if( allocated( nsphg ) ) deallocate( nsphg )
      if( allocated( nsphk ) ) deallocate( nsphk )

      allocate( nsphg(nfft3), nsphk(nfft3) )

      k3 = min( k3d, (nfft3-1)/2 )  ! needed if fftmesh smaller than
                                        ! star list

      nsphg = -1

      do i3 = -k3,k3
         if( nsphg_raw(i3) < 0 ) cycle
         k = mod( i3 + nfft3, nfft3 ) + 1
         nsphg(k) = nsphg_raw(i3)
      enddo

!     nsphk = -1

!     do i3 = -k3,k3
!        if( nsphk_raw(i3) < 0 ) cycle
!        k = mod( i3 + nfft3, nfft3 ) + 1
!        nsphk(k) = nsphk_raw(i3)
!     enddo

      return
      end subroutine fft_cutoffs
