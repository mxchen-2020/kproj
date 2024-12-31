      subroutine wf_to_fftmesh(f,nfft,nv,kv,z)

!     puts the wave function onto the fft mesh (for complex z)

      use mod_comp

      implicit none

      integer,     intent(in)    :: nfft(3)
      complex(DP), intent(inout) :: f(:,:,:)
      integer,     intent(in)    :: nv,kv(:,:)
      complex(DP), intent(in)    :: z(:)

      integer                         :: k,nf1,nf2,nf3

      f = cmplx( 0.0_DP, 0.0_DP, DP )
      do k=1,nv
         nf1 = modulo( kv(1,k)+nfft(1), nfft(1) ) + 1
         nf2 = modulo( kv(2,k)+nfft(2), nfft(2) ) + 1
         nf3 = modulo( kv(3,k)+nfft(3), nfft(3) ) + 1
         f(nf1,nf2,nf3) = z(k)
      enddo

      return
      end subroutine wf_to_fftmesh
