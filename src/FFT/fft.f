      subroutine fft(f,transform,nfft)

      use mod_comp
      use gfft

      implicit none

!     calls the (modified) gpfa fft package of temperton.
!
!     now makes use of the fact that only spheres in reciprocal
!     space are used. this requires cutoffs and a distinction between
!     k_max (wave functions) and g_max (density and potentials)
!
!     f is the complex function to be transformed and is overwritten
!        with the result
!
!     transform (character variable) determines the type of transform:
!
!        'g->r'  reciprocal space to direct     space (g_max)
!        'k->r'  reciprocal space to direct     space (k_max)
!        '2k->r' reciprocal space to direct     space (2*k_max)
!        'G->r'  reciprocal space to direct     space (no cutoffs)
!
!        'r->g'  direct     space to reciprocal space (g_max)
!        'r->k'  direct     space to reciprocal space (2*k_max)
!        'r->G'  direct     space to reciprocal space (no cutoffs)
!
!
!     isig determines  direction of transform
!
!     isig = 0 is a set up call (not needed)
!          = 1 is reciprocal space to real space
!          =-1 is real       space to reciprocal space
!

      integer,          intent(in)    :: nfft(3)
      complex(DP), intent(inout) :: f(nfft(1),nfft(2),nfft(3))
      character(len=*), intent(in)    :: transform

      real(DP)    :: fac,t1,t2
      integer          :: i,id,nf2,is,isig,nn(3),npts

      nn(1) = nfft(1)
      nn(2) = nfft(2)
      nn(3) = nfft(3)
      id    = nfft(1)   ! can change to pad first dimension if needed
      nf2   = nfft(2)

      select case( transform )
         case( 'g->r' )  ! g_max sphere cutoff
            isig = 1
            call gpf3d_cut(f,id,nf2,nn,isig,nsphg)
         case( 'k->r' )  ! k_max sphere cutoff
            isig = 1
            call gpf3d_cut(f,id,nf2,nn,isig,nsphk)
!        case( '2k->r' )  ! 2*k_max sphere cutoff - for valence density
!           isig = 1
!           call gpf3d_cut(f,id,nf2,nn,isig,nsph2k)
         case( 'G->r' )  ! no use of cutoffs
            isig = 1
            call gpf3d(f,id,nf2,nn,isig)
         case( 'r->g' )  ! g_max sphere cutoff
            isig = -1
            call gpf3d_cut(f,id,nf2,nn,isig,nsphg)
!        case( 'r->k' )  ! 2*k_max sphere cutoff - for valence density
!           isig = -1
!           call gpf3d_cut(f,id,nf2,nn,isig,nsph2k)
         case( 'r->G' )  ! no use of cutoffs
            isig = -1
            call gpf3d(f,id,nf2,nn,isig)
         case default
            stop 'fft: invalid transform value'
      end select

      if (isig.eq.-1) then
         npts = nn(1)*nn(2)*nn(3)
         fac = 1.0_DP / npts
         f = fac*f
      endif

      return
      end subroutine fft
