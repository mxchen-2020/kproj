      module mod_fft
!
!     setting up fft mesh -- MW
!
      use mod_comp

!  3-d star mapping functions etc. and fft

!  number of 3-d stars
      integer, save :: ng3

!  dimensions for mapping function
      integer, save :: k1d,k2d,k3d

!  representative of star (lattice units)
      integer, save, allocatable :: kv3(:,:)  !  kv3(3,ng3)

!  length of star
      real(DP),    save, allocatable :: sk3(:) !  sk3(ng3)

!        cutoffs in reciprocal space to avoid calculating zeroes
!        this information is correct for gmax,kmax cutoffs
      integer, save, allocatable, private :: nsphg_raw(:)
      integer, save, allocatable, private :: nsphk_raw(:)
      integer, save,              private :: k3d_max


!     fft lengths
!         lengths corresponding to gmax (density, potential, etc.)
      integer, save, public               :: nfft_g(3) 

!         lengths to be used at any given time
      integer, save, public               :: nfft(3)
!         lengths used previously
      integer, save, private              :: nfft_old(3) = 0

!         star index such that g < 2*kmax
      integer, save, public               :: ng3_k

!        cutoffs - can change depending on the fft mesh used
      integer, save, allocatable, private :: nsphg(:),nsphk(:)



      contains


      subroutine strgn3
!********************************************************************
!     3-d star generator for general symmetry
!
!     for present purposes, assume no symmetry since not dealing with
!     stars, etc.
!     MW
!********************************************************************

      use mod_comp        
      use mod_lattice
      use mod_wave, only : gmax,NGX,NGY,NGZ

      implicit none

      integer               :: i,j,k,nng1,nng2,nng3,n0
      real(DP)              :: aa1,aa2,aa3,g,gcut2


!     real(DP), parameter   :: eps=1.e-6_DP

      integer,  allocatable :: ktmp(:,:)
      real(DP), allocatable :: stmp(:)

!--->    find optimal ng1, ng2, and ng3 (general case!)
!     aa1 = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
!     aa2 = sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 )
!     aa3 = sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 )

!     k1d = int( aa1*gmax/TWOPI + eps ) + 1
!     k2d = int( aa2*gmax/TWOPI + eps ) + 1
!     k3d = int( aa3*gmax/TWOPI + eps ) + 1

      k1d = NGX; k2d = NGY; k3d = NGZ 

      n0 = (2*k1d+1)*(2*k2d+1)*(2*k3d+1)                ! temporary arrays
      allocate( ktmp(3,n0), stmp(n0) )

!--->    allocate fft reciprocal space cutoff information
      allocate( nsphg_raw(-k3d:k3d) )
      allocate( nsphk_raw(-k3d:k3d) )
      nsphg_raw = -1
      nsphk_raw = -1

      gcut2 = gmax**2

!--->    loop over possible 3-d vectors

      ng3  = 0        !  number of stars
      nng1 = 0        !  max. index for b_1
      nng2 = 0        !  max. index for b_2
      nng3 = 0        !  max. index for b_3

      do i = k1d,-k1d,-1
         do j = k2d,-k2d,-1
            do k = k3d,-k3d,-1

!--->    determine length of candidate vector: g = i*b_1 + j*b_2 + k*b_3

               g = ( i*bmat(1,1) + j*bmat(2,1) + k*bmat(3,1) )**2        &
     &           + ( i*bmat(1,2) + j*bmat(2,2) + k*bmat(3,2) )**2        &
     &           + ( i*bmat(1,3) + j*bmat(2,3) + k*bmat(3,3) )**2

               if( g > gcut2 ) cycle  ! skip to next vector if too large

!--->    new vector (representative)

               ng3 = ng3+1
               stmp( ng3 ) = sqrt(g)
               ktmp(:,ng3) = (/ i,j,k /)
               nng1 = max( abs(i) , nng1 )
               nng2 = max( abs(j) , nng2 )
               nng3 = max( abs(k) , nng3 )
               nsphg_raw(k) = max( nsphg_raw(k) , abs(j) )

            enddo   ! k loop
         enddo      ! j loop
      enddo         ! i loop

!     do k = k3d,-k3d,-1
!       write(IUO,'(5X,A,2I4)')"k,nsphg_raw:",k,nsphg_raw(k)
!     enddo

!--->    generated all stars; now need to allocate other arrays and
!        collect all the needed information
      allocate( kv3(3,ng3) )
      allocate( sk3(ng3) )

      kv3(:,1:ng3) = ktmp(:,1:ng3)
      sk3(1:ng3) = stmp(1:ng3)

      write(IUO,'(/,3X,"Star generation:",                              &
     &         /,5X,"gmax              =",                              &
     &         f12.6,"  a.u.^-1    ",/,                                 &
     &         5X,"total g-vectors   =",i8)') gmax,ng3


!--->    fft lengths (nfft_g: density/potential; nfft_conv: convolutions)
      nfft_g(1)=noptimal(2*nng1+1)
      nfft_g(2)=noptimal(2*nng2+1)
      nfft_g(3)=noptimal(2*nng3+1)

      nfft = nfft_g

      ng3_k = ng3


      write(IUO,'(5X,"dimensions of g vectors       =",3i6)') k1d,k2d,k3d
      write(IUO,'(5X,"fft mesh: minimum mesh size   =",3i6)')           &
     &           2*nng1+1,2*nng2+1,2*nng3+1
      write(IUO,'(16x,"actual mesh size   =",3i6)')                        &
     &             (nfft_g(i),i=1,3)
      write(IUO,'(5X,A,I8)')"ng3_k:",ng3_k


      end subroutine strgn3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine set_nsphk(NPL_TMP,GV)

!     calculates the basis vectors corresponding to this k-point using
!     the same ordering as vasp

      use mod_comp

      implicit none

      integer,  intent(in)  :: NPL_TMP
      integer,  intent(in)  :: GV(3,NPL_TMP)       ! G-vectors

      integer               :: n,k,k1,k2,k3,n1,n2,n3
      real(DP)              :: g(3),g2,kcut2

!     write(IUO,'(/5X,A,3I5)')"NFFT MESH:",nfft(1:3)

      nsphk_raw = -1 ! reset the k-point limits

      do n=1,NPL_TMP
         k2 = GV(2,n)
         k3 = GV(3,n)
         nsphk_raw(k3) = max( nsphk_raw(k3), abs(k2) )
      enddo

      if( nfft(3) < 1 ) return
      if( allocated( nsphk) ) deallocate( nsphk )
      allocate( nsphk( nfft(3) ) )
      k3 = min( k3d, (nfft(3)-1)/2 )  ! needed if fftmesh smaller than
                                      ! list
      nsphk = -1

      do n3 = -k3,k3
         if( nsphk_raw(n3) < 0 ) cycle
         k = mod( n3 + nfft(3), nfft(3) ) + 1
         nsphk(k) = nsphk_raw(n3)
!       write(IUO,'(5X,A,2I4)')"k,nsphk:",k,nsphk(k)
      enddo

      return

      end subroutine set_nsphk

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      function noptimal(n)
!
!     determines the smallest integer m=(2**ip)*(3**iq)*(5**ir)
!     such that n <= m. this is needed for the fft meshes
!
      use mod_comp

      implicit none

      integer              :: noptimal
      integer, intent(in)  :: n
      integer              :: n2,n2m,n3,n3m,n5,n5m,nt2,nt23,ntest

      n2m = int( log(real(n,fp)) / log(2.0_DP) ) + 1
      n3m = int( log(real(n,fp)) / log(3.0_DP) ) + 1
      n5m = int( log(real(n,fp)) / log(5.0_DP) ) + 1

      noptimal=(2**n2m)

      do n2=n2m-1,0,-1
         nt2 = (2**n2)
         do n3=0,n3m
            nt23 = nt2*(3**n3)
            if( nt23 > noptimal ) exit
            do n5=0,n5m
               ntest = nt23*(5**n5)
               if( ntest > noptimal ) exit
               if( ntest >= n ) noptimal = min( noptimal, ntest )
            enddo
         enddo
      enddo

      return

      end function noptimal

      include 'FFT/fft_allocate.f'

      include 'FFT/fft_cutoffs.f'

      include 'FFT/fft.f'

      include 'FFT/wf_to_fftmesh.f'

      include 'FFT/cp_g_to_fftmesh.f'

      end module mod_fft
