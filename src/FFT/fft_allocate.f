      subroutine fft_allocate(f,nfft,nsp)

!     allocates (and initializes to zero) the (3d) complex array f using
!     the fft dimensions nfft there is a check to make sure that these
!     values are products of 2,3, and 5 as required by the gpfa; if they
!     are not, they are converted to be valid dimensions.

!     note: lower dimensions of array = 1 (not 0)

!     the arrays for the g-space cutoffs are checked and set here as
!     necessary - another reason for allocating the arrays here

      use mod_comp
      use gfft,         only : l_valid_dim

      implicit none

      complex(kind=fp),  allocatable   :: f(:,:,:,:)
      integer,           intent(inout) :: nfft(3),nsp

      integer :: i,n
!     integer :: noptimal    ! function

!--->    check whether these dimensions are valid: done anytime the fft
!--->    dimensions change
      if( any( nfft /= nfft_old ) .or.                                   &
     &  (.not. allocated(nsphg) ) .or. (.not. allocated(nsphk) ) ) then
         do i=1,3
            if( .not. l_valid_dim( nfft(i) ) ) then
               n       = noptimal( nfft(i) )
               nfft(i) = n
            endif
            nfft_old(i) = nfft(i)
         enddo
         call fft_cutoffs( nfft(3) )
      endif

!--->    allocate and zero the array
      allocate( f( nfft(1), nfft(2), nfft(3), nsp ) )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      return
      end subroutine fft_allocate
