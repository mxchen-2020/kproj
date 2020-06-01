      module gfft

!     fft using the gpfa of temperton

      use mod_comp

      private

      public :: gpf3d, gpf3d_cut, gpf2d, gpf2d_cut, l_valid_dim

!     note that the subroutines of the gpfa are outside of the module so
!     that type cheating (complex/real) is allowed. if these are wanted
!     in the module, then uncomment the reshape/transfer and use the
!     calls to gpfa with cr

!  array of trig functions/twiddle factors for gfft and dimensioning
      real(DP), save, private, allocatable :: trigs(:,:)
      integer,       save, private              :: ntrigs = -1
      integer,       save, private              :: nnold(3) = (/ 0,0,0 /)

      contains

!====================================================================

      subroutine GPF3D(c,id,nf2x,nn,is)
!     subroutine GPF3D(c,id,nf2x,nn,is,trigs,nfmx)
!
!     GPF3D: basic 3D complex FFT using the Generalized Prime Factor
!     Algorithm (Temperton, 1992).
!
!
!     Arguments
!
!     complex c(id,nf2x,nf3x) data array
!     id     first  dimension of data array.
!     nf2x   second dimension of data array
!     nn(3)  dimensions of transforms wanted.
!     is Forward (+1)/Reverse (-1)
!----------------------------------------------------------------------
!
!        definition of transform
!        -----------------------
!
!        x(j) = sum(k=0,...,n-1)(c(k)*exp(is*2*i*j*k*pi/n))
!
!----------------------------------------------------------------------
!     This version includes a number of modifications; in particular,
!     possible aliasing is removed. Note also that the ordering of
!     transforms is changed to allow for decoupling of dimensions
!     and the length of the transforms.                  m.w. 3-97

      use mod_comp

      implicit none

      integer,          intent(in)    :: id,nf2x,nn(3),is
      complex(DP), intent(inout) :: c(id,nf2x,nn(3))

      integer       :: i

!        if placed in module with rest the GPF terms, need to explicitly
!        change shape and transfer complex to real; not needed if module
!        ends before gpfa
!     real(DP) :: cr(2*id*nf2x*nn(3))
!     cr = reshape( transfer( c,cr ) , (/ 2*id*nf2x*nn(3) /) )

      if( any( nnold /= nn ) ) then
         ntrigs = ntrig_dim( nn )
         if( allocated(trigs) ) deallocate( trigs )
         allocate( trigs( ntrigs,3 ) )
         nnold = nn
         do i=1,3
            call setgpfa(trigs(1,i), nn(i))
         enddo
      endif

      if(is.eq.0) return

!
!  The following three loops are prime candidates for parallelism
!  as they call independent multi-1D ffts.
!

      do i=1,nn(3)
         call gpfa(c(1,1,i),trigs(1,1),2,2*id,nn(1),nn(2),is)
!        call gpfa(cr(2*((i-1)*id*nf2x)+1),trigs(1,1),2,2*id,nn(1),nn(2),is)
      enddo

      do i=1,nn(3)
         call gpfa(c(1,1,i),trigs(1,2),2*id,2,nn(2),nn(1),is)
!        call gpfa(cr(2*((i-1)*id*nf2x)+1),trigs(1,1),2,2*id,nn(1),nn(2),is)
      enddo

      do i=1,nn(2)
         call gpfa(c(1,i,1),trigs(1,3),2*nf2x*id,2,nn(3),nn(1),is)
!        call gpfa(cr(2*((i-1)*id)+1),trigs(1,1),2,2*id,nn(1),nn(2),is)
      enddo

      return
      end subroutine GPF3D

!====================================================================

      subroutine GPF3D_cut(c,id,nf2x,nn,is,nfz)
!     subroutine GPF3D_cut(c,id,nf2x,nn,is,trigs,nfmx,nfz)
!
!     GPF3D: basic 3D complex FFT using the Generalized Prime Factor
!     Algorithm (Temperton, 1992).
!
!     9-99 (mw): allows cutoffs in reciprocal space to avoid
!                transforming vectors of zeroes.
!
!     Arguments
!
!     complex c(id,nf2x,nf3x) data array
!     id     first  dimension of data array.
!     nf2x   second dimension of data array
!     nn(3)  dimensions of transforms wanted.
!     is Forward (+1)/Reverse (-1)
!
!     nfz(nf3x)  reciprocal space cutoffs for each value of third index
!
!----------------------------------------------------------------------
!
!        definition of transform
!        -----------------------
!
!        x(j) = sum(k=0,...,n-1)(c(k)*exp(is*2*i*j*k*pi/n))
!
!----------------------------------------------------------------------
!     This version includes a number of modifications; in particular,
!     possible aliasing is removed. Note also that the ordering of
!     transforms is changed to allow for decoupling of dimensions
!     and the length of the transforms.                  m.w. 3-97

      use mod_comp

      implicit none

      integer,          intent(in)    :: id,nf2x,nn(3),is
      complex(DP), intent(inout) :: c(id,nf2x,nn(3))

!--->   array giving cutoff in reciprocal space for each third index
      integer,          intent(in)    :: nfz(*)

      integer       :: i,j,m2

!        if placed in module with rest the GPF terms, need to explicitly
!        change shape and transfer complex to real; not needed if module
!        ends before gpfa
!     real(DP) :: cr(2*id*nf2x*nn(3))
!     cr = reshape( transfer( c,cr ) , (/ 2*id*nf2x*nn(3) /) )

      if( any( nnold /= nn ) ) then
         ntrigs = ntrig_dim( nn )
         if( allocated(trigs) ) deallocate( trigs )
         allocate( trigs( ntrigs,3 ) )
         nnold = nn
         do i=1,3
            call setgpfa(trigs(1,i), nn(i))
         enddo
      endif

      if(is.eq.0) return

!
!  The following three loops are prime candidates for parallelism
!  as they call independent multi-1D ffts.
!

      if(is.eq.1) then  ! =========>> reciprocal to real space

      do i=1,nn(3)

         if(nfz(i).ge.0) then

!!       do j = 1,nfz(i)+1  ! nn(2)
!!          call gpfa(c(1,j,i),trigs(1,1),2,2,nn(1),1,is)
!!       enddo
         m2 = min( nn(2), nfz(i)+1 )
         call gpfa(c(1,1,i),trigs(1,1),2,2*id,nn(1),m2,is)
!        call gpfa(cr(2*((i-1)*id*nf2x)+1),trigs(1,1),2,2*id,nn(1),nfz(i)+1,is)

!!       do j = nn(2)-nfz(i)+1,nn(2)
!!          call gpfa(c(1,j,i),trigs(1,1),2,2,nn(1),1,is)
!!       enddo
         if( nfz(i) > 0 .and. m2 < nn(2) ) then
!           j = nn(2)-nfz(i)+1
            j = max( nn(2) - nfz(i) + 1, m2 + 1 )
            m2 = nn(2) - j + 1
            call gpfa(c(1,j,i),trigs(1,1),2,2*id,nn(1),m2,is)
!           call gpfa(cr(2*(i-1)*id*nf2x+2*(j-1)*id+1),trigs(1,1),2,2*id,nn(1),nfz(i),is)
         endif

!!       do j = 1,nn(1)
!!          call gpfa(c(j,1,i),trigs(1,2),2*id,2,nn(2),1,is)
!!       enddo
         call gpfa(c(1,1,i),trigs(1,2),2*id,2,nn(2),nn(1),is)
!        call gpfa(cr(2*(i-1)*id*nf2x+1),trigs(1,2),2*id,2,nn(2),nn(1),is)

      endif

      enddo

      do i=1,nn(2)
!!       do j = 1,nn(1)
!!          call gpfa(c(j,i,1),trigs(1,3),2*nf2x*id,2,nn(3),1,is)
!!       enddo
         call gpfa(c(1,i,1),trigs(1,3),2*nf2x*id,2,nn(3),nn(1),is)
!        call gpfa(cr(2*(i-1)*id+1),trigs(1,3),2*nf2x*id,2,nn(3),nn(1),is)
      enddo


      else  ! ==================>> real to reciprocal space

!--->    reverse order of transforms
      do i=1,nn(2)
!!       do j = 1,nn(1)
!!          call gpfa(c(j,i,1),trigs(1,3),2*nf2x*id,2,nn(3),1,is)
!!       enddo
         call gpfa(c(1,i,1),trigs(1,3),2*nf2x*id,2,nn(3),nn(1),is)
!        call gpfa(cr(2*(i-1)*id+1),trigs(1,3),2*nf2x*id,2,nn(3),nn(1),is)
      enddo

      do i=1,nn(3)
         if( nfz(i).lt.0 ) cycle

!!       do j = 1,nn(1)
!!          call gpfa(c(j,1,i),trigs(1,2),2*id,2,nn(2),1,is)
!!       enddo
         call gpfa(c(1,1,i),trigs(1,2),2*id,2,nn(2),nn(1),is)
!        call gpfa(cr(2*(i-1)*id*nf2x+1),trigs(1,2),2*id,2,nn(2),nn(1),is)


!!       do j = 1,nfz(i)+1
!!          call gpfa(c(1,j,i),trigs(1,1),2,2,nn(1),1,is)
!!       enddo
         m2 = min( nn(2), nfz(i)+1 )
         call gpfa(c(1,1,i),trigs(1,1),2,2*id,nn(1),m2,is)
!        call gpfa(cr(2*(i-1)*id*nf2x+1),trigs(1,1),2,2*id,nn(1),nfz(i)+1,is)

!!       do j = nn(2)-nfz(i)+1,nn(2)
!!          call gpfa(c(1,j,i),trigs(1,1),2,2,nn(1),1,is)
!!       enddo
         if( nfz(i) > 0 .and. m2 < nn(2) ) then
!           j = nn(2)-nfz(i)+1
            j = max( nn(2) - nfz(i) + 1, m2 + 1 )
            m2 = nn(2) - j + 1
            call gpfa(c(1,j,i),trigs(1,1),2,2*id,nn(1),m2,is)
!           call gpfa(cr(2*(i-1)*id*nf2x+2*(j-1)*id+1),trigs(1,1),2,2*id,nn(1),nfz(i),is)
         endif

      enddo

      endif

      return
      end subroutine GPF3D_cut

!====================================================================

      subroutine gpf2d(c,id,nf2x,nn,is)
!
!     basic 2D complex FFT using the Generalized Prime Factor
!     Algorithm (Temperton, 1992).
!     (modified from gpf3d)
!
!     Arguments
!
!     complex c(id,nf2x) data array
!     id     first  dimension of data array.
!     nf2x   second dimension of data array
!     nn(2)  dimensions of transforms wanted.
!     is Forward (+1)/Reverse (-1)
!----------------------------------------------------------------------

      use mod_comp

      implicit none

      integer,          intent(in)    :: id,nf2x,nn(2),is
      complex(DP), intent(inout) :: c(id,nf2x)

      integer       :: i

!        if placed in module with rest the GPF terms, need to explicitly
!        change shape and transfer complex to real; not needed if module
!        ends before gpfa
!     real(DP) :: cr(2*id*nf2x)
!     cr = reshape( transfer( c,cr ) , (/ 2*id*nf2x /) )

      if( any( nnold(1:2) /= nn(1:2) ) ) then
         ntrigs = ntrig_dim( nn )
         if( allocated(trigs) ) deallocate( trigs )
         allocate( trigs( ntrigs,3 ) )
         nnold(1:2) = nn(1:2)
         nnold(3)   = 0
         do i=1,2
            call setgpfa(trigs(1,i), nn(i))
         enddo
      endif

      if(is.eq.0) return

      call gpfa(c,trigs(1,1),2,2*id,nn(1),nn(2),is)
      call gpfa(c,trigs(1,2),2*id,2,nn(2),nn(1),is)

      return
      end subroutine gpf2d

!====================================================================

      subroutine gpf2d_cut(c,id,nf2x,nn,is,nfz)
!
!     basic 2D complex FFT using the Generalized Prime Factor
!     Algorithm (Temperton, 1992). modification of gpf3d.
!
!     9-99 (mw): allows cutoffs in reciprocal space to avoid
!                transforming vectors of zeroes.
!
!     Arguments
!
!     complex c(id,nf2x) data array
!     id     first  dimension of data array.
!     nf2x   second dimension of data array
!     nn(2)  dimensions of transforms wanted.
!     is Forward (+1)/Reverse (-1)
!
!     nfz  reciprocal space cutoff (must be > 0 or else not 2-D stars)

      use mod_comp

      implicit none

      integer,          intent(in)    :: id,nf2x,nn(2),is
      complex(DP), intent(inout) :: c(id,nf2x)

!--->   array giving cutoff in reciprocal space for each third index
      integer,          intent(in)    :: nfz

      integer       :: i,j,m2

!        if placed in module with rest the GPF terms, need to explicitly
!        change shape and transfer complex to real; not needed if module
!        ends before gpfa
!     real(DP) :: cr(2*id*nf2x))
!     cr = reshape( transfer( c,cr ) , (/ 2*id*nf2x /) )

      if( nfz < 0 ) stop 'gpf2d_cut error: nfz < 0'

      if( any( nnold(1:2) /= nn(1:2) ) ) then
         ntrigs = ntrig_dim( nn )
         if( allocated(trigs) ) deallocate( trigs )
         allocate( trigs( ntrigs,3 ) )
         nnold(1:2) = nn(1:2)
         nnold(3)   = 0
         do i=1,2
            call setgpfa(trigs(1,i), nn(i))
         enddo
      endif

      if(is.eq.0) return


      if(is.eq.1) then  ! =========>> reciprocal to real space


         m2 = min( nfz + 1, nn(2) )
         call gpfa(c,trigs(1,1),2,2*id,nn(1),m2,is)
         if( nfz > 0 .and. m2 < nn(2) ) then
            j = max( nn(2)-nfz+1, m2 + 1 )
            m2 = nn(2) - j + 1
            call gpfa(c(1,j),trigs(1,1),2,2*id,nn(1),m2,is)
         endif
         call gpfa(c,trigs(1,2),2*id,2,nn(2),nn(1),is)


      else  ! ==================>> real to reciprocal space

!--->    reverse order of transforms

         call gpfa(c,trigs(1,2),2*id,2,nn(2),nn(1),is)
         m2 = min( nfz + 1, nn(2) )
         call gpfa(c,trigs(1,1),2,2*id,nn(1),m2,is)
         if( nfz > 0 .and. m2 < nn(2) ) then
            j = max( nn(2)-nfz+1, m2 + 1 )
            m2 = nn(2) - j + 1
            call gpfa(c(1,j),trigs(1,1),2,2*id,nn(1),m2,is)
         endif


      endif

      return
      end subroutine gpf2d_cut

!====================================================================

      subroutine setgpfa(trigs,n)
!
!        subroutine 'setgpfa'
!        setup routine for self-sorting in-place
!            generalized prime factor (complex) fft [gpfa]
!
!        call setgpfa(trigs,n)
!
!        input :
!        -----
!        n is the length of the transforms. n must be of the form:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!
!        output:
!        ------
!        trigs is a table of twiddle factors,
!          of length 2*ipqr (real) words, where:
!          --------------------------------------
!            ipqr = (2**ip) + (3**iq) + (5**ir)
!          --------------------------------------
!
!        written by clive temperton 1990
!
!----------------------------------------------------------------------
!
      use mod_comp

      implicit none

      real(DP), intent(out) :: trigs(*)
      integer,       intent(in)  :: n

      real(DP) :: angle,del
      integer       :: i,ifac,ip,iq,ir,irot,k,kink,kk,ll,ni,nj(3),nn

      character(len=80) :: errormessage

!     real(DP), parameter :: twopi = 6.283185307179586477_DP
!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do 30 ll = 1 , 3
      kk = 0
   10 continue
      if (mod(nn,ifac).ne.0) go to 20
      kk = kk + 1
      nn = nn / ifac
      go to 10
   20 continue
      nj(ll) = kk
      ifac = ifac + ll
   30 continue
!
      if (nn.ne.1) then
         write(*,'(" gpfa:",i0," not a legal value of n")') n
         stop  'gpfa: illegal value of n'
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute list of rotated twiddle factors
!     ---------------------------------------
      nj(1) = 2**ip
      nj(2) = 3**iq
      nj(3) = 5**ir
!
      i = 1
!
      do 60 ll = 1 , 3
      ni = nj(ll)
      if (ni.eq.1) go to 60
!
      del = TWOPI / real(ni,DP)
      irot = n / ni
      kink = mod(irot,ni)
      kk = 0
!
      do 50 k = 1 , ni
!     angle = kk * del
      angle = (kk * TWOPI)/real(ni,DP)

      trigs(i) = cos(angle)
      trigs(i+1) = sin(angle)
      i = i + 2
      kk = kk + kink
      if (kk.gt.ni) kk = kk - ni
   50 continue
   60 continue
!
      return
      end subroutine setgpfa

!====================================================================


      function ntrig_dim(n)
!
!     determines the decomposition m=(2**ip)*(3**iq)*(5**ir) of n
!     and determines the required dimension for twiddle factors (trigs)
!

      implicit none

      integer              :: ntrig_dim
      integer, intent(in)  :: n(:)
      integer              :: nn
      integer              :: i,ip,iq,ir

      ntrig_dim = 0
      do i=1,size(n)  ! loop over the fft dimensions

      nn = n(i)
      
!--->    get powers of 2
      ip = 0
      do 
         if( mod( nn,2) /= 0 ) exit
         ip = ip + 1
         nn = nn/2
      enddo
      
!--->    get powers of 3
      iq = 0
      do 
         if( mod( nn,3) /= 0 ) exit
         iq = iq + 1
         nn = nn/3
      enddo
      
!--->    get powers of 5
      ir = 0
      do 
         if( mod( nn,5) /= 0 ) exit
         ir = ir + 1
         nn = nn/5
      enddo

      if( nn /= 1 ) then
         stop "error setgpfa: n not a product of 2, 3, 5"
      endif

      ntrig_dim = max( ntrig_dim, 2* (2**ip + 3**iq + 5**ir ) )

      enddo ! loop over fft directions

      return

      end function ntrig_dim

!====================================================================


      function l_valid_dim(n)
!
!     determines whether n is an integer of the form
!     n=(2**ip)(3**iq)(5**ir) as required for gpfa
!

      implicit none

      logical              :: l_valid_dim
      integer, intent(in)  :: n
      integer              :: nn

      l_valid_dim = .false.

      if( n <= 0 ) return

      nn = n
      
!--->    get powers of 2
      do 
         if( mod( nn,2) /= 0 ) exit
         nn = nn/2
      enddo
      
!--->    get powers of 3
      do 
         if( mod( nn,3) /= 0 ) exit
         nn = nn/3
      enddo
      
!--->    get powers of 5
      do 
         if( mod( nn,5) /= 0 ) exit
         nn = nn/5
      enddo

      if( nn == 1 ) then
         l_valid_dim = .true.
      endif

      return
      end function l_valid_dim

      end module gfft


!====================================================================

      subroutine gpfa(c,trigs,inc,jump,n,lot,isign)
!
!        subroutine 'gpfa'
!        self-sorting in-place generalized prime factor (complex) fft
!
!        *** this is the all-fortran version ***
!            -------------------------------
!
!        was:  call gpfa(a,b,trigs,inc,jump,n,lot,isign)
!        now:  call gpfa(c,trigs,inc,jump,n,lot,isign)
!
!        a is first real input/output vector
!        b is first imaginary input/output vector
!        ****  replaced with single vector c *****
!
!        trigs is a table of twiddle factors, precalculated
!              by calling subroutine 'setgpfa'
!        inc is the increment within each data vector
!        jump is the increment between data vectors
!        n is the length of the transforms:
!          -----------------------------------
!            n = (2**ip) * (3**iq) * (5**ir)
!          -----------------------------------
!        lot is the number of transforms
!        isign = +1 for forward transform
!              = -1 for inverse transform
!
!        written by clive temperton
!        recherche en prevision numerique
!        atmospheric environment service, canada
!
!----------------------------------------------------------------------
!
!        definition of transform
!        -----------------------
!
!        x(j) = sum(k=0,...,n-1)(c(k)*exp(isign*2*i*j*k*pi/n))
!
!---------------------------------------------------------------------
!
!        for a mathematical development of the algorithm used,
!        see:
!
!        c temperton : "a generalized prime factor fft algorithm
!          for any n = (2**p)(3**q)(5**r)",
!          siam j. sci. stat. comp., may 1992.
!
!----------------------------------------------------------------------
!        modified version (see gpf3d.f)                 m.w.  3-97
!
      use mod_comp

      implicit none

      real(DP), intent(inout) :: c(*)
      real(DP), intent(in)    :: trigs(*)
      integer,       intent(in)    :: inc,jump,n,lot,isign

      integer :: nj(3)
      integer :: i,ifac,ip,iq,ir,kk,ll,nn

      character(len=80) :: errormessage

!
!     decompose n into factors 2,3,5
!     ------------------------------
      nn = n
      ifac = 2
!
      do 30 ll = 1 , 3
      kk = 0
   10 continue
      if (mod(nn,ifac).ne.0) go to 20
      kk = kk + 1
      nn = nn / ifac
      go to 10
   20 continue
      nj(ll) = kk
      ifac = ifac + ll
   30 continue
!
      if (nn.ne.1) then
         write(*,'(" gpfa:",i0," not a legal value of n")') n
         stop  'gpfa: illegal value of n'
      endif
!
      ip = nj(1)
      iq = nj(2)
      ir = nj(3)
!
!     compute the transform
!     ---------------------
      i = 1
      if (ip.gt.0) then
         call gpfa2f(c,trigs,inc,jump,n,ip,lot,isign)
         i = i + 2 * ( 2**ip)
      endif
      if (iq.gt.0) then
         call gpfa3f(c,trigs(i),inc,jump,n,iq,lot,isign)
         i = i + 2 * (3**iq)
      endif
      if (ir.gt.0) then
         call gpfa5f(c,trigs(i),inc,jump,n,ir,lot,isign)
      endif
!
      return
      end subroutine gpfa

      subroutine gpfa2f(c,trigs,inc,jump,n,mm,lot,isign)
!     subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
!
!     fortran version of *gpfa2* -
!     radix-2 section of self-sorting, in-place, generalized pfa
!     central radix-2 and radix-8 passes included
!      so that transform length can be any power of 2
!
!-------------------------------------------------------------------
!
      use mod_comp

      implicit none

      real(DP), intent(inout) :: c(*)
      real(DP), intent(in)    :: trigs(*)
      integer,       intent(in)    :: inc,isign,jump,lot,mm,n

      real(DP) :: aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji,ajj
      real(DP) :: ajk,ajl,ajm,ajn,ajo,ajp,bja,bjb,bjc,bjd
      real(DP) :: bje,bjf,bjg,bjh,bji,bjj,bjk,bjl,bjm,bjn
      real(DP) :: bjo,bjp,c1,c2,c3,co1,co2,co3,co4,co5,co6,co7
      real(DP) :: s,si1,si2,si3,si4,si5,si6,si7,ss
      real(DP) :: t0,t1,t2,t3,u0,u1,u2,u3
      integer :: ink,inq,ipass,istart,j,ja,jb,jc,jd,je
      integer :: jf,jg,jh,ji,jj,jjj,jk,jl,jm,jn,jo,jp,jstep,jstepl
      integer :: jstepx,k,kk,l,la,laincl,left,ll,m,m2
      integer :: m8,mh,mu,n2,nb,nblox,ninc,nu,nvex

      integer, parameter :: lvr=64
!
!     ***************************************************************
!     *                                                             *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *                                                             *
!     ***************************************************************
!
      n2 = 2**mm
      inq = n/n2
      jstepx = (n2-n) * inc
      ninc = n * inc
      ink = inc * inq
!
      m2 = 0
      m8 = 0
      if (mod(mm,2).eq.0) then
         m = mm/2
      else if (mod(mm,4).eq.1) then
         m = (mm-1)/2
         m2 = 1
      else if (mod(mm,4).eq.3) then
         m = (mm-3)/2
         m8 = 1
      endif
      mh = (m+1)/2
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = real(isign,DP)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type I radix-4 passes
!  -----------------------------
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0_DP
      if (mu.eq.3) ss = -1.0_DP
!
      if (mh.eq.0) go to 200
!
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 4*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 110 l = 1 , nvex
      aja = c(ja+j)
      ajc = c(jc+j)
      t0 = aja + ajc
      t2 = aja - ajc
      ajb = c(jb+j)
      ajd = c(jd+j)
      t1 = ajb + ajd
      t3 = ss * ( ajb - ajd )
      bja = c(ja+j+1)
      bjc = c(jc+j+1)
      u0 = bja + bjc
      u2 = bja - bjc
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u1 = bjb + bjd
      u3 = ss * ( bjb - bjd )
      c(ja+j)   = t0 + t1
      c(jc+j)   = t0 - t1
      c(ja+j+1) = u0 + u1
      c(jc+j+1) = u0 - u1
      c(jb+j)   = t2 - u3
      c(jd+j)   = t2 + u3
      c(jb+j+1) = u2 + t3
      c(jd+j+1) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
!
!  finished if n2 = 4
!  ------------------
      if (n2.eq.4) go to 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 4*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 130 l = 1 , nvex
      aja = c(ja+j)
      ajc = c(jc+j)
      t0 = aja + ajc
      t2 = aja - ajc
      ajb = c(jb+j)
      ajd = c(jd+j)
      t1 = ajb + ajd
      t3 = ss * ( ajb - ajd )
      bja = c(ja+j+1)
      bjc = c(jc+j+1)
      u0 = bja + bjc
      u2 = bja - bjc
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u1 = bjb + bjd
      u3 = ss * ( bjb - bjd )
      c(ja+j)   = t0 + t1
      c(ja+j+1) = u0 + u1
      c(jb+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jb+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jc+j)   = co2*(t0-t1) - si2*(u0-u1)
      c(jc+j+1) = si2*(t0-t1) + co2*(u0-u1)
      c(jd+j)   = co3*(t2+u3) - si3*(u2-t3)
      c(jd+j+1) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  130 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 4*la
  160 continue
!-----( end of loop on type I radix-4 passes)
!
!  central radix-2 pass
!  --------------------
  200 continue
      if (m2.eq.0) go to 300
!
      jstep = (n*inc) / (2*la)
      jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 220 jjj = 0 , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 215 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 210 l = 1 , nvex
      aja = c(ja+j)
      ajb = c(jb+j)
      t0 = aja - ajb
      c(ja+j) = aja + ajb
      c(jb+j) = t0
      bja = c(ja+j+1)
      bjb = c(jb+j+1)
      u0 = bja - bjb
      c(ja+j+1) = bja + bjb
      c(jb+j+1) = u0
      j = j + jump
  210 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  215 continue
  220 continue
!
!  finished if n2=2
!  ----------------
      if (n2.eq.2) go to 490
!
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 260 k = ink , jstep - ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
!
!  loop along transforms
!  ---------------------
      do 250 jjj = k , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 245 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (kk.eq.n2/2) then
      do 230 l = 1 , nvex
      aja = c(ja+j)
      ajb = c(jb+j)
      t0 = ss * ( aja - ajb )
      c(ja+j) = aja + ajb
      bjb = c(jb+j+1)
      bja = c(ja+j+1)
      c(jb+j)   = ss * ( bjb - bja )
      c(ja+j+1) = bja + bjb
      c(jb+j+1) = t0
      j = j + jump
  230 continue
!
      else
!
      do 240 l = 1 , nvex
      aja = c(ja+j)
      ajb = c(jb+j)
      t0 = aja - ajb
      c(ja+j)   = aja + ajb
      bja = c(ja+j+1)
      bjb = c(jb+j+1)
      u0 = bja - bjb
      c(ja+j+1) = bja + bjb
      c(jb+j)   = co1*t0 - si1*u0
      c(jb+j+1) = si1*t0 + co1*u0
      j = j + jump
  240 continue
!
      endif
!
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  245 continue
  250 continue
!-----(end of loop along transforms)
      kk = kk + 2 * la
  260 continue
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
      la = 2 * la
      go to 400
!
!  central radix-8 pass
!  --------------------
  300 continue
      if (m8.eq.0) go to 400
      jstep = (n*inc) / (8*la)
      jstepl = jstep - ninc
      mu = mod(inq,8)
      if (isign.eq.-1) mu = 8 - mu
      c1 = 1.0_DP
      if (mu.eq.3.or.mu.eq.7) c1 = -1.0_DP
      c2 = sqrt(0.5_DP)
      if (mu.eq.3.or.mu.eq.5) c2 = -c2
      c3 = c1 * c2
!
!  stage 1
!  -------
      do 320 k = 0 , jstep - ink , ink
      do 315 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 312 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
      do 310 l = 1 , nvex
      aja = c(ja+j)
      aje = c(je+j)
      t0 = aja - aje
      c(ja+j) = aja + aje
      ajc = c(jc+j)
      ajg = c(jg+j)
      t1 = c1 * ( ajc - ajg )
      c(je+j) = ajc + ajg
      ajb = c(jb+j)
      ajf = c(jf+j)
      t2 = ajb - ajf
      c(jc+j) = ajb + ajf
      ajd = c(jd+j)
      ajh = c(jh+j)
      t3 = ajd - ajh
      c(jg+j) = ajd + ajh
      c(jb+j) = t0
      c(jf+j) = t1
      c(jd+j) = c2 * ( t2 - t3 )
      c(jh+j) = c3 * ( t2 + t3 )
      bja = c(ja+j+1)
      bje = c(je+j+1)
      u0 = bja - bje
      c(ja+j+1) = bja + bje
      bjc = c(jc+j+1)
      bjg = c(jg+j+1)
      u1 = c1 * ( bjc - bjg )
      c(je+j+1) = bjc + bjg
      bjb = c(jb+j+1)
      bjf = c(jf+j+1)
      u2 = bjb - bjf
      c(jc+j+1) = bjb + bjf
      bjd = c(jd+j+1)
      bjh = c(jh+j+1)
      u3 = bjd - bjh
      c(jg+j+1) = bjd + bjh
      c(jb+j+1) = u0
      c(jf+j+1) = u1
      c(jd+j+1) = c2 * ( u2 - u3 )
      c(jh+j+1) = c3 * ( u2 + u3 )
      j = j + jump
  310 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  312 continue
  315 continue
  320 continue
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
      do 330 jjj = 0 , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 328 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
      do 325 l = 1 , nvex
      aja = c(ja+j)
      aje = c(je+j)
      t0 = aja + aje
      t2 = aja - aje
      ajc = c(jc+j)
      ajg = c(jg+j)
      t1 = ajc + ajg
      t3 = c1 * ( ajc - ajg )
      bja = c(ja+j+1)
      bje = c(je+j+1)
      u0 = bja + bje
      u2 = bja - bje
      bjc = c(jc+j+1)
      bjg = c(jg+j+1)
      u1 = bjc + bjg
      u3 = c1 * ( bjc - bjg )
      c(ja+j)   = t0 + t1
      c(je+j)   = t0 - t1
      c(ja+j+1) = u0 + u1
      c(je+j+1) = u0 - u1
      c(jc+j)   = t2 - u3
      c(jg+j)   = t2 + u3
      c(jc+j+1) = u2 + t3
      c(jg+j+1) = u2 - t3
      ajb = c(jb+j)
      ajd = c(jd+j)
      t0 = ajb + ajd
      t2 = ajb - ajd
      ajf = c(jf+j)
      ajh = c(jh+j)
      t1 = ajf - ajh
      t3 = ajf + ajh
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u0 = bjb + bjd
      u2 = bjb - bjd
      bjf = c(jf+j+1)
      bjh = c(jh+j+1)
      u1 = bjf - bjh
      u3 = bjf + bjh
      c(jb+j)   = t0 - u3
      c(jh+j)   = t0 + u3
      c(jb+j+1) = u0 + t3
      c(jh+j+1) = u0 - t3
      c(jd+j)   = t2 + u1
      c(jf+j)   = t2 - u1
      c(jd+j+1) = u2 - t1
      c(jf+j+1) = u2 + t1
      j = j + jump
  325 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  328 continue
  330 continue
!
      if (n2.eq.8) go to 490
!
!  loop on nonzero k
!  -----------------
      kk = 2 * la
!
      do 350 k = ink , jstep - ink , ink
!
      co1 = trigs(kk+1)
      si1 = s * trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s * trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s * trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s * trigs(4*kk+2)
      co5 = trigs(5*kk+1)
      si5 = s * trigs(5*kk+2)
      co6 = trigs(6*kk+1)
      si6 = s * trigs(6*kk+2)
      co7 = trigs(7*kk+1)
      si7 = s * trigs(7*kk+2)
!
      do 345 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 342 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
      do 340 l = 1 , nvex
      aja = c(ja+j)
      aje = c(je+j)
      t0 = aja + aje
      t2 = aja - aje
      ajc = c(jc+j)
      ajg = c(jg+j)
      t1 = ajc + ajg
      t3 = c1 * ( ajc - ajg )
      bja = c(ja+j+1)
      bje = c(je+j+1)
      u0 = bja + bje
      u2 = bja - bje
      bjc = c(jc+j+1)
      bjg = c(jg+j+1)
      u1 = bjc + bjg
      u3 = c1 * ( bjc - bjg )
      c(ja+j)   = t0 + t1
      c(ja+j+1) = u0 + u1
      c(je+j)   = co4*(t0-t1) - si4*(u0-u1)
      c(je+j+1) = si4*(t0-t1) + co4*(u0-u1)
      c(jc+j)   = co2*(t2-u3) - si2*(u2+t3)
      c(jc+j+1) = si2*(t2-u3) + co2*(u2+t3)
      c(jg+j)   = co6*(t2+u3) - si6*(u2-t3)
      c(jg+j+1) = si6*(t2+u3) + co6*(u2-t3)
      ajb = c(jb+j)
      ajd = c(jd+j)
      t0 = ajb + ajd
      t2 = ajb - ajd
      ajf = c(jf+j)
      ajh = c(jh+j)
      t1 = ajf - ajh
      t3 = ajf + ajh
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u0 = bjb + bjd
      u2 = bjb - bjd
      bjf = c(jf+j+1)
      bjh = c(jh+j+1)
      u1 = bjf - bjh
      u3 = bjf + bjh
      c(jb+j)   = co1*(t0-u3) - si1*(u0+t3)
      c(jb+j+1) = si1*(t0-u3) + co1*(u0+t3)
      c(jh+j)   = co7*(t0+u3) - si7*(u0-t3)
      c(jh+j+1) = si7*(t0+u3) + co7*(u0-t3)
      c(jd+j)   = co3*(t2+u1) - si3*(u2-t1)
      c(jd+j+1) = si3*(t2+u1) + co3*(u2-t1)
      c(jf+j)   = co5*(t2-u1) - si5*(u2+t1)
      c(jf+j+1) = si5*(t2-u1) + co5*(u2+t1)
      j = j + jump
  340 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  342 continue
  345 continue
      kk = kk + 2 * la
  350 continue
!
      la = 8 * la
!
!  loop on type II radix-4 passes
!  ------------------------------
  400 continue
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0_DP
      if (mu.eq.3) ss = -1.0_DP
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 4*jstep
!
      do 420 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 410 l = 1 , nvex
      aja = c(ja+j)
      ajc = c(jc+j)
      t0 = aja + ajc
      t2 = aja - ajc
      ajb = c(jb+j)
      ajd = c(jd+j)
      t1 = ajb + ajd
      t3 = ss * ( ajb - ajd )
      aji = c(ji+j)
      ajc =  aji
      bja = c(ja+j+1)
      bjc = c(jc+j+1)
      u0 = bja + bjc
      u2 = bja - bjc
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u1 = bjb + bjd
      u3 = ss * ( bjb - bjd )
      aje = c(je+j)
      ajb =  aje
      c(ja+j)   = t0 + t1
      c(ji+j)   = t0 - t1
      c(ja+j+1) = u0 + u1
      bjc =  u0 - u1
      bjm = c(jm+j+1)
      bjd =  bjm
      c(je+j)   = t2 - u3
      ajd =  t2 + u3
      bjb =  u2 + t3
      c(jm+j+1) = u2 - t3
!----------------------
      ajg = c(jg+j)
      t0 = ajb + ajg
      t2 = ajb - ajg
      ajf = c(jf+j)
      ajh = c(jh+j)
      t1 = ajf + ajh
      t3 = ss * ( ajf - ajh )
      ajj = c(jj+j)
      ajg =  ajj
      bje = c(je+j+1)
      bjg = c(jg+j+1)
      u0 = bje + bjg
      u2 = bje - bjg
      bjf = c(jf+j+1)
      bjh = c(jh+j+1)
      u1 = bjf + bjh
      u3 = ss * ( bjf - bjh )
      c(je+j+1) = bjb
      c(jb+j)   = t0 + t1
      c(jj+j)   = t0 - t1
      bjj = c(jj+j+1)
      bjg =  bjj
      c(jb+j+1) = u0 + u1
      c(jj+j+1) = u0 - u1
      c(jf+j)   = t2 - u3
      ajh =  t2 + u3
      c(jf+j+1) = u2 + t3
      bjh =  u2 - t3
!----------------------
      ajk = c(jk+j)
      t0 = ajc + ajk
      t2 = ajc - ajk
      ajl = c(jl+j)
      t1 = ajg + ajl
      t3 = ss * ( ajg - ajl )
      bji = c(ji+j+1)
      bjk = c(jk+j+1)
      u0 = bji + bjk
      u2 = bji - bjk
      ajo = c(jo+j)
      ajl =  ajo
      bjl = c(jl+j+1)
      u1 = bjg + bjl
      u3 = ss * ( bjg - bjl )
      c(ji+j+1) = bjc
      c(jc+j)   = t0 + t1
      c(jk+j)   = t0 - t1
      bjo = c(jo+j+1)
      bjl =  bjo
      c(jc+j+1) = u0 + u1
      c(jk+j+1) = u0 - u1
      c(jg+j)   = t2 - u3
      c(jo+j)   = t2 + u3
      c(jg+j+1) = u2 + t3
      c(jo+j+1) = u2 - t3
!----------------------
      ajm = c(jm+j)
      t0 = ajm + ajl
      t2 = ajm - ajl
      ajn = c(jn+j)
      ajp = c(jp+j)
      t1 = ajn + ajp
      t3 = ss * ( ajn - ajp )
      c(jm+j) = ajd
      u0 = bjd + bjl
      u2 = bjd - bjl
      bjn = c(jn+j+1)
      bjp = c(jp+j+1)
      u1 = bjn + bjp
      u3 = ss * ( bjn - bjp )
      c(jn+j)   = ajh
      c(jd+j)   = t0 + t1
      c(jl+j)   = t0 - t1
      c(jd+j+1) = u0 + u1
      c(jl+j+1) = u0 - u1
      c(jn+j+1) = bjh
      c(jh+j)   = t2 - u3
      c(jp+j)   = t2 + u3
      c(jh+j+1) = u2 + t3
      c(jp+j+1) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 4*jstep
!
      do 450 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 440 l = 1 , nvex
      aja = c(ja+j)
      ajc = c(jc+j)
      t0 = aja + ajc
      t2 = aja - ajc
      ajb = c(jb+j)
      ajd = c(jd+j)
      t1 = ajb + ajd
      t3 = ss * ( ajb - ajd )
      aji = c(ji+j)
      ajc =  aji
      bja = c(ja+j+1)
      bjc = c(jc+j+1)
      u0 = bja + bjc
      u2 = bja - bjc
      bjb = c(jb+j+1)
      bjd = c(jd+j+1)
      u1 = bjb + bjd
      u3 = ss * ( bjb - bjd )
      aje = c(je+j)
      ajb =  aje
      c(ja+j)   = t0 + t1
      c(ja+j+1) = u0 + u1
      c(je+j)   = co1*(t2-u3) - si1*(u2+t3)
      bjb =  si1*(t2-u3) + co1*(u2+t3)
      bjm = c(jm+j+1)
      bjd =  bjm
      c(ji+j)   = co2*(t0-t1) - si2*(u0-u1)
      bjc =  si2*(t0-t1) + co2*(u0-u1)
      ajd =  co3*(t2+u3) - si3*(u2-t3)
      c(jm+j+1) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      ajg = c(jg+j)
      t0 = ajb + ajg
      t2 = ajb - ajg
      ajf = c(jf+j)
      ajh = c(jh+j)
      t1 = ajf + ajh
      t3 = ss * ( ajf - ajh )
      ajj = c(jj+j)
      ajg =  ajj
      bje = c(je+j+1)
      bjg = c(jg+j+1)
      u0 = bje + bjg
      u2 = bje - bjg
      bjf = c(jf+j+1)
      bjh = c(jh+j+1)
      u1 = bjf + bjh
      u3 = ss * ( bjf - bjh )
      c(je+j+1) = bjb
      c(jb+j)   = t0 + t1
      c(jb+j+1) = u0 + u1
      bjj = c(jj+j+1)
      bjg =  bjj
      c(jf+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jf+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jj+j)   = co2*(t0-t1) - si2*(u0-u1)
      c(jj+j+1) = si2*(t0-t1) + co2*(u0-u1)
      ajh =  co3*(t2+u3) - si3*(u2-t3)
      bjh =  si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      ajk = c(jk+j)
      t0 = ajc + ajk
      t2 = ajc - ajk
      ajl = c(jl+j)
      t1 = ajg + ajl
      t3 = ss * ( ajg - ajl )
      bji = c(ji+j+1)
      bjk = c(jk+j+1)
      u0 = bji + bjk
      u2 = bji - bjk
      ajo = c(jo+j)
      ajl =  ajo
      bjl = c(jl+j+1)
      u1 = bjg + bjl
      u3 = ss * ( bjg - bjl )
      c(ji+j+1) = bjc
      c(jc+j)   = t0 + t1
      c(jc+j+1) = u0 + u1
      bjo = c(jo+j+1)
      bjl =  bjo
      c(jg+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jg+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jk+j)   = co2*(t0-t1) - si2*(u0-u1)
      c(jk+j+1) = si2*(t0-t1) + co2*(u0-u1)
      c(jo+j)   = co3*(t2+u3) - si3*(u2-t3)
      c(jo+j+1) = si3*(t2+u3) + co3*(u2-t3)
!----------------------------------------
      ajm = c(jm+j)
      t0 = ajm + ajl
      t2 = ajm - ajl
      ajn = c(jn+j)
      ajp = c(jp+j)
      t1 = ajn + ajp
      t3 = ss * ( ajn - ajp )
      c(jm+j) = ajd
      u0 = bjd + bjl
      u2 = bjd - bjl
      c(jn+j) = ajh
      bjn = c(jn+j+1)
      bjp = c(jp+j+1)
      u1 = bjn + bjp
      u3 = ss * ( bjn - bjp )
      c(jn+j+1) = bjh
      c(jd+j)   = t0 + t1
      c(jd+j+1) = u0 + u1
      c(jh+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jh+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jl+j)   = co2*(t0-t1) - si2*(u0-u1)
      c(jl+j+1) = si2*(t0-t1) + co2*(u0-u1)
      c(jp+j)   = co3*(t2+u3) - si3*(u2-t3)
      c(jp+j+1) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 4*la
  480 continue
!-----( end of loop on type II radix-4 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa2f

      subroutine gpfa3f(c,trigs,inc,jump,n,mm,lot,isign)
!     subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)
!
!     fortran version of *gpfa3* -
!     radix-3 section of self-sorting, in-place
!        generalized PFA
!
!-------------------------------------------------------------------
!
      use mod_comp

      implicit none

      real(DP), intent(inout) :: c(*)
      real(DP), intent(in)    :: trigs(*)
      integer,       intent(in)    :: inc,isign,jump,lot,mm,n

      real(DP) :: aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji,bja
      real(DP) :: bjb,bjc,bjd,bje,bjf,bjg,bjh,bji,c1,co1,co2
      real(DP) :: s,si1,si2,t1,t2,t3,u1,u2,u3
      integer :: ink,inq,ipass,istart,j,ja,jb,jc,jd,je,jf,jg,jh,ji
      integer :: jjj,jstep,jstepl,jstepx,k,kk,l,la,laincl,left,ll
      integer :: m,mh,mu,n3,nb,nblox,ninc,nu,nvex

      integer,       parameter :: lvr=64
      real(DP), parameter :: sin60=0.866025403784438647_DP
!
!     ***************************************************************
!     *                                                             *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *                                                             *
!     ***************************************************************
!
      n3 = 3**mm
      inq = n/n3
      jstepx = (n3-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,3)
      if (isign.eq.-1) mu = 3-mu
      m = mm
      mh = (m+1)/2
      s = real(isign,DP)
      c1 = sin60
      if (mu.eq.2) c1 = -c1
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = real(isign,DP)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type I radix-3 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 110 l = 1 , nvex
      ajb = c(jb+j)
      ajc = c(jc+j)
      t1 = ajb + ajc
      aja = c(ja+j)
      t2 = aja - 0.5_DP * t1
      t3 = c1 * ( ajb - ajc )
      bjb = c(jb+j+1)
      bjc = c(jc+j+1)
      u1 = bjb + bjc
      bja = c(ja+j+1)
      u2 = bja - 0.5_DP * u1
      u3 = c1 * ( bjb - bjc )
      c(ja+j)   = aja + t1
      c(ja+j+1) = bja + u1
      c(jb+j)   = t2 - u3
      c(jb+j+1) = u2 + t3
      c(jc+j)   = t2 + u3
      c(jc+j+1) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
!
!  finished if n3 = 3
!  ------------------
      if (n3.eq.3) go to 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 3*jstep
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 130 l = 1 , nvex
      ajb = c(jb+j)
      ajc = c(jc+j)
      t1 = ajb + ajc
      aja = c(ja+j)
      t2 = aja - 0.5_DP * t1
      t3 = c1 * ( ajb - ajc )
      bjb = c(jb+j+1)
      bjc = c(jc+j+1)
      u1 = bjb + bjc
      bja = c(ja+j+1)
      u2 = bja - 0.5_DP * u1
      u3 = c1 * ( bjb - bjc )
      c(ja+j)   = aja + t1
      c(ja+j+1) = bja + u1
      c(jb+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jb+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jc+j)   = co2*(t2+u3) - si2*(u2-t3)
      c(jc+j+1) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  130 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 3*la
  160 continue
!-----( end of loop on type I radix-3 passes)
!
!  loop on type II radix-3 passes
!  ------------------------------
  400 continue
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
      laincl = la*ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 3*jstep
!
      do 420 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 410 l = 1 , nvex
      ajb = c(jb+j)
      ajc = c(jc+j)
      t1 = ajb + ajc
      aja = c(ja+j)
      t2 = aja - 0.5_DP * t1
      t3 = c1 * ( ajb - ajc )
      ajd = c(jd+j)
      ajb =  ajd
      bjb = c(jb+j+1)
      bjc = c(jc+j+1)
      u1 = bjb + bjc
      bja = c(ja+j+1)
      u2 = bja - 0.5_DP * u1
      u3 = c1 * ( bjb - bjc )
      bjd = c(jd+j+1)
      bjb =  bjd
      c(ja+j)   = aja + t1
      c(ja+j+1) = bja + u1
      c(jd+j)   = t2 - u3
      c(jd+j+1) = u2 + t3
      ajc =  t2 + u3
      bjc =  u2 - t3
!----------------------
      aje = c(je+j)
      ajf = c(jf+j)
      t1 = aje + ajf
      t2 = ajb - 0.5_DP * t1
      t3 = c1 * ( aje - ajf )
      ajh = c(jh+j)
      ajf =  ajh
      bje = c(je+j+1)
      bjf = c(jf+j+1)
      u1 = bje + bjf
      u2 = bjb - 0.5_DP * u1
      u3 = c1 * ( bje - bjf )
      bjh = c(jh+j+1)
      bjf =  bjh
      c(jb+j)   = ajb + t1
      c(jb+j+1) = bjb + u1
      c(je+j)   = t2 - u3
      c(je+j+1) = u2 + t3
      c(jh+j)   = t2 + u3
      c(jh+j+1) = u2 - t3
!----------------------
      aji = c(ji+j)
      t1 = ajf + aji
      ajg = c(jg+j)
      t2 = ajg - 0.5_DP * t1
      t3 = c1 * ( ajf - aji )
      t1 = ajg + t1
      c(jg+j)   = ajc
      bji = c(ji+j+1)
      u1 = bjf + bji
      bjg = c(jg+j+1)
      u2 = bjg - 0.5_DP * u1
      u3 = c1 * ( bjf - bji )
      u1 = bjg + u1
      c(jg+j+1) = bjc
      c(jc+j)   = t1
      c(jc+j+1) = u1
      c(jf+j)   = t2 - u3
      c(jf+j+1) = u2 + t3
      c(ji+j)   = t2 + u3
      c(ji+j+1) = u2 - t3
      j = j + jump
  410 continue
!-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass.eq.m) go to 490
!
      kk = 2*la
!
!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 3*jstep
!
      do 450 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
!
!  "transverse" loop
!  -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      do 440 l = 1 , nvex
      ajb = c(jb+j)
      ajc = c(jc+j)
      t1 = ajb + ajc
      aja = c(ja+j)
      t2 = aja - 0.5_DP * t1
      t3 = c1 * ( ajb - ajc )
      ajd = c(jd+j)
      ajb =  ajd
      bjb = c(jb+j+1)
      bjc = c(jc+j+1)
      u1 = bjb + bjc
      bja = c(ja+j+1)
      u2 = bja - 0.5_DP * u1
      u3 = c1 * ( bjb - bjc )
      bjd = c(jd+j+1)
      bjb =  bjd
      c(ja+j)   = aja + t1
      c(ja+j+1) = bja + u1
      c(jd+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jd+j+1) = si1*(t2-u3) + co1*(u2+t3)
      ajc =  co2*(t2+u3) - si2*(u2-t3)
      bjc =  si2*(t2+u3) + co2*(u2-t3)
!----------------------
      aje = c(je+j)
      ajf = c(jf+j)
      t1 = aje + ajf
      t2 = ajb - 0.5_DP * t1
      t3 = c1 * ( aje - ajf )
      ajh = c(jh+j)
      ajf =  ajh
      bje = c(je+j+1)
      bjf = c(jf+j+1)
      u1 = bje + bjf
      u2 = bjb - 0.5_DP * u1
      u3 = c1 * ( bje - bjf )
      bjh = c(jh+j+1)
      bjf =  bjh
      c(jb+j)   = ajb + t1
      c(jb+j+1) = bjb + u1
      c(je+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(je+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(jh+j)   = co2*(t2+u3) - si2*(u2-t3)
      c(jh+j+1) = si2*(t2+u3) + co2*(u2-t3)
!----------------------
      aji = c(ji+j)
      t1 = ajf + aji
      ajg = c(jg+j)
      t2 = ajg - 0.5_DP * t1
      t3 = c1 * ( ajf - aji )
      t1 = ajg + t1
      c(jg+j)   = ajc
      bji = c(ji+j+1)
      u1 = bjf + bji
      bjg = c(jg+j+1)
      u2 = bjg - 0.5_DP * u1
      u3 = c1 * ( bjf - bji )
      u1 = bjg + u1
      c(jg+j+1) = bjc
      c(jc+j)   = t1
      c(jc+j+1) = u1
      c(jf+j)   = co1*(t2-u3) - si1*(u2+t3)
      c(jf+j+1) = si1*(t2-u3) + co1*(u2+t3)
      c(ji+j)   = co2*(t2+u3) - si2*(u2-t3)
      c(ji+j+1) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  440 continue
!-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 3*la
  480 continue
!-----( end of loop on type II radix-3 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa3f


      subroutine gpfa5f(c,trigs,inc,jump,n,mm,lot,isign)
!     subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)
!     fortran version of *gpfa5* -
!     radix-5 section of self-sorting, in-place,
!        generalized pfa
!
!-------------------------------------------------------------------
!
      use mod_comp

      implicit none

      real(DP), intent(inout) :: c(*)
      real(DP), intent(in)    :: trigs(*)
      integer,       intent(in)    :: inc,isign,jump,lot,mm,n

      real(DP) :: aja,ajb,ajc,ajd,aje,ajf,ajg,ajh,aji,ajj
      real(DP) :: ajk,ajl,ajm,ajn,ajo,ajp,ajq,ajr,ajs,ajt
      real(DP) :: aju,ajv,ajw,ajx,ajy,ax
      real(DP) :: bja,bjb,bjc,bjd,bje,bjf,bjg,bjh,bji,bjj
      real(DP) :: bjk,bjl,bjm,bjn,bjo,bjp,bjq,bjr,bjs,bjt
      real(DP) :: bju,bjv,bjw,bjx,bjy,bx
      real(DP) :: c1,c2,c3,co1,co2,co3,co4
      real(DP) :: s,si1,si2,si3,si4
      real(DP) :: t1,t10,t11,t2,t3,t4,t5,t6,t7,t8,t9
      real(DP) :: u1,u10,u11,u2,u3,u4,u5,u6,u7,u8,u9
      integer :: ink,inq,ipass,istart,j,ja,jb,jc,jd,je,jf,jg,jh
      integer :: ji,jj,jjj,jk,jl,jm,jn,jo,jp,jq,jr,js
      integer :: jstep,jstepl,jstepx,jt,ju,jv,jw,jx,jy,k,kk,l,la
      integer :: laincl,left,ll,m,mh,mu,n5,nb,nblox,ninc,nu,nvex

      real(DP), parameter :: sin36=0.587785252292473129_DP
      real(DP), parameter :: sin72=0.951056516295153572_DP
      real(DP), parameter ::  qrt5=0.559016994374947424_DP
      integer,       parameter :: lvr=64
!
!     ***************************************************************
!     *                                                             *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *                                                             *
!     ***************************************************************
!
      n5 = 5 ** mm
      inq = n / n5
      jstepx = (n5-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,5)
      if (isign.eq.-1) mu = 5 - mu
!
      m = mm
      mh = (m+1)/2
      s = real(isign,DP)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu.eq.2.or.mu.eq.3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      endif
      if (mu.eq.3.or.mu.eq.4) c2 = -c2
      if (mu.eq.2.or.mu.eq.4) c3 = -c3
!
      nblox = 1 + (lot-1)/lvr
      left = lot
      s = real(isign,DP)
      istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
      do 500 nb = 1 , nblox
!
      if (left.le.lvr) then
         nvex = left
      else if (left.lt.(2*lvr)) then
         nvex = left/2
         nvex = nvex + mod(nvex,2)
      else
         nvex = lvr
      endif
      left = left - nvex
!
      la = 1
!
!  loop on type I radix-5 passes
!  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      kk = 0
!
!  loop on k
!  ---------
      do 150 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  loop along transform
!  --------------------
      do 140 jjj = k , (n-1)*inc , 5*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
      do 110 l = 1 , nvex
      ajb = c(jb+j)
      aje = c(je+j)
      t1 = ajb + aje
      ajc = c(jc+j)
      ajd = c(jd+j)
      t2 = ajc + ajd
      t3 = ajb - aje
      t4 = ajc - ajd
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aja = c(ja+j)
      t7 = aja - 0.25_DP * t5
      c(ja+j) = aja + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjb = c(jb+j+1)
      bje = c(je+j+1)
      u1 = bjb + bje
      bjc = c(jc+j+1)
      bjd = c(jd+j+1)
      u2 = bjc + bjd
      u3 = bjb - bje
      u4 = bjc - bjd
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bja = c(ja+j+1)
      u7 = bja - 0.25_DP * u5
      c(ja+j+1) = bja + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jb+j)   = t8 - u11
      c(jb+j+1) = u8 + t11
      c(je+j)   = t8 + u11
      c(je+j+1) = u8 - t11
      c(jc+j)   = t9 - u10
      c(jc+j+1) = u9 + t10
      c(jd+j)   = t9 + u10
      c(jd+j+1) = u9 - t10
      j = j + jump
  110 continue
!
      else
!
      do 130 l = 1 , nvex
      ajb = c(jb+j)
      aje = c(je+j)
      t1 = ajb + aje
      ajc = c(jc+j)
      ajd = c(jd+j)
      t2 = ajc + ajd
      t3 = ajb - aje
      t4 = ajc - ajd
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aja = c(ja+j)
      t7 = aja - 0.25_DP * t5
      c(ja+j) = aja + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjb = c(jb+j+1)
      bje = c(je+j+1)
      u1 = bjb + bje
      bjc = c(jc+j+1)
      bjd = c(jd+j+1)
      u2 = bjc + bjd
      u3 = bjb - bje
      u4 = bjc - bjd
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bja = c(ja+j+1)
      u7 = bja - 0.25_DP * u5
      c(ja+j+1) = bja + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jb+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(jb+j+1) = si1*(t8-u11) + co1*(u8+t11)
      c(je+j)   = co4*(t8+u11) - si4*(u8-t11)
      c(je+j+1) = si4*(t8+u11) + co4*(u8-t11)
      c(jc+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jc+j+1) = si2*(t9-u10) + co2*(u9+t10)
      c(jd+j)   = co3*(t9+u10) - si3*(u9-t10)
      c(jd+j+1) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  130 continue
!
      endif
!
!-----( end of loop across transforms )
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 5*la
  160 continue
!-----( end of loop on type I radix-5 passes)
!
      if (n.eq.5) go to 490
!
!  loop on type II radix-5 passes
!  ------------------------------
  400 continue
!
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
      kk = 0
!
!     loop on k
!     ---------
      do 470 k = 0 , jstep-ink , ink
!
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 5*jstep
!
      do 450 jjj = ll , (n-1)*inc , 5*la*ink
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = ja + laincl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jf + laincl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = jl + jstepl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jk + laincl
      if (jp.lt.istart) jp = jp + ninc
      jq = jp + jstepl
      if (jq.lt.istart) jq = jq + ninc
      jr = jq + jstepl
      if (jr.lt.istart) jr = jr + ninc
      js = jr + jstepl
      if (js.lt.istart) js = js + ninc
      jt = js + jstepl
      if (jt.lt.istart) jt = jt + ninc
      ju = jp + laincl
      if (ju.lt.istart) ju = ju + ninc
      jv = ju + jstepl
      if (jv.lt.istart) jv = jv + ninc
      jw = jv + jstepl
      if (jw.lt.istart) jw = jw + ninc
      jx = jw + jstepl
      if (jx.lt.istart) jx = jx + ninc
      jy = jx + jstepl
      if (jy.lt.istart) jy = jy + ninc
      j = 0
!
!  loop across transforms
!  ----------------------
      if (k.eq.0) then
!
      do 410 l = 1 , nvex
      ajb = c(jb+j)
      aje = c(je+j)
      t1 = ajb + aje
      ajc = c(jc+j)
      ajd = c(jd+j)
      t2 = ajc + ajd
      t3 = ajb - aje
      t4 = ajc - ajd
      ajf = c(jf+j)
      ajb =  ajf
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aja = c(ja+j)
      t7 = aja - 0.25_DP * t5
      c(ja+j) = aja + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajk = c(jk+j)
      ajc =  ajk
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjb = c(jb+j+1)
      bje = c(je+j+1)
      u1 = bjb + bje
      bjc = c(jc+j+1)
      bjd = c(jd+j+1)
      u2 = bjc + bjd
      u3 = bjb - bje
      u4 = bjc - bjd
      bjf = c(jf+j+1)
      bjb =  bjf
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bja = c(ja+j+1)
      u7 = bja - 0.25_DP * u5
      c(ja+j+1) = bja + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjk = c(jk+j+1)
      bjc =  bjk
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jf+j)   = t8 - u11
      c(jf+j+1) = u8 + t11
      aje =  t8 + u11
      bje =  u8 - t11
      c(jk+j)   = t9 - u10
      c(jk+j+1) = u9 + t10
      ajd =  t9 + u10
      bjd =  u9 - t10
!----------------------
      ajg = c(jg+j)
      ajj = c(jj+j)
      t1 = ajg + ajj
      ajh = c(jh+j)
      aji = c(ji+j)
      t2 = ajh + aji
      t3 = ajg - ajj
      t4 = ajh - aji
      ajl = c(jl+j)
      ajh =  ajl
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = ajb - 0.25_DP * t5
      c(jb+j) = ajb + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajq = c(jq+j)
      aji =  ajq
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjg = c(jg+j+1)
      bjj = c(jj+j+1)
      u1 = bjg + bjj
      bjh = c(jh+j+1)
      bji = c(ji+j+1)
      u2 = bjh + bji
      u3 = bjg - bjj
      u4 = bjh - bji
      bjl = c(jl+j+1)
      bjh =  bjl
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = bjb - 0.25_DP * u5
      c(jb+j+1) = bjb + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjq = c(jq+j+1)
      bji =  bjq
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jg+j)   = t8 - u11
      c(jg+j+1) = u8 + t11
      ajj =  t8 + u11
      bjj =  u8 - t11
      c(jl+j)   = t9 - u10
      c(jl+j+1) = u9 + t10
      c(jq+j)   = t9 + u10
      c(jq+j+1) = u9 - t10
!----------------------
      ajo = c(jo+j)
      t1 = ajh + ajo
      ajm = c(jm+j)
      ajn = c(jn+j)
      t2 = ajm + ajn
      t3 = ajh - ajo
      t4 = ajm - ajn
      ajr = c(jr+j)
      ajn =  ajr
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = ajc - 0.25_DP * t5
      c(jc+j) = ajc + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajw = c(jw+j)
      ajo =  ajw
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjo = c(jo+j+1)
      u1 = bjh + bjo
      bjm = c(jm+j+1)
      bjn = c(jn+j+1)
      u2 = bjm + bjn
      u3 = bjh - bjo
      u4 = bjm - bjn
      bjr = c(jr+j+1)
      bjn =  bjr
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = bjc - 0.25_DP * u5
      c(jc+j+1) = bjc + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjw = c(jw+j+1)
      bjo =  bjw
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jh+j)   = t8 - u11
      c(jh+j+1) = u8 + t11
      c(jw+j)   = t8 + u11
      c(jw+j+1) = u8 - t11
      c(jm+j)   = t9 - u10
      c(jm+j+1) = u9 + t10
      c(jr+j)   = t9 + u10
      c(jr+j+1) = u9 - t10
!----------------------
      ajt = c(jt+j)
      t1 = aji + ajt
      ajs = c(js+j)
      t2 = ajn + ajs
      t3 = aji - ajt
      t4 = ajn - ajs
      ajx = c(jx+j)
      ajt =  ajx
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      ajp = c(jp+j)
      t7 = ajp - 0.25_DP * t5
      ax = ajp + t5
      t8 = t7 + t6
      t9 = t7 - t6
      c(jp+j) = ajd
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      c(jd+j) = ax
      bjt = c(jt+j+1)
      u1 = bji + bjt
      bjs = c(js+j+1)
      u2 = bjn + bjs
      u3 = bji - bjt
      u4 = bjn - bjs
      bjx = c(jx+j+1)
      bjt =  bjx
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bjp = c(jp+j+1)
      u7 = bjp - 0.25_DP * u5
      bx = bjp + u5
      u8 = u7 + u6
      u9 = u7 - u6
      c(jp+j+1) = bjd
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jd+j+1) = bx
      c(ji+j)   = t8 - u11
      c(ji+j+1) = u8 + t11
      c(jx+j)   = t8 + u11
      c(jx+j+1) = u8 - t11
      c(jn+j)   = t9 - u10
      c(jn+j+1) = u9 + t10
      c(js+j)   = t9 + u10
      c(js+j+1) = u9 - t10
!----------------------
      ajv = c(jv+j)
      ajy = c(jy+j)
      t1 = ajv + ajy
      t2 = ajo + ajt
      t3 = ajv - ajy
      t4 = ajo - ajt
      c(jv+j) = ajj
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aju = c(ju+j)
      t7 = aju - 0.25_DP * t5
      ax = aju + t5
      t8 = t7 + t6
      t9 = t7 - t6
      c(ju+j) = aje
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      c(je+j) = ax
      bjv = c(jv+j+1)
      bjy = c(jy+j+1)
      u1 = bjv + bjy
      u2 = bjo + bjt
      u3 = bjv - bjy
      u4 = bjo - bjt
      c(jv+j+1) = bjj
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bju = c(ju+j+1)
      u7 = bju - 0.25_DP * u5
      bx = bju + u5
      u8 = u7 + u6
      u9 = u7 - u6
      c(ju+j+1) = bje
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(je+j+1) = bx
      c(jj+j)   = t8 - u11
      c(jj+j+1) = u8 + t11
      c(jy+j)   = t8 + u11
      c(jy+j+1) = u8 - t11
      c(jo+j)   = t9 - u10
      c(jo+j+1) = u9 + t10
      c(jt+j)   = t9 + u10
      c(jt+j+1) = u9 - t10
      j = j + jump
  410 continue
!
      else
!
      do 440 l = 1 , nvex
      ajb = c(jb+j)
      aje = c(je+j)
      t1 = ajb + aje
      ajc = c(jc+j)
      ajd = c(jd+j)
      t2 = ajc + ajd
      t3 = ajb - aje
      t4 = ajc - ajd
      ajf = c(jf+j)
      ajb =  ajf
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aja = c(ja+j)
      t7 = aja - 0.25_DP * t5
      c(ja+j) = aja + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajk = c(jk+j)
      ajc =  ajk
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjb = c(jb+j+1)
      bje = c(je+j+1)
      u1 = bjb + bje
      bjc = c(jc+j+1)
      bjd = c(jd+j+1)
      u2 = bjc + bjd
      u3 = bjb - bje
      u4 = bjc - bjd
      bjf = c(jf+j+1)
      bjb =  bjf
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bja = c(ja+j+1)
      u7 = bja - 0.25_DP * u5
      c(ja+j+1) = bja + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjk = c(jk+j+1)
      bjc =  bjk
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jf+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(jf+j+1) = si1*(t8-u11) + co1*(u8+t11)
      aje =  co4*(t8+u11) - si4*(u8-t11)
      bje =  si4*(t8+u11) + co4*(u8-t11)
      c(jk+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jk+j+1) = si2*(t9-u10) + co2*(u9+t10)
      ajd =  co3*(t9+u10) - si3*(u9-t10)
      bjd =  si3*(t9+u10) + co3*(u9-t10)
!----------------------
      ajg = c(jg+j)
      ajj = c(jj+j)
      t1 = ajg + ajj
      ajh = c(jh+j)
      aji = c(ji+j)
      t2 = ajh + aji
      t3 = ajg - ajj
      t4 = ajh - aji
      ajl = c(jl+j)
      ajh =  ajl
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = ajb - 0.25_DP * t5
      c(jb+j) = ajb + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajq = c(jq+j)
      aji =  ajq
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjg = c(jg+j+1)
      bjj = c(jj+j+1)
      u1 = bjg + bjj
      bjh = c(jh+j+1)
      bji = c(ji+j+1)
      u2 = bjh + bji
      u3 = bjg - bjj
      u4 = bjh - bji
      bjl = c(jl+j+1)
      bjh =  bjl
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = bjb - 0.25_DP * u5
      c(jb+j+1) = bjb + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjq = c(jq+j+1)
      bji =  bjq
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jg+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(jg+j+1) = si1*(t8-u11) + co1*(u8+t11)
      ajj =  co4*(t8+u11) - si4*(u8-t11)
      bjj =  si4*(t8+u11) + co4*(u8-t11)
      c(jl+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jl+j+1) = si2*(t9-u10) + co2*(u9+t10)
      c(jq+j)   = co3*(t9+u10) - si3*(u9-t10)
      c(jq+j+1) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      ajo = c(jo+j)
      t1 = ajh + ajo
      ajm = c(jm+j)
      ajn = c(jn+j)
      t2 = ajm + ajn
      t3 = ajh - ajo
      t4 = ajm - ajn
      ajr = c(jr+j)
      ajn =  ajr
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = ajc - 0.25_DP * t5
      c(jc+j) = ajc + t5
      t8 = t7 + t6
      t9 = t7 - t6
      ajw = c(jw+j)
      ajo =  ajw
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      bjo = c(jo+j+1)
      u1 = bjh + bjo
      bjm = c(jm+j+1)
      bjn = c(jn+j+1)
      u2 = bjm + bjn
      u3 = bjh - bjo
      u4 = bjm - bjn
      bjr = c(jr+j+1)
      bjn =  bjr
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = bjc - 0.25_DP * u5
      c(jc+j+1) = bjc + u5
      u8 = u7 + u6
      u9 = u7 - u6
      bjw = c(jw+j+1)
      bjo =  bjw
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jh+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(jh+j+1) = si1*(t8-u11) + co1*(u8+t11)
      c(jw+j)   = co4*(t8+u11) - si4*(u8-t11)
      c(jw+j+1) = si4*(t8+u11) + co4*(u8-t11)
      c(jm+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jm+j+1) = si2*(t9-u10) + co2*(u9+t10)
      c(jr+j)   = co3*(t9+u10) - si3*(u9-t10)
      c(jr+j+1) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      ajt = c(jt+j)
      t1 = aji + ajt
      ajs = c(js+j)
      t2 = ajn + ajs
      t3 = aji - ajt
      t4 = ajn - ajs
      ajx = c(jx+j)
      ajt =  ajx
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      ajp = c(jp+j)
      t7 = ajp - 0.25_DP * t5
      ax = ajp + t5
      t8 = t7 + t6
      t9 = t7 - t6
      c(jp+j) = ajd
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      c(jd+j) = ax
      bjt = c(jt+j+1)
      u1 = bji + bjt
      bjs = c(js+j+1)
      u2 = bjn + bjs
      u3 = bji - bjt
      u4 = bjn - bjs
      bjx = c(jx+j+1)
      bjt =  bjx
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bjp = c(jp+j+1)
      u7 = bjp - 0.25_DP * u5
      bx = bjp + u5
      u8 = u7 + u6
      u9 = u7 - u6
      c(jp+j+1) = bjd
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(jd+j+1) = bx
      c(ji+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(ji+j+1) = si1*(t8-u11) + co1*(u8+t11)
      c(jx+j)   = co4*(t8+u11) - si4*(u8-t11)
      c(jx+j+1) = si4*(t8+u11) + co4*(u8-t11)
      c(jn+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jn+j+1) = si2*(t9-u10) + co2*(u9+t10)
      c(js+j)   = co3*(t9+u10) - si3*(u9-t10)
      c(js+j+1) = si3*(t9+u10) + co3*(u9-t10)
!----------------------
      ajv = c(jv+j)
      ajy = c(jy+j)
      t1 = ajv + ajy
      t2 = ajo + ajt
      t3 = ajv - ajy
      t4 = ajo - ajt
      c(jv+j) = ajj
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      aju = c(ju+j)
      t7 = aju - 0.25_DP * t5
      ax = aju + t5
      t8 = t7 + t6
      t9 = t7 - t6
      c(ju+j) = aje
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      c(je+j) = ax
      bjv = c(jv+j+1)
      bjy = c(jy+j+1)
      u1 = bjv + bjy
      u2 = bjo + bjt
      u3 = bjv - bjy
      u4 = bjo - bjt
      c(jv+j+1) = bjj
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      bju = c(ju+j+1)
      u7 = bju - 0.25_DP * u5
      bx = bju + u5
      u8 = u7 + u6
      u9 = u7 - u6
      c(ju+j+1) = bje
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      c(je+j+1) = bx
      c(jj+j)   = co1*(t8-u11) - si1*(u8+t11)
      c(jj+j+1) = si1*(t8-u11) + co1*(u8+t11)
      c(jy+j)   = co4*(t8+u11) - si4*(u8-t11)
      c(jy+j+1) = si4*(t8+u11) + co4*(u8-t11)
      c(jo+j)   = co2*(t9-u10) - si2*(u9+t10)
      c(jo+j+1) = si2*(t9-u10) + co2*(u9+t10)
      c(jt+j)   = co3*(t9+u10) - si3*(u9-t10)
      c(jt+j+1) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  440 continue
!
      endif
!
!-----(end of loop across transforms)
!
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 5*la
  480 continue
!-----( end of loop on type II radix-5 passes )
!-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
      return
      end subroutine gpfa5f

!====================================================================
