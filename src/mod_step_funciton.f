      Module mod_step_funciton
        use mod_comp
        use mod_input
        use mod_wave
        use mod_fft
        use mod_lattice
      
      contains

      subroutine step_zlayer (f,nfft_tmp,zlay1_tmp,zlay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between zlay1 and zlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: zlay1_tmp,zlay2_tmp

      integer      :: k3,k3m,nf3
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k3m = min( k3d, (nfft_tmp(3)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k3=-k3m,k3m
         nf3 = modulo( k3+nfft_tmp(3), nfft_tmp(3) ) + 1

!--->    layer part: k1=0, k2=0
         if( k3 == 0 ) then
            ff = (zlay2_tmp - zlay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k3*zlay1_tmp
            arg2 = TWOPI*k3*zlay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k3*TWOPI)
         endif

         f(nf3) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine

       subroutine step_ylayer (f,nfft_tmp,ylay1_tmp,ylay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between zlay1 and zlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: ylay1_tmp,ylay2_tmp

      integer      :: k2,k2m,nf2
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k2m = min( k2d, (nfft_tmp(2)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k2=-k2m,k2m
         nf2 = modulo( k2+nfft_tmp(2), nfft_tmp(2) ) + 1

!--->    layer part: k1=0, k2=0
         if( k2 == 0 ) then
            ff = (ylay2_tmp - ylay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k2*ylay1_tmp
            arg2 = TWOPI*k2*ylay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k2*TWOPI)
         endif

         f(nf2) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine

       subroutine step_xlayer (f,nfft_tmp,xlay1_tmp,xlay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between xlay1 and xlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: xlay1_tmp,xlay2_tmp

      integer      :: k1,k1m,nf1
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k1m = min( k1d, (nfft_tmp(1)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k1=-k1m,k1m
         nf1 = modulo( k1+nfft_tmp(1), nfft_tmp(1) ) + 1

!--->    layer part: k1=0, k2=0
         if( k1 == 0 ) then
            ff = (xlay2_tmp - xlay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k1*xlay1_tmp
            arg2 = TWOPI*k1*xlay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k1*TWOPI)
         endif

         f(nf1) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine

       subroutine invert_3x3(mat1,inv_mat1)
         !!! invert of matrix
         REAL(DP)        :: mat1(3,3),inv_mat1(3,3),det

         inv_mat1(1,1)=mat1(2,2)*mat1(3,3)-mat1(3,2)*mat1(2,3)
         inv_mat1(1,2)=mat1(3,2)*mat1(1,3)-mat1(1,2)*mat1(3,3)
         inv_mat1(1,3)=mat1(1,2)*mat1(2,3)-mat1(2,2)*mat1(1,3)
         inv_mat1(2,1)=mat1(2,3)*mat1(3,1)-mat1(3,3)*mat1(2,1)
         inv_mat1(2,2)=mat1(3,3)*mat1(1,1)-mat1(1,3)*mat1(3,1)
         inv_mat1(2,3)=mat1(1,3)*mat1(2,1)-mat1(2,3)*mat1(1,1)
         inv_mat1(3,1)=mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2)
         inv_mat1(3,2)=mat1(3,1)*mat1(1,2)-mat1(1,1)*mat1(3,2)
         inv_mat1(3,3)=mat1(1,1)*mat1(2,2)-mat1(2,1)*mat1(1,2)
         call determ(mat1,det)
         inv_mat1(1:3,1:3)=inv_mat1(1:3,1:3)/det
       end subroutine invert_3x3

       subroutine determ(a,b)
       ! returnerer determinaten af a 3x3 matrix
       REAL(DP),intent(in)  :: a(3,3)
       real(DP),intent(out) :: b
         b=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
     & a(1,3)*a(3,2)*a(2,1)-a(1,3)*a(2,2)*a(3,1)- &
     & a(2,3)*a(3,2)*a(1,1)-a(3,3)* &
         a(1,2)*a(2,1)
       END subroutine determ
        !
        !-------------------------------------------------------------
        ! reciprocal representation to cartessian representation
        !-------------------------------------------------------------
        !
        subroutine dire2car(A,B,re_cell)
          real(DP)     :: A(3), B(3), re_cell(3,3)
          B(1)=A(1)*re_cell(1,1)+A(2)*re_cell(1,2)+A(3)*re_cell(1,3)
          B(2)=A(1)*re_cell(2,1)+A(2)*re_cell(2,2)+A(3)*re_cell(2,3)
          B(3)=A(1)*re_cell(3,1)+A(2)*re_cell(3,2)+A(3)*re_cell(3,3)
        end subroutine dire2car

        SUBROUTINE LATTIC(A,B)
           real(DP)   ::  A(3,3),B(3,3),ANORM(3),BNORM(3),OMEGA

           CALL EXPRO(B(1,1),A(1,2),A(1,3))
           CALL EXPRO(B(1,2),A(1,3),A(1,1))
           CALL EXPRO(B(1,3),A(1,1),A(1,2))

           OMEGA =ABS(B(1,1)*A(1,1)+B(2,1)*A(2,1)+B(3,1)*A(3,1))

           do I=1,3
             do J=1,3
                B(I,J)=TWOPI*B(I,J)/OMEGA
             enddo
           enddo

         END SUBROUTINE LATTIC
         !
         !----------------- SUBROUTINE EXPRO   ---------------------------------
         ! EXPRO
         ! caclulates the x-product of two vectors
         !
         !----------------------------------------------------------------------
         !
         SUBROUTINE EXPRO(H,U1,U2)
           real(DP)    :: H(3),U1(3),U2(3)

           H(1)=U1(2)*U2(3)-U1(3)*U2(2)
           H(2)=U1(3)*U2(1)-U1(1)*U2(3)
           H(3)=U1(1)*U2(2)-U1(2)*U2(1)
         END SUBROUTINE EXPRO

      subroutine out_bands
      !
      use mod_input
      use mod_wave
      use mod_lattice

      implicit none

!              wavecar file variables
      real(DP),allocatable :: e(:),occ(:)
      real(DP)             :: vkpt(3),rnv
      integer              :: nv

!              k-projected information
      real(DP),allocatable :: bkpts(:,:)

!              for HSE calculation
      real(DP),allocatable :: proj_wt_HSE(:,:,:),bkpts_HSE(:,:)
      real(DP),allocatable :: CELTOT_HSE(:,:,:),FERTOT_HSE(:,:,:)
      real(DP),allocatable :: KPTVEC_HSE(:,:)
      integer              :: NK_EFF ! the number of effctive kpoints

      logical,allocatable  :: l_good_kpt(:)

      integer              :: irec
      integer              :: i,k,l,n,nb,jspin,ikpt,ios,kpt,nkpts0
      real(DP)             :: emin_ef,emax_ef,ef

!     file name and handle for band data
      integer,   parameter :: bskfh = 73
      character(len=30)    :: bskfn = 'bs_projected.dat'

      character(len=256)   :: line

      integer, allocatable :: segment_endpoint(:),epts(:)
      real(DP),allocatable :: x_kpt(:)
      real(DP)             :: dk(3),dk_old(3),x0
      integer              :: kpt_lines,nseg,k0,max_segments

        if(LHSE)then
        NK_EFF = NKPTS - NK_HSE
        allocate(proj_wt_HSE(NBANDS,NK_EFF,ISPIN))
        allocate(CELTOT_HSE(NBANDS,NK_EFF,ISPIN))
        allocate(FERTOT_HSE(NBANDS,NK_EFF,ISPIN))
        allocate(KPTVEC_HSE(3,NK_EFF))

        do jspin=1,ISPIN
          do ikpt = 1,NK_EFF
            proj_wt_HSE(1:NBANDS,ikpt,jspin) = &
     &      proj_wt(1:NBANDS,ikpt+NK_HSE,jspin)

            CELTOT_HSE(1:NBANDS,ikpt,jspin) =  &
     &      CELTOT(1:NBANDS,ikpt+NK_HSE,jspin)

            FERTOT_HSE(1:NBANDS,ikpt,jspin) =  &
     &      FERTOT(1:NBANDS,ikpt+NK_HSE,jspin)

            KPTVEC_HSE(1:3,ikpt) = KPTVEC(1:3,ikpt+NK_HSE)
          enddo
        enddo

        NKPTS = NK_EFF

        do jspin=1,ISPIN
          do ikpt = 1,NK_EFF
            proj_wt(1:NBANDS,ikpt,jspin) = &
     &      proj_wt_HSE(1:NBANDS,ikpt,jspin)

            CELTOT(1:NBANDS,ikpt,jspin) =  &
     &      CELTOT_HSE(1:NBANDS,ikpt,jspin)

            FERTOT(1:NBANDS,ikpt,jspin) =  &
     &      FERTOT_HSE(1:NBANDS,ikpt,jspin)

            KPTVEC(1:3,ikpt) = KPTVEC_HSE(1:3,ikpt)
          enddo
        enddo
      endif

            allocate( e(NBANDS),occ(NBANDS))
      emin_ef = -1000
      emax_ef =  1000
      do jspin=1,ISPIN
        do ikpt=1,NKPTS
           e(1:NBANDS)   = CELTOT(1:NBANDS,ikpt,jspin)
           occ(1:NBANDS) = FERTOT(1:NBANDS,ikpt,jspin)
!  get a better estimate of the fermi energy
         occ(1:nbands) = occ(1:nbands)/occ(1) ! assumes 1 occupied state
         emin_ef = max( maxval( e, mask = ( occ >= 0.5_DP ) ) , emin_ef )
         emax_ef = min( minval( e, mask = ( occ <  0.5_DP ) ) , emax_ef )
        enddo
      enddo
      ef = ( emin_ef + emax_ef )/2
      ef = 0.0_DP
      write(IUO,'(/,5X,"Fermi level obtained from bands = ",f0.5," eV")') ef
      write(IUO,'("    (need to change in bs_projected.dat)")')


!--->    break up the k-points into lines, assuming that the order is correct
      allocate( l_good_kpt(NKPTS) )
      allocate( bkpts(3,NKPTS) )
      bkpts = KPTVEC
      l_good_kpt = .true.
      max_segments = 8
      allocate( segment_endpoint(max_segments) )
      nseg = 0
      dk_old = bkpts(:,2) - bkpts(:,1)
      dk_old = dk_old/maxval( abs( dk_old ) ) ! normalize the direction
      do ikpt = 2,NKPTS-1
         dk = bkpts(:,ikpt+1) - bkpts(:,ikpt)
         if( all( abs( dk ) < 1.e-05_DP ) ) then ! vasp sometimes does
            l_good_kpt(ikpt) = .false.           ! k-points multiple times
            cycle
         endif
         dk = dk/maxval( abs( dk ) )
         if( any( abs( dk - dk_old ) > 1.e-04_DP ) ) then
            nseg = nseg + 1
            if( nseg > max_segments ) then
               allocate( epts( 2* nseg ) )
               epts = 0
               epts(1:nseg-1) = segment_endpoint(1:nseg-1)
               deallocate( segment_endpoint )
               allocate( segment_endpoint(2*nseg) )
               segment_endpoint = epts
               deallocate( epts )
            endif

            segment_endpoint(nseg) = ikpt
            dk_old = dk
         endif
      enddo
      nseg = nseg + 1
      if( nseg > max_segments ) then
         allocate( epts( nseg ) )
         epts = 0
         epts(1:nseg-1) = segment_endpoint(1:nseg-1)
         deallocate( segment_endpoint )
         allocate( segment_endpoint(nseg) )
         segment_endpoint = epts
         deallocate( epts )
      endif
      segment_endpoint(nseg) = NKPTS
      kpt_lines = nseg

      nkpts0 = count( l_good_kpt )
      if( nkpts0 /= NKPTS ) then
       write(IUO,'(/,5X,i0," repeated kpts; now total of ",i0/)') NKPTS-nkpts0,nkpts0
      endif

      allocate( x_kpt(NKPTS) )

      dk_old = bkpts(:,1)
      k0 = 1
      x0 = 0
      x_kpt(1) = x0
      kpt = 1
      do nseg = 1,kpt_lines
         do k = k0+1,segment_endpoint(nseg)
            if( .not. l_good_kpt(k) ) cycle
            kpt = kpt + 1
            dk = bkpts(:,k) - dk_old
            x_kpt(k) = x0 + sqrt( dot_product( dk, matmul( bbmat, dk ) ) )

            ! change to A^{-1}
            !x_kpt(k)=x_kpt(k)/AUTOA
         enddo
         x0 = x_kpt(segment_endpoint(nseg))
         k0 = segment_endpoint(nseg)
         dk_old = bkpts(:,segment_endpoint(nseg))
      enddo

      !--->    write to bs_projected.dat

      open(bskfh,file=bskfn,form='formatted',status='replace')
      write(bskfh,'(i0,2x,i0,f20.8,"  ! NKPTS, jspins, ef (eV)")')      &
     &       nkpts0,ISPIN,ef
      write(bskfh,'(i0,"  ! number of segments; positions:")') kpt_lines
      write(bskfh,'(6f12.6)')                                           &
     & x_kpt( segment_endpoint(1:kpt_lines) )/AUTOA

      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &      kpt,NBANDS,jspin,bkpts(:,ikpt)+KSHIFT(:),x_kpt(ikpt)/AUTOA
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo
      close( bskfh )

      end subroutine out_bands  
      end module 
