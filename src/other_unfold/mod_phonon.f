Module mod_phonon
use mod_comp   
use mod_input
   REAL(DP),ALLOCATABLE     :: proj_wt(:,:,:),CELTOT(:,:,:)
   REAL(DP),ALLOCATABLE     :: KPTVEC(:,:),distance(:)
   INTEGER                  :: NKPTS,NBANDS,npath
   INTEGER,ALLOCATABLE      :: seg(:)

   contains

   SUBROUTINE unfold_phonon
   INTEGER                  :: I , J , K  , M , Q , W , E , R  ,T ,P , I_DIRECT             !NORMA LOOP INDEX
   COMPLEX(DP)              :: WEIGHT_TEMP,S_TEMP
   !-----------PHONON
   COMPLEX(DP),ALLOCATABLE  :: C_PHONON(:,:,:,:,:)
   INTEGER                  :: IDIR,ions_tot,natom
   REAL(DP),ALLOCATABLE     :: POSCAR_C(:,:),CONTCAR_C(:,:),POSCAR_D(:,:),CONTCAR_D(:,:),proj_wt_type(:,:,:,:)
   CHARACTER(LEN=1)         :: tchar
   logical                  :: expos,excont
   REAL(DP)                 :: b(3,3),b_inver(3,3)
   REAL                     :: real_vector,im_vector
   REAL(DP),ALLOCATABLE     :: NORMAL(:,:)
   NBAND=3*sum(ions_n)
   ALLOCATE(C_PHONON(NKPOINT,NBAND,ions_type,maxval(ions_n),3))
   ALLOCATE(proj_wt(NBAND,NKPOINT,1),CELTOT(NBAND,NKPOINT,1),proj_wt_type(ions_type,NBAND,NKPOINT,1))
   ALLOCATE(KPTVEC(3,NKPOINT),distance(NKPOINT))

   !interface to out_bands!
   NBANDS=NBAND           !
   NKPTS=NKPOINT          !
   !----------------------!

!======================================

!  inquire(file='POSCAR',exist=expos)
!  inquire(file='CONTCAR',exist=excont)
!  IF(expos==.false..or.excont==.false.) THEN
!     print*,'POSCAR or CONTCAR does not exist!'
!     stop
!  ENDIF

  ions_tot=sum(ions_n)
!  ALLOCATE(POSCAR_C(ions_tot,3),POSCAR_D(ions_tot,3))
!  ALLOCATE(CONTCAR_C(ions_tot,3),CONTCAR_D(ions_tot,3))
!  OPEN(UNIT=13,FILE='POSCAR')
!  READ(13,*)
!  READ(13,*)FACTOR
!  DO I =1,3,1
!  READ(13,*)CELL(I,:)
!  ENDDO
!  CELL=CELL*FACTOR
!  READ(13,*)
!  READ(13,*)
!  READ(13,*)tchar
!
!  IF (tchar =='S') THEN
!     READ(13,*)tchar
!     IF(tchar == 'D') THEN
!        DO I =1,ions_tot
!        READ(13,*)POSCAR_D(I,:)
!        ENDDO
!     ELSE IF(tchar=='C') THEN
!        DO I =1,ions_tot
!        READ(13,*)POSCAR_C(I,:)
!        ENDDO
!     ELSE
!     PRINT*,'ERRO WHEN READING POSCAR'
!     STOP
!     ENDIF!select dy
!  ELSE IF(tchar == 'D') THEN
!     DO I =1,ions_tot
!     READ(13,*)POSCAR_D(I,:)
!     ENDDO
!  ELSE IF(tchar=='C') THEN
!     DO I =1,ions_tot
!     READ(13,*)POSCAR_C(I,:)
!     ENDDO
!  ELSE
!     PRINT*,'ERRO WHEN READING POSCAR'
!     STOP
!  ENDIF
!  CLOSE(13)
!!---------------------------------------
!  OPEN(UNIT=13,FILE='CONTCAR')
!  DO I =1,7
!  READ(13,*)
!  ENDDO
!
!  READ(13,*)tchar
!
!  IF (tchar =='S') THEN
!     READ(13,*)tchar
!     IF(tchar == 'D') THEN
!        DO I =1,ions_tot
!        READ(13,*)CONTCAR_D(I,:)
!        ENDDO
!     ELSE IF(tchar=='C') THEN
!        DO I =1,ions_tot
!        READ(13,*)CONTCAR_C(I,:)
!        ENDDO
!     ELSE
!     PRINT*,'ERRO WHEN READING CONTCAR'
!     STOP
!     ENDIF!select dy
!  ELSE IF(tchar == 'D') THEN
!     DO I =1,ions_tot
!     READ(13,*)CONTCAR_D(I,:)
!     ENDDO
!  ELSE IF(tchar=='C') THEN
!     DO I =1,ions_tot
!     READ(13,*)CONTCAR_C(I,:)
!     ENDDO
!  ELSE
!     PRINT*,'ERRO WHEN READING CONTCAR'
!     STOP
!  ENDIF
!  CLOSE(13)
!
!  DO I=1,ions_tot
!     POSCAR_C(I,1)=CELL(1,1)*POSCAR_D(I,1)+CELL(2,1)*POSCAR_D(I,2)+CELL(3,1)*POSCAR_D(I,3)
!     POSCAR_C(I,2)=CELL(1,2)*POSCAR_D(I,1)+CELL(2,2)*POSCAR_D(I,2)+CELL(3,2)*POSCAR_D(I,3)
!     POSCAR_C(I,3)=CELL(1,3)*POSCAR_D(I,1)+CELL(2,3)*POSCAR_D(I,2)+CELL(3,3)*POSCAR_D(I,3)
!     CONTCAR_C(I,1)=CELL(1,1)*CONTCAR_D(I,1)+CELL(2,1)*CONTCAR_D(I,2)+CELL(3,1)*CONTCAR_D(I,3)
!     CONTCAR_C(I,2)=CELL(1,2)*CONTCAR_D(I,1)+CELL(2,2)*CONTCAR_D(I,2)+CELL(3,2)*CONTCAR_D(I,3)
!     CONTCAR_C(I,3)=CELL(1,3)*CONTCAR_D(I,1)+CELL(2,3)*CONTCAR_D(I,2)+CELL(3,3)*CONTCAR_D(I,3)
!  ENDDO
!  K=1
!  DO I=1,ions_type
!     DO J=1,ions_n(I)
!        POS_MAP_CARTESIAN(I,J,:)=POSCAR_C(K,:)-CONTCAR_C(K,:)
!        K=K+1
!     ENDDO
!  ENDDO

!======================================Get coeff

open(unit=17,file="band.yaml")
read(17,"(8x,I)")NKPTS
read(17,"(6x,I)")npath

allocate(seg(npath))

read(17,*)

do i=1,npath
  read(17,"(1X,I)")seg(i)
  if(i>1)then
     seg(i)=seg(i-1)+seg(i)
  endif
enddo

do i=1,(npath+1)+4,1
read(17,*)
enddo
read(17,"(6x,I)")natom
do i=1,5+(3*natom)+2
read(17,*)
enddo
do IKPOINT=1,NKPOINT
read(17,"(18X,F10.7,1X,F10.7,1X,F10.7)")KPTVEC(1,IKPOINT),KPTVEC(2,IKPOINT),KPTVEC(3,IKPOINT)
read(17,"(11X,F)")distance(IKPOINT)
read(17,*)

do IBAND=1,NBAND
read(17,*)
read(17,"(14X,F)")CELTOT(IBAND,IKPOINT,1)
read(17,*)

do I=1,ions_type
   do J=1,ions_n(I)
      read(17,*)
      do IDIR=1,3
      read(17,"(10X,F17.14,1X,F)")real_vector,im_vector
      C_PHONON(IKPOINT,IBAND,I,J,IDIR)=cmplx(real_vector,im_vector)
      enddo
   enddo
enddo

enddo
read(17,*)
enddo

CLOSE(17)

!NORM  !phonopy has normalized the coeff
!ALLOCATE(NORMAL(NKPOINT,NBAND))
!NORMAL=0
!DO IKPOINT =1,NKPOINT
!   DO IBAND = 1,NBAND
!      DO I=1,ions_type
!         DO J=1,ions_n(I)
!            DO IDIR=1,3
!               NORMAL(IKPOINT,IBAND)=C_PHONON(IKPOINT,IBAND,I,J,IDIR)*CONJG(C_PHONON(IKPOINT,IBAND,I,J,IDIR))+NORMAL(IKPOINT,IBAND)
!            ENDDO
!         ENDDO
!      ENDDO
!   ENDDO
!ENDDO
!
!print*,NORMAL
!DO IKPOINT=1,NKPOINT
!   DO IBAND=1,NBAND
!      C_PHONON(IKPOINT,IBAND,:,:,:)=C_PHONON(IKPOINT,IBAND,:,:,:)/NORMAL(NKPOINT,NBAND)**0.5
!   ENDDO
!ENDDO

!start phonon unfolding 
   DO IKPOINT =1,NKPOINT
      DO IBAND = 1,NBAND 
         DO I = 1,ions_type
            DO IDIR = 1, 3 !xyz
               DO J= 1, ions_n(I)
                  WEIGHT_TEMP=WEIGHT_TEMP+C_PHONON(IKPOINT,IBAND,I,J,IDIR)!*exp(-TWOPI*im*dot_product(KPOINT_CARTESIAN(IKPOINT,:),POS_MAP_CARTESIAN(I,J,:)))
               ENDDO
               proj_wt_type(I,IBAND,IKPOINT,1) = WEIGHT_TEMP*CONJG(WEIGHT_TEMP)+proj_wt_type(I,IBAND,IKPOINT,1)
               !proj_wt(IBAND,IKPOINT,1) = WEIGHT_TEMP*CONJG(WEIGHT_TEMP)+proj_wt(IBAND,IKPOINT,1)  
               WEIGHT_TEMP = 0
            ENDDO
         ENDDO
      ENDDO
   ENDDO    !K  
   CLOSE(16)
   IF (ALL(proj_type==0)) THEN
   DO I=1,ions_type
   proj_wt(:,:,1)=proj_wt_type(I,:,:,1)+proj_wt(:,:,1)
   ENDDO
   ELSE
   DO I=1,ions_type
   IF(proj_type(I)/=0) THEN
   proj_wt(:,:,1)=proj_wt_type(proj_type(I),:,:,1)+proj_wt(:,:,1)
   ENDIF
   ENDDO
   ENDIF
   proj_wt=proj_wt/maxval(proj_wt)
  
   call out_bands_phonon

END SUBROUTINE 

subroutine out_bands_phonon
      !
      use mod_comp
      use mod_metric, only:bbmat
      implicit none

!              wavecar file variables
      real(DP),allocatable :: e(:),occ(:)
      real(DP)             :: vkpt(3),rnv
      integer              :: nv

!              k-projected information
      real(DP),allocatable :: bkpts(:,:)

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
      integer              :: nseg,k0,max_segments,ISPIN,kpt_lines

      ISPIN=1
      allocate( e(NBANDS),occ(NBANDS))
      emin_ef = -1000
      emax_ef =  1000
      do jspin=1,ISPIN
        do ikpt=1,NKPTS
           e(1:NBANDS)   = CELTOT(1:NBANDS,ikpt,jspin)
           occ(1:NBANDS) = 1
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
         if( any( abs( dk - dk_old ) > 1.e-1_DP ) ) then
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

      write(bskfh,'(i0,"  ! number of segments; positions:")') npath
      write(bskfh,'(6f12.6)') distance( seg(1:npath) )

      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &      kpt,NBANDS,jspin,bkpts(:,ikpt),distance(ikpt)
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo
      close( bskfh )
end subroutine out_bands_phonon
        
SUBROUTINE inverse_l(b_,b_inverse)
     real(DP) :: omtil,b1_(3),b2_(3),b3_(3)
     real(DP),intent(out)::b_inverse(3,3)
     real(DP) ::b_(3,3)
     integer :: i,j
     b1_=b_(1:3,1)
     b2_=b_(1:3,2)
     b3_=b_(1:3,3)
     det_b_ = b1_(1)*b2_(2)*b3_(3) + b2_(1)*b3_(2)*b1_(3) +          &
     &         b3_(1)*b1_(2)*b2_(3) - b1_(3)*b2_(2)*b3_(1) -         &
     &         b2_(3)*b3_(2)*b1_(1) - b3_(3)*b1_(2)*b2_(1)
     ! b_=transpose(b_)
     b_inverse(1,1)=  b_(2,2)*b_(3,3)-b_(2,3)*b_(3,2)
     b_inverse(1,2)=-(b_(2,1)*b_(3,3)-b_(2,3)*b_(3,1))
     b_inverse(1,3)=  b_(2,1)*b_(3,2)-b_(2,2)*b_(3,1)
     b_inverse(2,1)=-(b_(1,2)*b_(3,3)-b_(1,3)*b_(3,2))
     b_inverse(2,2)=  b_(1,1)*b_(3,3)-b_(1,3)*b_(3,1)
     b_inverse(2,3)=-(b_(1,1)*b_(3,2)-b_(1,2)*b_(3,1))
     b_inverse(3,1)=  b_(1,2)*b_(2,3)-b_(1,3)*b_(2,2)
     b_inverse(3,2)=-(b_(1,1)*b_(2,3)-b_(1,3)*b_(2,1))
     b_inverse(3,3)=  b_(1,1)*b_(2,2)-b_(2,1)*b_(1,2)
     b_inverse=b_inverse/det_b_
     do i=1,3
       do j=1,3
       if(abs(b_inverse(i,j))<1E-7) b_inverse(i,j)=0
       end do
     end do
END SUBROUTINE inverse_l

END MODULE

 

