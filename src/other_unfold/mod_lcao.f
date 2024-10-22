Module mod_lcao
   use mod_comp
   use mod_input
   REAL(DP),ALLOCATABLE         :: proj_wt(:,:,:),CELTOT(:,:,:)
   REAL(DP),ALLOCATABLE         :: KPTVEC_Car(:,:),KPTVEC(:,:)
   INTEGER                      :: NKPTS,NBANDS
   contains

   subroutine lcao
   use mod_metric ,only:bmat
   COMPLEX(DP),ALLOCATABLE  :: C(:,:), C_high(:,:,:,:,:), S(:,:),S_high(:,:,:,:,:,:,:,:) !COEFFICIENT AND OVERLAP MATRIX 
   INTEGER                  :: I , J , K , M , Q , W , E , R  ,T ,P , I_DIRECT , I_SHIFT            !NORMA LOOP INDEX
   INTEGER                  :: NDIMENSION                                          !TOT NUMBER
   INTEGER                  :: IBAND , IKPOINT , ISPIN                             !SPECIAL LOOP INDEX
   CHARACTER (LEN = 1000)   :: S_filename,C_filename,KPOINT_CHARACTER_S,KPOINT_CHARACTER_C !ABACUS OUTPUT FILENAME
   INTEGER                  :: IKPOINT_TEMP                                        
   REAL(DP)                 :: C_REAL , C_IM                                       ! COEFFICIENT(REAL,IM)
   REAL(DP),ALLOCATABLE     :: NORMAL(:)
   COMPLEX(DP)              :: S_TEMP,WEIGHT_TEMP
!------------------------------------------------------------------------------------
   COMPLEX(DP)              :: PAULI_X(2,2),PAULI_Y(2,2),PAULI_Z(2,2)    !PAULI: 1:up 2:down PAULI_Y(1,1) UP UP 
   COMPLEX(DPC),ALLOCATABLE :: SPIN(:,:,:)
   REAL(DP)                 :: a(3,3),b(3,3)
   REAL(DP),ALLOCATABLE     :: proj_wt_type(:,:,:,:),w_temp(:,:,:,:)
   REAL(DP),ALLOCATABLE     :: KSHIFT_DIR(:,:)
   REAL(DP),ALLOCATABLE     :: KSHIFT_CAR(:,:)

!----------------------------------------   
   OPEN(UNIT=12,FILE="LOWF_K_1.dat")
   READ(12,*)
   READ(12,*)
   READ(12,*)NBAND
   CLOSE(12)
   !---------------!interface to out_bands
   NBANDS= NBAND   !
   NKPTS = NKPOINT*det !
   !---------------!
   ALLOCATE(CELTOT(NBAND,NKPOINT*det,1))
   ALLOCATE(SPIN(NKPOINT,NBAND,3))
   ALLOCATE(proj_wt(NBAND,NKPOINT*det,1),proj_wt_type(ions_type,NBAND,NKPOINT*det,1))
   ALLOCATE(w_temp(ions_type,NBAND,NKPOINT*det,1))
   ALLOCATE(KPTVEC(3,NKPOINT*det),KPTVEC_Car(3,NKPOINT*det))
   proj_wt=0
   proj_wt_type=0 
   KPOINT_DIRECT=0
   SPIN=0
! spin texture
   !x
   PAULI_X(1,1)=0
   PAULI_X(2,2)=0
   PAULI_X(1,2)=1
   PAULI_X(2,1)=1
   !y
   PAULI_Y(1,1)=0
   PAULI_Y(2,2)=0
   PAULI_Y(1,2)=-i
   PAULI_Y(2,1)=i
   !z
   PAULI_Z(1,1)=1
   PAULI_Z(2,2)=-1
   PAULI_Z(1,2)=0
   PAULI_Z(2,1)=0


! GET THE S-MATRIX AND COEFFICIENT

   DO IKPOINT = 1,NKPOINT!NKPOINT,1
      

      IKPOINT_TEMP=IKPOINT-1
      write(KPOINT_CHARACTER_S,*) IKPOINT_TEMP
      write(KPOINT_CHARACTER_C,*) IKPOINT
 
      S_filename = "data-"//trim(adjustl(KPOINT_CHARACTER_S))//"-S"  ! 0-nkpoint-1
      C_filename = "LOWF_K_"//trim(adjustl(KPOINT_CHARACTER_C))//".dat"  ! 1-nkpoint
      
      !=================== READ OVERLAP MATRIX =====================
      IF (spin_texture == .true.) THEN
      OPEN(UNIT=11,file=S_filename)
      READ(11,*) NDIMENSION
      CLOSE(11)
      ALLOCATE(S(NDIMENSION,NDIMENSION))
      S=0
      call READ_OVERLAP_MATRIX(IKPOINT,S,NDIMENSION,S_filename)
      ENDIF
     !
     !=================== READ COEFFICIENT ========================
      OPEN(UNIT=12,file=C_filename)
      READ(12,*)
      READ(12,*)KPOINT_CARTESIAN(IKPOINT,:)
     !
      READ(12,*)
      READ(12,*)orbital_tot
      CLOSE(12)
      ALLOCATE(C(NBAND,orbital_tot),NORMAL(NBAND)) 
      C=0
      call READ_COEFFICIENT(C_filename,NBAND,NKPOINT,orbital_tot,C,IKPOINT,CELTOT)
     !NORMAL
      NORMAL=0
      DO I =1, NBAND
         DO J = 1, orbital_tot
           NORMAL(I)=C(I,J)*CONJG(C(I,J))+NORMAL(I)
         ENDDO
      ENDDO
   
      DO I = 1,NBAND
         C(I,:)=C(I,:)/NORMAL(I)**0.5
      ENDDO
 
     !=================== END =====================================
     !
     !Ascending dimensions
     !
      ALLOCATE(C_high(NBAND,NSPIN,ions_type,maxval(ions_n),maxval(orbital_n)))
      C_high=0
     !
      DO IBAND = 1,NBAND!NBAND
         N=1   
         DO I = 1,ions_type
            DO J = 1 , ions_n(I)
               DO K = 1,orbital_n(I)
                  DO ISPIN = 1,NSPIN
                     C_high(IBAND,ISPIN,I,J,K)=C(IBAND,N) !Ascending dimensions
                     N=N+1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
     !
      IF (spin_texture == .true.) THEN
      ALLOCATE(S_high(NSPIN,ions_type,maxval(ions_n),maxval(orbital_n),NSPIN,ions_type,maxval(ions_n),maxval(orbital_n)))
      E=1
      S_high=0
      DO I =  1,ions_type ! 1,ions_type
         DO J = 1,ions_n(I)
            DO K = 1 , orbital_n(I)
               DO T=1,NSPIN
                  R=1! reset column_index,M*Q*W=NDIMENSION
                  DO M = 1 , ions_type
                     DO Q = 1 , ions_n(M)
                        DO W = 1 , orbital_n(M)
                           DO P =1,NSPIN
                           S_high(T,I,J,K,P,M,Q,W)=S(E,R) !Ascending dimensions
                           R=R+1  ! column_index+1
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  E=E+1 ! row_index+1
               ENDDO
            ENDDO
         ENDDO 
      ENDDO
      ENDIF
   ! 
   !start cal spin-texture
   !
   IF (spin_texture == .true.) then

   DO IBAND=1,NBAND
      DO I = 1,ions_type
         DO J = 1,ions_n(I)
            DO K = 1, orbital_n(I)
               DO  T= 1,NSPIN
                  DO M = 1, ions_type
                     DO Q = 1, ions_n(M)
                        DO W = 1, orbital_n(M)
                           DO P = 1,NSPIN
                           SPIN(IKPOINT,IBAND,1)=(C_high(IBAND,P,M,Q,W))*conjg(C_high(IBAND,T,I,J,K))*NORMAL(IBAND)*(S_high(T,I,J,K,T,M,Q,W))*PAULI_X(T,P)+SPIN(IKPOINT,IBAND,1)
                           SPIN(IKPOINT,IBAND,2)=(C_high(IBAND,P,M,Q,W))*conjg(C_high(IBAND,T,I,J,K))*NORMAL(IBAND)*(S_high(T,I,J,K,T,M,Q,W))*PAULI_Y(T,P)+SPIN(IKPOINT,IBAND,2)
                           SPIN(IKPOINT,IBAND,3)=(C_high(IBAND,P,M,Q,W))*conjg(C_high(IBAND,T,I,J,K))*NORMAL(IBAND)*(S_high(T,I,J,K,T,M,Q,W))*PAULI_Z(T,P)+SPIN(IKPOINT,IBAND,3)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ENDIF
  
   !
   !start bandunfolding 
   !
   IF (lcao_unfold == .true.) then
   S_TEMP = 0
   WEIGHT_TEMP = 0

   DO IBAND = 1,NBAND
      DO I = 1,ions_type
         DO K = 1 , orbital_n(I)
            DO ISPIN =1, NSPIN
               DO J= 1 , ions_n(I)
               WEIGHT_TEMP=WEIGHT_TEMP+C_high(IBAND,ISPIN,I,J,K)*exp(-TWOPI*im*dot_product(KPOINT_CARTESIAN(IKPOINT,:),POS_MAP_CARTESIAN(I,J,:)))
               ENDDO
               proj_wt_type(I,IBAND,IKPOINT,1)=proj_wt_type(I,IBAND,IKPOINT,1)+WEIGHT_TEMP*CONJG(WEIGHT_TEMP)
              !proj_wt(IBAND,IKPOINT,1)=WEIGHT_TEMP*CONJG(WEIGHT_TEMP)+proj_wt(IBAND,IKPOINT,1)
               WEIGHT_TEMP=0
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   ENDIF
 
   IF (spin_texture == .true.) then
      DEALLOCATE(S,S_high)
   ENDIF

   DEALLOCATE(C,C_high,NORMAL)

   ENDDO !KPOINT

   !OUT DATA
   IF (det/=1) THEN
    
   w_temp=proj_wt_type

   do i=2,det

      IF (i/2 == i/2.0) THEN
         do j=1,NKPOINT
            CELTOT(:,(i-1)*NKPOINT+j,1)=CELTOT(:,NKPOINT-j+1,1)
            proj_wt_type(:,:,(i-1)*NKPOINT+j,1)=w_temp(:,:,i*NKPOINT-j+1,1)
         enddo
      else
            CELTOT(:,1+(i-1)*NKPOINT:i*NKPOINT,1)=CELTOT(:,1:NKPOINT,1)
      endif

   enddo

   endif
   a(:,1)=bmat(:,1)
   a(:,2)=bmat(:,2)
   a(:,3)=bmat(:,3)
   call inverse_l(a,b)
   do i=1,NKPOINT*det
      KPTVEC(1,i)=KPTVEC_Car(1,i)*b(1,1)+KPTVEC_Car(2,i)*b(1,2)+KPTVEC_Car(3,i)*b(1,3)
      KPTVEC(2,i)=KPTVEC_Car(1,i)*b(2,1)+KPTVEC_Car(2,i)*b(2,2)+KPTVEC_Car(3,i)*b(2,3)
      KPTVEC(3,i)=KPTVEC_Car(1,i)*b(3,1)+KPTVEC_Car(2,i)*b(3,2)+KPTVEC_Car(3,i)*b(3,3)
   enddo

   IF (spin_texture == .true.) then
   OPEN(UNIT=14,FILE="spin_texture.dat")
   WRITE(14,*)"     IKPOINT        ","IBAND        ","SPIN_X        ","SPIN_Y        ","SPIN_Z        "
   DO IKPOINT=1,NKPOINT
      DO IBAND=1,NBAND
          WRITE(14,"(100F)")real(IKPOINT),real(IBAND),real(SPIN(IKPOINT,IBAND,1)),real(SPIN(IKPOINT,IBAND,2)),real(SPIN(IKPOINT,IBAND,3))
      ENDDO
   ENDDO
   CLOSE(14)
   ENDIF
   
   CELTOT=CELTOT*13.605826

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
   IF (lcao_unfold == .true.) then
      proj_wt=proj_wt/maxval(proj_wt)
      call out_bands_lcao
   ENDIF
 


   END subroutine lcao


SUBROUTINE READ_OVERLAP_MATRIX(IKPOINT,S,NDIMENSION,S_filename)
   INTEGER,PARAMETER                    :: DP =SELECTED_real_KIND(10)
   COMPLEX(DP),INTENT(OUT)              :: S(NDIMENSION,NDIMENSION)                         
   INTEGER                              :: I , J , K , N , M      !NORMA LOOP
   INTEGER                              :: NDIMENSION             !TOT NUMBER
   INTEGER,INTENT(IN)                   :: IKPOINT                !INDEX FOR LOOP 
   CHARACTER (LEN = 1000)               :: S_filename,KPOINT_CHARACTER
!
   S=0
   OPEN(UNIT=11,file=S_filename)

   READ(11,*)NDIMENSION,S(1,:)
   S(:,1)=CONJG(S(1,:))
   DO I = 2,NDIMENSION
   READ(11,*)S(I,I:NDIMENSION)
   S(I:NDIMENSION,I)=CONJG(S(I,I:NDIMENSION))
   ENDDO
   CLOSE(11)

END SUBROUTINE READ_OVERLAP_MATRIX

SUBROUTINE READ_COEFFICIENT(C_filename,NBAND,NKPOINT,orbital_tot,C,IKPOINT,EIGENV)
   INTEGER,PARAMETER                    :: DP =SELECTED_real_KIND(10)
   INTEGER                              :: I , J , K , N , M       !NORMA LOOP
   CHARACTER(LEN = 1000),INTENT(IN)     :: C_filename
   INTEGER                              :: ions_type , inos_n , orbital_type , orbital_n 
   INTEGER,INTENT(IN)                   :: NBAND, NKPOINT, orbital_tot, IKPOINT   
   INTEGER                              :: IBAND , I_TEN , I_MOD            !INDEX FOR LOOP 
   COMPLEX(DP),INTENT(OUT)              :: C(NBAND,orbital_tot)  
   REAL(DP)                             :: C_temp(2*orbital_tot)
   REAL(DP),INTENT(OUT)                 :: EIGENV(NBAND,NKPOINT,1)                          
  
   I_MOD=0
   I_TEN=0
   C_TEMP=0
   C=0
   OPEN(UNIT=12,file=C_filename)
   READ(12,*) !skip head
   READ(12,*) KPTVEC_Car(:,IKPOINT) 
   READ(12,*)
   READ(12,*)

   I_MOD = MOD(orbital_tot,5)
   I_TEN = (orbital_tot-I_MOD)/5
   
   DO IBAND = 1,NBAND
      READ(12,*)
      READ(12,*)EIGENV(IBAND,IKPOINT,1)
      READ(12,*)
      DO I = 1 ,I_TEN
         READ(12,*)C_temp((1+10*(I-1)):I*10)
      ENDDO
      IF(I_MOD/=0)THEN
      READ(12,*)C_temp((10*I_TEN+1):2*orbital_tot)
      ENDIF 
      DO I = 1, orbital_tot,1
         C(IBAND,I)=cmplx(C_temp(2*I-1),C_temp(2*I))
      ENDDO
   ENDDO
   CLOSE(12)

END SUBROUTINE READ_COEFFICIENT

subroutine out_bands_lcao
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
      integer              :: kpt_lines,nseg,k0,max_segments,ISPIN
      
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
     &   nkpts0,ISPIN,ef
      write(bskfh,'(i0,"  ! number of segments; positions:")') kpt_lines
      write(bskfh,'(6f12.6)')                                           &
     & x_kpt( segment_endpoint(1:kpt_lines) )/AUTOA

      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &      kpt,NBANDS,jspin,bkpts(:,ikpt),x_kpt(ikpt)/AUTOA
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo
      close( bskfh )

end subroutine out_bands_lcao

        
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
