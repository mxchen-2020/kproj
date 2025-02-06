PROGRAM main 
   INTEGER,PARAMETER        :: DP =SELECTED_real_KIND(10)
   INTEGER                  :: ions_type , orbital_tot
   INTEGER,ALLOCATABLE      :: ions_n(:) , orbital_n(:)                            !number of ions_x  , number of  ions_type_x' orbital
   COMPLEX(DP),ALLOCATABLE  :: C(:,:), C_high(:,:,:,:), S(:,:),S_high(:,:,:,:,:,:) !COEFFICIENT AND OVERLAP MATRIX 
   INTEGER                  :: I , J , K , N , M , Q , W , E , R  ,T               !NORMA LOOP INDEX
   INTEGER                  :: NBAND , NKPOINT , NDIMENSION                        !TOT NUMBER
   INTEGER                  :: IBAND , IKPOINT                                     !SPECIAL LOOP INDEX
   CHARACTER (LEN = 1000)   :: S_filename,C_filename,KPOINT_CHARACTER_S,KPOINT_CHARACTER_C !ABACUS OUTPUT FILENAME
   INTEGER                  :: IKPOINT_TEMP                                        
   REAL(DP),ALLOCATABLE     :: WEIGHT(:,:),EIGENV(:,:)
   REAL(DP),ALLOCATABLE     :: KPOINT_DIRECT(:,:),KPOINT_CARTESIAN(:,:)
   REAL(DP)                 :: C_REAL , C_IM                                       ! COEFFICIENT(REAL,IM)
   INTEGER,ALLOCATABLE      :: POS_MAP_DIRECT(:,:,:)
   REAL(DP),ALLOCATABLE     :: POS_MAP_CARTESIAN(:,:,:),NORMAL(:)
   REAL(DP)                 :: CELL(3,3),RECELL(3,3),POS(2,4,3)
   COMPLEX(DP)              :: im,WEIGHT_TEMP,S_TEMP
   REAL(DP),PARAMETER       :: PI = 3.141592653589793238462643383279502884197_DP
   REAL(DP),PARAMETER       :: TWOPI = 2._DP*PI

!caes for 2x2 graphene=================!INPUT
  
   CELL(1,:)=(/ 2.4638104439 ,  0.0000000000 ,  0.0000000000/) 
   CELL(2,:)=(/ 1.2345617743 ,  2.1352538856 ,  0.0000000000/)
   CELL(3,:)=(/ 0.0000000000 ,  0.0000000000 ,  15.4974455833/)

   ions_type = 2
   OPEN(UNIT=13,FILE="WEIGHT")   


   ALLOCATE(ions_n(ions_type))
   ions_n=(/24,25/)

   ALLOCATE(orbital_n(ions_type))
   orbital_n=(/13,13/)
   !这里有一个前提，被归为一类的原子应有一样的轨道
   ALLOCATE(POS_MAP_DIRECT(ions_type,maxval(ions_n),3),POS_MAP_CARTESIAN(ions_type,maxval(ions_n),3))


   OPEN(UNIT=9,FILE="pos_map.dat")
   DO I = 1,ions_type
     DO J = 1,ions_n(I)
       READ(9,*)POS_MAP_DIRECT(I,J,:)
       POS_MAP_CARTESIAN(I,J,1)=POS_MAP_DIRECT(I,J,1)*CELL(1,1)+POS_MAP_DIRECT(I,J,2)*CELL(2,1)+POS_MAP_DIRECT(I,J,3)*CELL(3,1)
       POS_MAP_CARTESIAN(I,J,2)=POS_MAP_DIRECT(I,J,1)*CELL(1,2)+POS_MAP_DIRECT(I,J,2)*CELL(2,2)+POS_MAP_DIRECT(I,J,3)*CELL(3,2)
       POS_MAP_CARTESIAN(I,J,3)=POS_MAP_DIRECT(I,J,1)*CELL(1,3)+POS_MAP_DIRECT(I,J,2)*CELL(2,3)+POS_MAP_DIRECT(I,J,3)*CELL(3,3)
     ENDDO
   ENDDO


!=======================================
   im=(0,1)
   NKPOINT=200
   NBAND=117
   ALLOCATE(EIGENV(NKPOINT,NBAND))
   ALLOCATE(WEIGHT(NKPOINT,NBAND))
   ALLOCATE(KPOINT_DIRECT(NKPOINT,3),KPOINT_CARTESIAN(NKPOINT,3))
   WEIGHT=0 
   KPOINT_DIRECT=0
  

! GET THE S-MATRIX AND COEFFICIENT

   DO IKPOINT = 1,NKPOINT,1


      IKPOINT_TEMP=IKPOINT-1
      write(KPOINT_CHARACTER_S,*) IKPOINT_TEMP
      write(KPOINT_CHARACTER_C,*) IKPOINT
 
      S_filename = "data-"//trim(adjustl(KPOINT_CHARACTER_S))//"-S"  ! 0-nkpoint-1
      C_filename = "LOWF_K_"//trim(adjustl(KPOINT_CHARACTER_C))//".dat"  ! 1-nkpoint
      
      !=================== READ OVERLAP MATRIX =====================
      OPEN(UNIT=11,file=S_filename)
      READ(11,*) NDIMENSION
      CLOSE(11)
      ALLOCATE(S(NDIMENSION,NDIMENSION))
      S=0
      call READ_OVERLAP_MATRIX(IKPOINT,S,NDIMENSION,S_filename)
     !
     !=================== READ COEFFICIENT ========================
      OPEN(UNIT=12,file=C_filename)
      READ(12,*)
      READ(12,*)KPOINT_CARTESIAN(IKPOINT,:)
     !
      READ(12,*)NBAND
      READ(12,*)orbital_tot
      CLOSE(12)
      ALLOCATE(C(NBAND,orbital_tot),NORMAL(NBAND)) 
      C=0
      call READ_COEFFICIENT(C_filename,NBAND,NKPOINT,orbital_tot,C,IKPOINT,EIGENV)

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
      ALLOCATE(C_high(NBAND,ions_type,maxval(ions_n),maxval(orbital_n)))
      C_high=0
     !
      DO IBAND = 1,NBAND
         N=1   
         DO I = 1,ions_type
            DO J = 1 , ions_n(I)
               DO K = 1,orbital_n(I)
                  C_high(IBAND,I,J,K)=C(IBAND,N) !Ascending dimensions
                  N=N+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
     !
      ALLOCATE(S_high(ions_type,maxval(ions_n),maxval(orbital_n),ions_type,maxval(ions_n),maxval(orbital_n)))
      E=1
      S_high=0
      DO I = 1,ions_type
         DO J = 1,ions_n(I)
            DO K = 1 , orbital_n(I)
               R=1! reset column_index,M*Q*W=NDIMENSION
               DO M = 1 , ions_type
                  DO Q = 1 , ions_n(M)
                     DO W = 1 , orbital_n(M)
                        S_high(I,J,K,M,Q,W)=S(E,R) !Ascending dimensions
                        R=R+1  ! column_index+1
                     ENDDO
                  ENDDO
               ENDDO
               E=E+1 ! row_index+1
            ENDDO
         ENDDO
      ENDDO 
   ! 
   !start k-projection
   WEIGHT_TEMP = 0
   S_TEMP = 0 
   DO IBAND = 1,NBAND 
      DO I = 1, ions_type
         DO K = 1 , orbital_n(I)
            DO J= 1 , ions_n(I) 
                  WEIGHT_TEMP=WEIGHT_TEMP+C_high(IBAND,I,J,K)*exp(-TWOPI*im*dot_product(KPOINT_CARTESIAN(IKPOINT,:),POS_MAP_CARTESIAN(I,J,:)))
            ENDDO
               WEIGHT(IKPOINT,IBAND)=WEIGHT_TEMP*CONJG(WEIGHT_TEMP)+WEIGHT(IKPOINT,IBAND)
               WEIGHT_TEMP=0
         ENDDO
      ENDDO
   ENDDO
    

   DEALLOCATE(C,S,C_high,S_high,NORMAL)
   ENDDO ! KPOINT
   

OPEN(UNIT=14,FILE="bs_project.dat")
WRITE(14,*)" @version 50120           "                  
WRITE(14,*)" @page size 792, 612"
WRITE(14,*)" @with line"
WRITE(14,*)" @    line on"
WRITE(14,*)" @    line loctype world"
WRITE(14,*)" @    line g0"
WRITE(14,*)" @    line  0, 0,    3.177394, 0"
WRITE(14,*)" @    line linewidth 1.0"
WRITE(14,*)" @    line linestyle 1"
WRITE(14,*)" @    line color 1"
WRITE(14,*)" @    line arrow 0"
WRITE(14,*)" @    line arrow type 0"
WRITE(14,*)" @    line arrow length 1.0"
WRITE(14,*)" @    line arrow layout 1.0, 1.0"
WRITE(14,*)" @line def"
WRITE(14,*)" @with g0"
WRITE(14,*)" @    world 0, -20,    3.177394, 10"
WRITE(14,*)" @    view  0.25, 0.20, 1.00, 0.85"
WRITE(14,*)" @    xaxis  on"
WRITE(14,*)" @    xaxis  tick off"
WRITE(14,*)" @    xaxis  ticklabel off"
WRITE(14,*)" @    yaxis  on"
WRITE(14,*)" @    yaxis  tic on"
WRITE(14,*)" @    yaxis  tick major 5"
WRITE(14,*)" @    yaxis  tick minor ticks 1"
WRITE(14,*)" @    s0 type xysize"
WRITE(14,*)" @    s0 symbol 1"
WRITE(14,*)" @    s0 symbol sie 0.5"
WRITE(14,*)" @    s0 symbol color 11"
WRITE(14,*)" @    s0 symbol linewidth 1.0"
WRITE(14,*)" @    s0 line type 0"
WRITE(14,*)" @    s0 line linestyle 1"
WRITE(14,*)" @    s0 line linewidth 1.0"
WRITE(14,*)" @target g0.s0"
WRITE(14,*)" @type xysize"
WRITE(14,*)" "

   DO IKPOINT =1,NKPOINT
      DO IBAND=1,NBAND
         WRITE(14,"(100F)")real(IKPOINT),EIGENV(IKPOINT,IBAND)*13.605826,WEIGHT(IKPOINT,IBAND)/50
      ENDDO
   ENDDO
   CLOSE(14)


END

 
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
   REAL(DP),INTENT(OUT)                 :: EIGENV(NKPOINT,NBAND)                          
  
   I_MOD=0
   I_TEN=0
   C_TEMP=0
   C=0
   OPEN(UNIT=12,file=C_filename)
   DO I = 1,4
      READ(12,*) !skip head
   ENDDO

   I_MOD = MOD(orbital_tot,5)
   I_TEN = (orbital_tot-I_MOD)/5
   
   DO IBAND = 1,NBAND
      READ(12,*)
      READ(12,*)EIGENV(IKPOINT,IBAND)
      READ(12,*)
      DO I = 1 ,I_TEN
         READ(12,*)C_temp((1+10*(I-1)):I*10)
      ENDDO
      READ(12,*)C_temp((10*I_TEN+1):2*orbital_tot)

      DO I = 1, orbital_tot,1
         C(IBAND,I)=cmplx(C_temp(2*I-1),C_temp(2*I))
      ENDDO
   ENDDO
   CLOSE(12)
END SUBROUTINE 
  


SUBROUTINE CAL_RECELL(CELL,RECELL)
INTEGER,PARAMETER        :: DP =SELECTED_real_KIND(10)
REAL(DP)                 :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),omtil
REAL(DP)                 :: CELL(3,3)
REAL(DP)                 :: RECELL(3,3)

      a1(1:3) = CELL(1,1:3)
      a2(1:3) = CELL(2,1:3)
      a3(1:3) = CELL(3,1:3)

!  volume
      omtil = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) +                    &
     &        a3(1)*a1(2)*a2(3) - a1(3)*a2(2)*a3(1) -                    &
     &        a2(3)*a3(2)*a1(1) - a3(3)*a1(2)*a2(1)

!  reciprocal lattice vectors (no 2\pi)
      b1(1) = (a2(2)*a3(3) - a2(3)*a3(2))/omtil
      b1(2) = (a2(3)*a3(1) - a2(1)*a3(3))/omtil
      b1(3) = (a2(1)*a3(2) - a2(2)*a3(1))/omtil
      b2(1) = (a3(2)*a1(3) - a3(3)*a1(2))/omtil
      b2(2) = (a3(3)*a1(1) - a3(1)*a1(3))/omtil
      b2(3) = (a3(1)*a1(2) - a3(2)*a1(1))/omtil
      b3(1) = (a1(2)*a2(3) - a1(3)*a2(2))/omtil
      b3(2) = (a1(3)*a2(1) - a1(1)*a2(3))/omtil
      b3(3) = (a1(1)*a2(2) - a1(2)*a2(1))/omtil

      RECELL(1,:)=b1(:)
      RECELL(2,:)=b2(:)
      RECELL(3,:)=b3(:)
END SUBROUTINE CAL_RECELL


