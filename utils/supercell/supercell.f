! Generate a supercell out of a generating cell
      PROGRAM SUPERCELL
      USE mod_comp
      USE mod_poscar

      IMPLICIT NONE

      REAL(DP), PARAMETER :: TINY=1.E-10_DP

      LOGICAL     :: LFOUND
      INTEGER     :: NCELLS,NCOUNT,NCOUNT_TOT
      INTEGER     :: NX,NY,NZ,I,J,I1,I2,I3
      REAL(DP)    :: SCALE,VOL,OMEGA,VOLO,OMEGAO
      REAL(DP)    :: AO1B1,AO2B1,AO3B1,AO1B2,AO2B2,AO3B2,AO1B3,AO2B3,AO3B3
      REAL(DP)    :: X(3),Y(3),AOLD(3,3),A(3,3),BOLD(3,3),B(3,3),POSOLD(3)
      REAL(DP)    :: M(3,3),Y_OLD(3) 
      REAL(DP)    :: SHIFT_SC_A(3),SHIFT_PC_A(3,3)!jx-add
      REAL(DP), ALLOCATABLE :: POS_TRANS_OLD(:,:)!jx-add
      REAL(DP), ALLOCATABLE :: POS_TRANS(:,:),POS_C(:,:)
      integer,allocatable :: id_map(:,:),ICOUNT(:)
      integer              :: Irange , iions
!read POSCAR for the primitive cell
      call poscar()
      AOLD(:,:) = latt(:,:)
      CALL RECIPS(1._DP,AOLD(:,1),AOLD(:,2),AOLD(:,3), &
                       BOLD(:,1),BOLD(:,2),BOLD(:,3))
      CALL CELVOL(AOLD(:,1),AOLD(:,2),AOLD(:,3),VOLO)
      OMEGAO = ABS(VOLO)
      !!--jx-add range of ele
      ALLOCATE(ICOUNT(ntype))
      DO I=1,ntype
      ICOUNT(I)=sum(nitype(1:I))
      ENDDO
      iions=0


!read transform matrix that transforms a primitive cell into a supercell
      WRITE(*,'(3X,A)') "Transformation matrix (must be integers):"
      WRITE(*,'(5X,A)') "M(1,1)  M(2,1)  M(3,1)"
      WRITE(*,'(5X,A)') "M(1,2)  M(2,2)  M(3,2)"
      WRITE(*,'(5X,A)') "M(1,3)  M(2,3)  M(3,3)"
      READ(*,*) M(1,1),M(2,1),M(3,1)
      READ(*,*) M(1,2),M(2,2),M(3,2)
      READ(*,*) M(1,3),M(2,3),M(3,3)

!get lattice vectors of the supercell
      A(1,1) = M(1,1)*AOLD(1,1) + M(2,1)*AOLD(1,2) + M(3,1)*AOLD(1,3)
      A(2,1) = M(1,1)*AOLD(2,1) + M(2,1)*AOLD(2,2) + M(3,1)*AOLD(2,3)
      A(3,1) = M(1,1)*AOLD(3,1) + M(2,1)*AOLD(3,2) + M(3,1)*AOLD(3,3)

      A(1,2) = M(1,2)*AOLD(1,1) + M(2,2)*AOLD(1,2) + M(3,2)*AOLD(1,3)
      A(2,2) = M(1,2)*AOLD(2,1) + M(2,2)*AOLD(2,2) + M(3,2)*AOLD(2,3)
      A(3,2) = M(1,2)*AOLD(3,1) + M(2,2)*AOLD(3,2) + M(3,2)*AOLD(3,3)

      A(1,3) = M(1,3)*AOLD(1,1) + M(2,3)*AOLD(1,2) + M(3,3)*AOLD(1,3)
      A(2,3) = M(1,3)*AOLD(2,1) + M(2,3)*AOLD(2,2) + M(3,3)*AOLD(2,3)
      A(3,3) = M(1,3)*AOLD(3,1) + M(2,3)*AOLD(3,2) + M(3,3)*AOLD(3,3)

      CALL RECIPS(1._DP,A(:,1),A(:,2),A(:,3),B(:,1),B(:,2),B(:,3))
      CALL CELVOL(A(:,1),A(:,2),A(:,3),VOL)
      OMEGA=ABS(VOL)
      IF (OMEGA < (TINY*TINY*TINY)) THEN
         WRITE(*,'(3X,A)') 'ERROR: Volume of new cell is zero!'
         STOP
      ENDIF
      IF (ABS(OMEGA/OMEGAO-NINT(OMEGA/OMEGAO)) > TINY) THEN
         WRITE(*,'(3X,A)') 'ERROR1: Incommensurable lattice vectors!'
         stop 
      ENDIF

      NCELLS=NINT(OMEGA/OMEGAO)
      WRITE(*,'(3X,A,I4,2X,A)') &
     & "The supercell contains",NCELLS,"primitive cells"

! Commensurability check:
      X(1)=A(1,1)*BOLD(1,1)+A(2,1)*BOLD(2,1)+A(3,1)*BOLD(3,1)
      X(2)=A(1,1)*BOLD(1,2)+A(2,1)*BOLD(2,2)+A(3,1)*BOLD(3,2)
      X(3)=A(1,1)*BOLD(1,3)+A(2,1)*BOLD(2,3)+A(3,1)*BOLD(3,3)
      IF ((ABS(X(1)-NINT(X(1))) > TINY) .OR. &
          (ABS(X(2)-NINT(X(2))) > TINY) .OR. &
          (ABS(X(3)-NINT(X(3))) > TINY)) THEN
         WRITE(*,'(3X,A)') 'ERROR2: Incommensurable lattice vectors!'
         WRITE(*,'(3X,A)') 'HINT: CHECK THE TRANSFORM MATRIX.'
         stop 
      ENDIF
      X(1)=A(1,2)*BOLD(1,1)+A(2,2)*BOLD(2,1)+A(3,2)*BOLD(3,1)
      X(2)=A(1,2)*BOLD(1,2)+A(2,2)*BOLD(2,2)+A(3,2)*BOLD(3,2)
      X(3)=A(1,2)*BOLD(1,3)+A(2,2)*BOLD(2,3)+A(3,2)*BOLD(3,3)
      IF ((ABS(X(1)-NINT(X(1))) > TINY) .OR. &
          (ABS(X(2)-NINT(X(2))) > TINY) .OR. &
          (ABS(X(3)-NINT(X(3))) > TINY)) THEN
         WRITE(*,'(3X,A)') 'ERROR3: Incommensurable lattice vectors!'
         WRITE(*,'(3X,A)') 'HINT: CHECK THE TRANSFORM MATRIX.'
         stop 
      ENDIF
      X(1)=A(1,3)*BOLD(1,1)+A(2,3)*BOLD(2,1)+A(3,3)*BOLD(3,1)
      X(2)=A(1,3)*BOLD(1,2)+A(2,3)*BOLD(2,2)+A(3,3)*BOLD(3,2)
      X(3)=A(1,3)*BOLD(1,3)+A(2,3)*BOLD(2,3)+A(3,3)*BOLD(3,3)
      IF ((ABS(X(1)-NINT(X(1))) > TINY) .OR. &
          (ABS(X(2)-NINT(X(2))) > TINY) .OR. &
          (ABS(X(3)-NINT(X(3))) > TINY)) THEN
         WRITE(*,'(3X,A)') 'ERROR4: Incommensurable lattice vectors!'
         WRITE(*,'(3X,A)') 'HINT: CHECK THE TRANSFORM MATRIX.'
         stop 
      ENDIF

      open(1,file="POSCAR_SUPER",status='unknown')
      open(3,file="POSCAR_SUPER.cart",status='unknown')
      open(9,file="pos_map.dat",status='unknown')
       WRITE(9,'(A)')"    1  !number of supercells"
       WRITE(9,'(I5,A)')nions,"  !number of atoms in each primitive cell"
       WRITE(9,'(I5,A)')nions*NCELLS,"  !number of atoms in each supercel"
!      WRITE(9,*)alatt
!      WRITE(9,*) M(1,1),M(2,1),M(3,1)
!      WRITE(9,*) M(1,2),M(2,2),M(3,2)
!      WRITE(9,*) M(1,3),M(2,3),M(3,3)
!      WRITE(9,*) A(:,1)
!      WRITE(9,*) A(:,2)
!      WRITE(9,*) A(:,3)

! read POSOLD and write to POSNEW
! All seems okay ...
      WRITE(1,'(A)') "supercell"
      WRITE(1,'(F22.16)') 1.0_DP
      WRITE(1,'(3F22.16)') A(1,1),A(2,1),A(3,1)
      WRITE(1,'(3F22.16)') A(1,2),A(2,2),A(3,2)
      WRITE(1,'(3F22.16)') A(1,3),A(2,3),A(3,3)
      if(vasp4)then
       WRITE(1,'(20I9)') (nitype(I)*NCELLS,I=1,ntype)
      else
      WRITE(1,'(20A4)') (ele_name(I),I=1,ntype)
      WRITE(1,'(20I9)') (nitype(I)*NCELLS,I=1,ntype)
      endif
      if(L_SELECT_DYN)write(1,'(A)')"Selective dynamics"
      WRITE(1,'(A)') 'Direct'

      WRITE(3,'(A)') "supercell"
      WRITE(3,'(F22.16)') 1.0_DP
      WRITE(3,'(3F22.16)') A(1,1),A(2,1),A(3,1)
      WRITE(3,'(3F22.16)') A(1,2),A(2,2),A(3,2)
      WRITE(3,'(3F22.16)') A(1,3),A(2,3),A(3,3)
      if(vasp4)then
        WRITE(3,'(20I4)') (nitype(I)*NCELLS,I=1,ntype)
      else
        WRITE(3,'(20A4)') (ele_name(I),I=1,ntype)
        WRITE(3,'(20I4)') (nitype(I)*NCELLS,I=1,ntype)
      endif
      if(L_SELECT_DYN)write(3,'(A)')"Selective dynamics"
      WRITE(3,'(A)') 'Cartessian'

! How many times must we repeat the cell in each direction?
      AO1B1 = AOLD(1,1)*B(1,1) + AOLD(2,1)*B(2,1) + AOLD(3,1)*B(3,1)
      AO2B1 = AOLD(1,1)*B(1,2) + AOLD(2,1)*B(2,2) + AOLD(3,1)*B(3,2)
      AO3B1 = AOLD(1,1)*B(1,3) + AOLD(2,1)*B(2,3) + AOLD(3,1)*B(3,3)
      AO1B2 = AOLD(1,2)*B(1,1) + AOLD(2,2)*B(2,1) + AOLD(3,2)*B(3,1)
      AO2B2 = AOLD(1,2)*B(1,2) + AOLD(2,2)*B(2,2) + AOLD(3,2)*B(3,2)
      AO3B2 = AOLD(1,2)*B(1,3) + AOLD(2,2)*B(2,3) + AOLD(3,2)*B(3,3)
      AO1B3 = AOLD(1,3)*B(1,1) + AOLD(2,3)*B(2,1) + AOLD(3,3)*B(3,1)
      AO2B3 = AOLD(1,3)*B(1,2) + AOLD(2,3)*B(2,2) + AOLD(3,3)*B(3,2)
      AO3B3 = AOLD(1,3)*B(1,3) + AOLD(2,3)*B(2,3) + AOLD(3,3)*B(3,3)
      IF (ABS(AO1B1) < TINY) AO1B1 = 1.E30_DP
      IF (ABS(AO1B2) < TINY) AO1B2 = 1.E30_DP
      IF (ABS(AO1B3) < TINY) AO1B3 = 1.E30_DP
      IF (ABS(AO2B1) < TINY) AO2B1 = 1.E30_DP
      IF (ABS(AO2B2) < TINY) AO2B2 = 1.E30_DP
      IF (ABS(AO2B3) < TINY) AO2B3 = 1.E30_DP
      IF (ABS(AO3B1) < TINY) AO3B1 = 1.E30_DP
      IF (ABS(AO3B2) < TINY) AO3B2 = 1.E30_DP
      IF (ABS(AO3B3) < TINY) AO3B3 = 1.E30_DP
      NX = NINT(MAX( MAX(1._DP/ABS(AO1B1),1._DP/ABS(AO2B1)), 1._DP/ABS(AO3B1)))
      NY = NINT(MAX( MAX(1._DP/ABS(AO1B2),1._DP/ABS(AO2B2)), 1._DP/ABS(AO3B2)))
      NZ = NINT(MAX( MAX(1._DP/ABS(AO1B3),1._DP/ABS(AO2B3)), 1._DP/ABS(AO3B3)))
      ALLOCATE(POS_TRANS(3,NX*NY*NZ),POS_TRANS_OLD(3,NX*NY*NZ))
      allocate(POS_C(3,NX*NY*NZ))
      allocate(id_map(3,NCELLS))

! main loop (implicit) for reading and writing atomic positions
      NCOUNT_TOT = 0
      do I=1,nions
! Get cartesian coordinates
         POSOLD(1:3) = posca(1:3,I)
         NCOUNT = 0
! Repeat positions ...
      DO I3=0,NZ-1
       DO I2=0,NY-1
        DO I1=0,NX-1
! ...by adding lattice vectors ...
         X(1) = POSOLD(1) + I1*AOLD(1,1) + I2*AOLD(1,2) + I3*AOLD(1,3)
         X(2) = POSOLD(2) + I1*AOLD(2,1) + I2*AOLD(2,2) + I3*AOLD(2,3)
         X(3) = POSOLD(3) + I1*AOLD(3,1) + I2*AOLD(3,2) + I3*AOLD(3,3)
! Back to direct coordinates
         Y(1) = X(1)*B(1,1) + X(2)*B(2,1) + X(3)*B(3,1)
         Y(2) = X(1)*B(1,2) + X(2)*B(2,2) + X(3)*B(3,2)
         Y(3) = X(1)*B(1,3) + X(2)*B(2,3) + X(3)*B(3,3)
! Periodic boundary conditions, move to interval [0.,1.)
         Y_OLD=Y
         Y(1)=MOD(MOD(Y(1),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
         Y(2)=MOD(MOD(Y(2),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
         Y(3)=MOD(MOD(Y(3),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
! Does this position already exist?
         LFOUND = .FALSE.
         DO J=1,NCOUNT
            LFOUND = LFOUND .OR. ((ABS(Y(1)-POS_TRANS(1,J)) < TINY) .AND. &
     &                            (ABS(Y(2)-POS_TRANS(2,J)) < TINY) .AND. &
     &                            (ABS(Y(3)-POS_TRANS(3,J)) < TINY))
         ENDDO
! If not store this position
         IF (.NOT.LFOUND) THEN
            NCOUNT = NCOUNT + 1
            POS_TRANS(1,NCOUNT) = Y(1)  
            POS_TRANS(2,NCOUNT) = Y(2)
            POS_TRANS(3,NCOUNT) = Y(3)
          !jx-add 
            SHIFT_SC_A(1)=Y_OLD(1)-Y(1)
            SHIFT_SC_A(2)=Y_OLD(2)-Y(2)
            SHIFT_SC_A(3)=Y_OLD(3)-Y(3)
          ! print*,SHIFT_SC_A
          !
            SHIFT_PC_A(1,1)= SHIFT_SC_A(1)*M(1,1)     ! nA1=n*(Ia1+Ja2+Ka3)
            SHIFT_PC_A(1,2)= SHIFT_SC_A(1)*M(2,1)
            SHIFT_PC_A(1,3)= SHIFT_SC_A(1)*M(3,1)
            
            SHIFT_PC_A(2,1)= SHIFT_SC_A(2)*M(1,2)
            SHIFT_PC_A(2,2)= SHIFT_SC_A(2)*M(2,2)
            SHIFT_PC_A(2,3)= SHIFT_SC_A(2)*M(3,2)

            SHIFT_PC_A(3,1)= SHIFT_SC_A(3)*M(1,3)
            SHIFT_PC_A(3,2)= SHIFT_SC_A(3)*M(2,3)
            SHIFT_PC_A(3,3)= SHIFT_SC_A(3)*M(3,3)
           ! print*,SHIFT_PC_A(1,:)
           ! print*,SHIFT_PC_A(2,:)
           ! print*,SHIFT_PC_A(3,:)
           ! PRINT*,"------------------"
            call dire2car(POS_TRANS(:,NCOUNT),POS_C(:,NCOUNT),A)

            id_map(1,NCOUNT)= I1-(nint(SHIFT_PC_A(1,1))+&
            nint(SHIFT_PC_A(2,1))+nint(SHIFT_PC_A(3,1)))

            id_map(2,NCOUNT)= I2-(nint(SHIFT_PC_A(1,2))+&
            nint(SHIFT_PC_A(2,2))+nint(SHIFT_PC_A(3,2)))

            id_map(3,NCOUNT)= I3-(nint(SHIFT_PC_A(1,3))+&
            nint(SHIFT_PC_A(2,3))+nint(SHIFT_PC_A(3,3)))
            SHIFT_SC_A=0
            SHIFT_PC_A=0
            Y_OLD    = 0
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      if(NCOUNT/=NCELLS)then
        write(*,'(3X,A)') &
     & "ERROR: I CANNOT FIND THE RIGHT NUMBER OS TRANSLATIONS"
        stop
      endif
! write all positions generated out of the initial position

     
      DO J=1,NCOUNT
         iions=iions+1
         IF (L_SELECT_DYN) THEN
            WRITE(1,'(3F22.16,3L4,I)') POS_TRANS(1:3,J),dynmics(1:3,I),iions
            WRITE(3,'(3F22.16,3L4,I)') POS_C(1:3,J),dynmics(1:3,I),iions
         ELSE
            WRITE(1,'(3F22.16,I)')POS_TRANS(1:3,J),iions
            WRITE(3,'(3F22.16,I)')POS_C(1:3,J),iions
         ENDIF
        ! WRITE(9,'(I4,3F20.14,2X,I4,3F20.14,2X,3I3)')I,POSOLD(1:3), &
     ! & NCOUNT_TOT+J,POS_C(1:3,J),id_map(1:3,J)  -JX-DELETE
         !print ele_name
         !
         IF (I <= ICOUNT(1))Then
         WRITE(9,'(I3,3X,3I3,3X,A,I)')I,id_map(1:3,J),ele_name(1),iions
         ENDIF
         DO Irange=2,ntype
         IF(ICOUNT(Irange-1) < I  .and. I <=  ICOUNT(Irange) )Then
         WRITE(9,'(I3,3X,3I3,3X,A,I)')I,id_map(1:3,J),ele_name(Irange),iions
         ENDIF
         ENDDO
      ENDDO
      NCOUNT_TOT = NCOUNT_TOT + NCOUNT
      enddo

      write(*,'(3X,A)') &
     &"Congratulations: the supercell is generated successfully !"
      close(1)
      close(3)
      close(9)
      END


      SUBROUTINE RECIPS(A,A1,A2,A3,B1,B2,B3)
      USE mod_comp

      IMPLICIT REAL(DP) (A-H,O-Z)
!                                                                      *
!***********************************************************************
!                                                                      *
!   Routine RECIPS sets up the reciprocal lattice vectors B1 ,B2 and . *
!   B3 for a given set of crystallographic vectors A1, A2, A3.         *
!                                                                      *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!                                                                      *
!      A1, A2 and A3 are the three crystallographic vectors of the     *
!        primitive unit cell (given in cartesian coordinates x,y,z)    *
!        of the real space lattice given in units of bohrs ('a').      *
!        A contains the 'lattice constant'.                            *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      B1, B2 and B3 are the three crystallographic vectors of the     *
!        primitive unit cell (given in cartesian coordinates x,y,z) of *
!        the k-space (reciprocal space) lattice in units of '2*pi/a'.  *
!                                                                      *
!***********************************************************************

      DIMENSION A1(3),A2(3),A3(3),B1(3),B2(3),B3(3)

      DEN=0._DP
      I=1
      J=2
      K=3
      S=1._DP
    1 DO 2 IPERM=1,3
         DEN=DEN+S*A1(I)*A2(J)*A3(K)
         L=I
         I=J
         J=K
         K=L
    2 CONTINUE
      I=2
      J=1
      K=3
      S=-S
      IF (S<0._DP) GOTO 1
      I=1
      J=2
      K=3
      DEN=A/ABS(DEN)
      DO 5 IR=1,3
         B1(IR)=DEN*(A2(J)*A3(K)-A2(K)*A3(J))
         B2(IR)=DEN*(A3(J)*A1(K)-A3(K)*A1(J))
         B3(IR)=DEN*(A1(J)*A2(K)-A1(K)*A2(J))
         L=I
         I=J
         J=K
         K=L
    5 CONTINUE
      RETURN
      END

      SUBROUTINE CELVOL(A1,A2,A3,OMEGA)
      USE mod_comp

      IMPLICIT REAL(DP) (A-H,O-Z)
!                                                                      *
!***********************************************************************
!                                                                      *
!   Routine CELVOL calculates the volume of the cell spanned by the    *
!   vectors A1, A2 and A3.                                             *
!                                                                      *
!                                                                      *
!   Input parameters:                                                  *
!   -----------------                                                  *
!                                                                      *
!      A1, A2 and A3 are the three crystallographic vectors of the     *
!        primitive unit cell (given in cartesian coordinates x,y,z)    *
!        of the real space lattice given in units of bohrs ('a').      *
!                                                                      *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!                                                                      *
!      OMEGA is the cell volume. (Warning: Can also be negative !!)    *
!                                                                      *
!***********************************************************************

      DIMENSION A1(3),A2(3),A3(3)

      OMEGA=0._DP
      S=1._DP
      I=1
      J=2
      K=3
    1 DO 2 IPERM=1,3
         OMEGA=OMEGA+S*A1(I)*A2(J)*A3(K)
         L=I
         I=J
         J=K
         K=L
    2 CONTINUE
      I=2
      J=1
      K=3
      S=-S
      IF (S<0._DP) GOTO 1
      RETURN
      END
