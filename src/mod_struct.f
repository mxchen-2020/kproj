       module mod_struct
         use mod_comp 

         real(DP),allocatable, public   :: pos(:,:),posca(:,:)
         real(DP),public                :: ACELL(3,3),BCELL(3,3)
         real(DP),public                :: latt(3,3),relatt(3,3)
         real(DP),public                :: alatt
         integer,public                 :: ntype,nions
         integer,allocatable, public    :: nitype(:),itype(:)

         contains
       
         subroutine poscar
         !-----------------------------------
         ! reads in structural parameters
         ! M. X. Chen -- 10/2012
         !-----------------------------------

           REAL(DP), PARAMETER :: TINY=1.E-10_DP

           real(DP)      :: omega,anorm(3),bnorm(3),V3
           character*40  :: namesys
           character*255 :: line
           character*1   :: tchar
           logical       :: ex,lerr,cartesian

           inquire (file='POSCAR', exist=ex)
           if (ex.eq..false.)then
           write(*,'(A)')'Error: file = POSCAR is not found, ABORTING!'
            stop
           endif           
           open(2,file="POSCAR",status="old")

           ACELL = 0.0_DP
           BCELL = 0.0_DP
           latt  = 0.0_DP
           relatt= 0.0_DP

           read(2,*) namesys
           read(2,*) alatt
           do I=1,3
              read(2,*) ACELL(1:3,I)
           enddo

           ! prepare for coordinate transform
           CALL recip(omega,ACELL,BCELL,anorm,bnorm)

           ! cell information
           if(alatt<0.0_DP)then
             CALL LATTIC(V3)
             alatt = ABS(alatt) / (ABS(V3)**(1._DP/3._DP))
           endif
           latt = alatt*ACELL
!          write(*,*)
!          write(*,'(3X,A)')"Lattice vectors:"
!          do i=1,3
!              write(*,'(5X,3f16.8)')latt(1:3,i)
!          enddo
!          write(*,*)
!          write(*,'(3X,A)')"ACELL and BCELL:"
!          do i=1,3
!              write(*,'(2(5X,3f16.8))')ACELL(1:3,i),BCELL(1:3,i)
!          enddo
           !
           ! calculate volume and reciprocal lattice
           !
           CALL recip(omega,latt,relatt,anorm,bnorm)
           !
!          write(*,*)
!          write(*,'(3X,A,f16.8)')"Volume of cell:",omega
           if(omega<=0.0)then
             write(*,'(3X,A)')"ERROR in lattice cell"
             write(*,'(A)')"ERROR in lattice cell. Stop!!!"
             stop
           endif
!          write(*,*)
           write(*,'(3X,A)')"Reciprocal lattice vectors:"
           do i=1,3
               write(*,'(5X,3f16.8)')relatt(1:3,i)
           enddo
           !
           ! find types of atoms
           !
           read(2,'(A)') line
           read(line,*) tchar
           if (tchar>='0' .and. tchar<='9') then ! if this is a number its vasp4.*
              call num_elements(line,255,ntype)  ! now the number of atom sorts in the POSCAR
           else                                  ! else it is vasp5.* and one line has to be discarded
              read(2,'(A)') line
              call num_elements(line,255,ntype)  ! now the number of atom sorts in the POSCAR
           endif
!          write(*,*)
           write(*,'(3X,A,I5,A)')'Found ',ntype,' atom type(s)'
           allocate(nitype(ntype))
           read(line,*)(nitype(i),i=1,ntype)
!          write(*,'(3X,A)')"Number of atom for each type:"
!          write(*,'(A)')"     Type   No."
!          do i=1,ntype
!             write(*,'(5X,I3,I6)')i,nitype(i)
!          enddo
           nions=0
           do i=1,ntype
              nions=nions+nitype(i)
           enddo
           write(*,'(5X,A,I5)')"Total No. atom:",nions
           allocate(itype(nions))
           do i=1,nitype(1)
              itype(i)=1
           enddo
           if(ntype>=2)then
             m=nitype(1)
             j=1
             do i=2,ntype
                m=m+nitype(i)
                j=j+nitype(i-1)
                do k=j,m
                   itype(k)=i
                enddo
             enddo
           endif
!          write(*,'(A)')"     Atom   Type"
!          do i=1,nions
!             write(*,'(5X,I3,I6)')i,itype(i)
!          enddo
           !
           ! read in postions of atoms
           !
           cartesian = .FALSE.
!          write(*,*)
           read(2,'(A1)') tchar
           if(tchar=='S')then
             read(2,'(A1)') tchar
             if (tchar=='K'.OR.tchar=='k'.OR. &
     &        tchar=='C'.OR.tchar=='c') then
              tchar='K'
              cartesian = .TRUE.
!             write(*,'(3X,A)') &
!    & 'original pos in cartesian coordinates'
             else
              tchar='D'
!             write(*,'(3X,A)')'original pos in direct lattice'
             endif

           elseif(tchar=='K'.OR.tchar=='k'.OR. &
     &        tchar=='C'.OR.tchar=='c') then
              tchar='K'
              cartesian = .TRUE.
!             write(*,'(3X,A)') &
!    & 'original pos in cartesian coordinates'
           else
              tchar='D'
!             write(*,'(3X,A)')'original pos in direct lattice'
           endif

           allocate(pos(3,nions),posca(3,nions))
           do i=1,nions
              read(2,*)pos(1:3,i)
           enddo

           ! transfrom from cartesian to direct coordinates
           if(cartesian) CALL KARDIR(nions,pos,relatt)
           do i=1,nions
              pos(1,i)=MOD(MOD(pos(1,i),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
              pos(2,i)=MOD(MOD(pos(2,i),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
              pos(3,i)=MOD(MOD(pos(3,i),1.0_DP)+60.0_DP+TINY/2,1.0_DP)-TINY/2-0.0_DP
!             write(*,'(3f14.10)')pos(1:3,i)
           enddo
           
           !---> global cartesian coordinates
!          write(*,'(5X,A)')"          direct     coordinates       &
!    &        global cartesian coordinates"
           do i=1,nions
              call dire2car(pos(:,i),posca(:,i),latt)
!             write(*,'(5X,I3,3f12.6,3X,3f12.6)')i,pos(1:3,i),posca(1:3,i)
           enddo

           ! sanity check: atoms must be distinct
           lerr = .false.
           do n=1,nions
             do i=n+1,nions
                if( all( abs( posca(:,n) - posca(:,i) ) < eps ) ) then
!                 write(*,'(/," Error in atomic positions: atoms",i3,&
!    &              " and",i3," are the same")') n,i
                  lerr = .true.
                endif
            enddo
          enddo
          if( lerr ) stop

          write(*,'(3X,A)')"POSCAR OK!"
          
           close(2)
         end subroutine poscar
         !
         !---------------- subroutine recip --------------------------
         !
         !  subroutine for calculating the reciprocal lattice from 
         !  the direct lattice
         !  in addition the norm of the lattice-vectors and the volume 
         !  of the basis-cell is calculated
         !
         !--------------------------------------------------------------
         !
         subroutine recip(omega,A,B,anorm,bnorm)
           real(DP)   ::  A(3,3),B(3,3),anorm(3),bnorm(3),omega

           CALL EXPRO(B(1,1),A(1,2),A(1,3))
           CALL EXPRO(B(1,2),A(1,3),A(1,1))
           CALL EXPRO(B(1,3),A(1,1),A(1,2))

           omega =ABS(B(1,1)*A(1,1)+B(2,1)*A(2,1)+B(3,1)*A(3,1))

           do I=1,3
             do J=1,3
                B(I,J)=B(I,J)/omega
             enddo
           enddo

           do I=1,3
              anorm(I)=SQRT(A(1,I)*A(1,I)+A(2,I)*A(2,I)+A(3,I)*A(3,I))
              bnorm(I)=SQRT(B(1,I)*B(1,I)+B(2,I)*B(2,I)+B(3,I)*B(3,I))
           enddo
         end subroutine recip
         !
         !----------------- SUBROUTINE EXPRO   ---------------------------------
         ! EXPRO
         ! caclulates the x-product of two vectors
         !
         !----------------------------------------------------------------------
         !
         SUBROUTINE EXPRO(H,U1,U2)
!          use nrtype
           real(DP)    :: H(3),U1(3),U2(3)

           H(1)=U1(2)*U2(3)-U1(3)*U2(2)
           H(2)=U1(3)*U2(1)-U1(1)*U2(3)
           H(3)=U1(1)*U2(2)-U1(2)*U2(1)
         END SUBROUTINE EXPRO
         !
         !---------------------------------------------------------------------
         !
         ! The subroutine num_elements(charin,charlen,nel)
         ! counts the number of float elements in the string charin and returns it via nel.
         !
         ! INPUT
         ! charin: (CHARACTER) that is the character string from which the number elements should
         ! be extracted
         ! charlen: (INTEGER) gives the length of the character
         !
         ! OUTPUT
         ! nel: (INTEGER) gives the number of elements
         !
         !---------------------------------------------------------------------
         !
         subroutine num_elements(charin,charl,nel)
         IMPLICIT NONE
         INTEGER :: nel,i,charl
         CHARACTER(charl):: charin
         LOGICAL ::last_int=.FALSE.
         nel =0
         do i=1,charl
            if (charin(i:i)>='0' .and. charin(i:i)<='9' .or. charin(i:i).eq.'.') then
               last_int=.TRUE.
            else if (last_int.eqv..TRUE.) then
               nel=nel+1
               last_int=.FALSE.
           endif
        end do
        end subroutine num_elements
        !
        !-------------------------------------------------------------
        ! direct representation to cartessian representation
        !-------------------------------------------------------------
        !
        subroutine dire2car(A,B,cell)
          real(DP)     :: A(3), B(3), cell(3,3)
          B(1)=A(1)*cell(1,1)+A(2)*cell(1,2)+A(3)*cell(1,3)
          B(2)=A(1)*cell(2,1)+A(2)*cell(2,2)+A(3)*cell(2,3)
          B(3)=A(1)*cell(3,1)+A(2)*cell(3,2)+A(3)*cell(3,3)
        end subroutine dire2car

!**************** SUBROUTINE KARDIR ************************************
! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE KARDIR(NMAX,V,BASIS)
      IMPLICIT REAL(DP) (A-H,O-Z)
      DIMENSION V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE

!**************** SUBROUTINE LATTIC  ***********************************
! RCS:  $Id: lattice.F,v 1.2 2001/04/05 10:34:10 kresse Exp $
!
!  subroutine for calculating the reciprocal lattice from the direct
!  lattice
!  in addition the norm of the lattice-vectors and the volume of
!  the basis-cell is calculated
!***********************************************************************

      SUBROUTINE LATTIC(VOLM)

      IMPLICIT NONE

      REAL(DP) Omega, VOLM, ANORM(3),BNORM(3)
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(BCELL(1:3,1),ACELL(1:3,2),ACELL(1:3,3))
      CALL EXPRO(BCELL(1:3,2),ACELL(1:3,3),ACELL(1:3,1))
      CALL EXPRO(BCELL(1:3,3),ACELL(1:3,1),ACELL(1:3,2))

      Omega =BCELL(1,1)*ACELL(1,1)+BCELL(2,1)*ACELL(2,1) &
     &      +BCELL(3,1)*ACELL(3,1)

!     DO I=1,3
!     DO J=1,3
!       BCELL(I,J)=BCELL(I,J)/Omega
!     ENDDO
!     ENDDO

!     DO I=1,3
!       ANORM(I)=SQRT(SUM(ACELL(:,I)*ACELL(:,I)))
!       BNORM(I)=SQRT(SUM(BCELL(:,I)*BCELL(:,I)))
!     ENDDO
      VOLM = Omega
      RETURN
      END SUBROUTINE

       end module mod_struct
