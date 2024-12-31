      module mod_lattice

      use mod_comp

!---> (supercell) lattice vectors in full units (a.u.)
      real(DP), save :: a1(3),a2(3),a3(3)

!===> lattice vectors, including all lattice constants and 2pi

      ! reciprocal lattice vectors
      real(DP),save :: b1(3),b2(3),b3(3)

!  real space (column) vectors
      real(DP), save :: amat(3,3)
!  reciprocal space (row) vectors
      real(DP), save :: bmat(3,3)

!  metric for real space
      real(DP), save :: aamat(3,3)
!  metric for reciprocal space
      real(DP), save :: bbmat(3,3)

!  unit cell volume in full Cartesian units
      real(DP), save :: omtil

      contains


      subroutine set_grid
!*********************************************************************
!     set up lattice quantities and matrices
!*********************************************************************

      use mod_comp
      use mod_wave, only : gmax,CELL,NGX,NGY,NGZ,exvasp

      implicit none
      real(DP)  :: aa1,aa2,aa3

      ! Angstrom to a.u.
!------------------------------------------------------------    
      IF(exvasp==.true.) THEN    !vasp--  A
      a1(1:3) = CELL(1:3,1)/AUTOA
      a2(1:3) = CELL(1:3,2)/AUTOA
      a3(1:3) = CELL(1:3,3)/AUTOA
      ELSE
      a1(1:3) = CELL(1:3,1)      !qe  -- .a.u
      a2(1:3) = CELL(1:3,2)
      a3(1:3) = CELL(1:3,3)
      END IF
!------------------------------------------------------------
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

!  matrices of lattice vectors in full Cartesian units
!-----------------------------------------
      IF(exvasp==.true.) THEN
      amat(:,:) = CELL(:,:)/AUTOA
      ELSE
      amat(:,:) = CELL(:,:)
      END IF
      bmat(1,:) = TWOPI * b1(:)
      bmat(2,:) = TWOPI * b2(:)
      bmat(3,:) = TWOPI * b3(:)
!------------------------------------------
      aamat = matmul( transpose(amat) , amat )
      bbmat = matmul( bmat, transpose(bmat) )

      aa1 = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
      aa2 = sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 )
      aa3 = sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 )

      NGX = int( aa1*gmax/TWOPI + eps ) + 1
      NGY = int( aa2*gmax/TWOPI + eps ) + 1
      NGZ = int( aa3*gmax/TWOPI + eps ) + 1

      write(IUO,'(5X,A,3I5)')"NGX,NGY,NGZ: ",NGX,NGY,NGZ
      return
      end subroutine set_grid

      subroutine gen_gvect(VK,NG_IN,NG,GV)
        use mod_wave, only : kmax,NGX,NGY,NGZ

        real(DP),intent(in):: VK(3)
        integer,intent(in) :: NG_IN
        integer,intent(out):: NG,GV(3,NG_IN)
        integer            :: j,ig1,ig2,ig3,ig1p,ig2p,ig3p
        real(DP)           :: g(3),g2,kcut2

        kcut2 = kmax**2
        NG = 0

        do ig3=0,2*NGZ
            ig3p=ig3
            if (ig3 > NGZ) ig3p=ig3-2*NGZ-1
            do ig2=0,2*NGY
               ig2p=ig2
               if (ig2 > NGY) ig2p=ig2-2*NGY-1
               do ig1=0,2*NGX
                  ig1p=ig1
                  if (ig1 > NGX) ig1p=ig1-2*NGX-1

                  ! change to Michael's implementation
                  g = [ ig1p+VK(1), ig2p+VK(2), ig3p+VK(3) ]
                  g2 = dot_product( g, matmul( bbmat, g ) )
                  if( g2 > kcut2 ) cycle
                  NG = NG +1
                  if(NG<=NG_IN)then
                    GV(1,NG)=ig1p
                    GV(2,NG)=ig2p
                    GV(3,NG)=ig3p
                  endif
               enddo
            enddo
         enddo
      end subroutine gen_gvect
      end module mod_lattice
