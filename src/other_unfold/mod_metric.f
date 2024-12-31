      module mod_metric

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


      subroutine metric
!*********************************************************************
!     set up lattice quantities and matrices
!*********************************************************************

      use mod_comp
      use mod_input
      implicit none
      real(DP)  :: aa1,aa2,aa3

!------------------------------------------------------------    
      a1(1:3) = CELL(1,:)
      a2(1:3) = CELL(2,:)
      a3(1:3) = CELL(3,:)
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
      amat(:,:) = CELL(:,:)
     
      bmat(1,:) =  b1(:)
      bmat(2,:) =  b2(:)
      bmat(3,:) =  b3(:)
!------------------------------------------
      aamat = matmul( transpose(amat) , amat )
      bbmat = matmul( bmat, transpose(bmat) )*TWOPI

      aa1 = sqrt( a1(1)**2 + a1(2)**2 + a1(3)**2 )
      aa2 = sqrt( a2(1)**2 + a2(2)**2 + a2(3)**2 )
      aa3 = sqrt( a3(1)**2 + a3(2)**2 + a3(3)**2 )

     end subroutine

      end module mod_metric
