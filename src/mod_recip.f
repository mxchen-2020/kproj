      module mod_recip

      use mod_comp
      use mod_struct


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


      subroutine get_recip
!*********************************************************************
!     set up lattice quantities and matrices
!*********************************************************************

      implicit none
      real(DP)  :: aa1,aa2,aa3

      ! Angstrom to a.u.
      a1(1:3) = latt(1:3,1)
      a2(1:3) = latt(1:3,2)
      a3(1:3) = latt(1:3,3)

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
      bmat(1,:) = TWOPI * b1(:)
      bmat(2,:) = TWOPI * b2(:)
      bmat(3,:) = TWOPI * b3(:)

      write(*,*)
      write(*,'(3X,3f12.6)')bmat(1,1:3)/TWOPI
      write(*,'(3X,3f12.6)')bmat(2,1:3)/TWOPI
      write(*,'(3X,3f12.6)')bmat(3,1:3)/TWOPI

      bbmat = matmul( bmat, transpose(bmat) )

      return
      end subroutine get_recip

      end module mod_recip
