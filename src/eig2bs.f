       program eig2bs
        use mod_comp
        ! put eigenvalues in EIGENVAL into the format of bs_band.dat 
        ! POSCAR and EIGENVAL are needed

          real(DP)              :: CELL(3,3),bbmat(3,3)
          real(DP),allocatable  :: proj_wt(:,:,:)
          real(DP),allocatable  :: CELTOT(:,:,:),KPTVEC(:,:)
          integer               :: ISPIN,NKPTS,NBANDS
        
          call initial(ISPIN,NKPTS,NBANDS)

          allocate(proj_wt(NBANDS,NKPTS,ISPIN))
          proj_wt(:,:,:) = 1.0_DP

          allocate(CELTOT(NBANDS,NKPTS,ISPIN))
          allocate(KPTVEC(3,NKPTS))

          call read_eig(ISPIN,NKPTS,NBANDS,CELTOT,KPTVEC)

          call read_poscar_header(CELL)

          call set_lattice(CELL,bbmat)

          call out_bands(ISPIN,NKPTS,NBANDS,CELTOT,KPTVEC,bbmat,proj_wt)
        end

      subroutine initial(ISPIN,NKPTS,NBANDS)
        integer       :: ISPIN,NKPTS,NBANDS
        integer       :: ADUM,BDUM,CDUM,NELECT,i
        logical       :: ex

        inquire (file='EIGENVAL', exist=ex)
        if (ex.eq..false.)then
          write(*,'(A)')'Error: file = EIGENVAL is not found, ABORTING!'
          stop
        endif
        open(1,file="EIGENVAL",status='old')
        read(1,*)ADUM,BDUM,CDUM,ISPIN
        do i=1,4
           read(1,*)
        enddo
        read(1,*)NELECT,NKPTS,NBANDS
        close(1)

      end subroutine initial

      subroutine read_eig(ISPIN,NKPTS,NBANDS,CELTOT,KPTVEC)
        use mod_comp

        integer       :: ISPIN,NKPTS,NBANDS
        real(DP)      :: CELTOT(NBANDS,NKPTS,ISPIN),KPTVEC(3,NKPTS) 

        character*5   :: ADUM
        integer       :: BDUM,CDUM,NELECT,i,IK,IB,ISP

        open(1,file="EIGENVAL",status='old')
        do i=1,6
           read(1,*)
        enddo

        do IK=1,NKPTS 
           read(1,*)
           read(1,*)KPTVEC(1:3,IK)
           do IB=1,NBANDS
              read(1,*)ADUM,(CELTOT(IB,IK,ISP),ISP=1,ISPIN)
           enddo
        enddo

        close(1)

      end subroutine read_eig

      subroutine read_poscar_header(CELL)
        use mod_comp

        real(DP)      :: CELL(3,3)
        character*40  :: title
        logical       :: ex
        real(DP)      :: V3,alatt
        integer       :: I

        inquire (file='POSCAR', exist=ex)
        if (ex.eq..false.)then
          write(*,'(A)')'Error: file = POSCAR is not found, ABORTING!'
          stop
        endif
        open(1,file="POSCAR",status='old')    

        read(1,*)title
        read(1,*)alatt
        do I=1,3
          read(1,*) CELL(1:3,I)
        enddo
        if(alatt<0.0_DP)then
           CALL LATTIC(CELL,V3)
           alatt = ABS(alatt) / (ABS(V3)**(1._DP/3._DP))
        endif
        CELL(:,:)=CELL(:,:)*alatt

      end subroutine read_poscar_header

      subroutine set_lattice(ACELL,bbmat)

      use mod_comp

      implicit none
      real(DP)  :: ACELL(3,3),bbmat(3,3)
      real(DP)  :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
      real(DP)  :: bmat(3,3),omtil

      a1(1:3)=ACELL(1:3,1); a2(1:3)=ACELL(1:3,2); a3(1:3)=ACELL(1:3,3)
      a1=a1/AUTOA; a2=a2/AUTOA; a3=a3/AUTOA

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
      bbmat = matmul( bmat, transpose(bmat) )

      return
      end subroutine set_lattice


      SUBROUTINE EXPRO(H,U1,U2)
        use mod_comp

        real(DP)    :: H(3),U1(3),U2(3)

        H(1)=U1(2)*U2(3)-U1(3)*U2(2)
        H(2)=U1(3)*U2(1)-U1(1)*U2(3)
        H(3)=U1(1)*U2(2)-U1(2)*U2(1)
      END SUBROUTINE EXPRO

      SUBROUTINE LATTIC(ACELL,VOLM)
      use mod_comp

      IMPLICIT NONE

      real(DP)  :: ACELL(3,3), VOLM
      real(DP)  :: BCELL(3,3)
      REAL(DP) Omega,  ANORM(3),BNORM(3)
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(BCELL(1:3,1),ACELL(1:3,2),ACELL(1:3,3))
      CALL EXPRO(BCELL(1:3,2),ACELL(1:3,3),ACELL(1:3,1))
      CALL EXPRO(BCELL(1:3,3),ACELL(1:3,1),ACELL(1:3,2))

      Omega =BCELL(1,1)*ACELL(1,1)+BCELL(2,1)*ACELL(2,1) &
     &      +BCELL(3,1)*ACELL(3,1)

      VOLM = Omega
      RETURN
      END SUBROUTINE

     subroutine out_bands(ISPIN,NKPTS,NBANDS,CELTOT,bkpts,bbmat,proj_wt)
      use mod_comp

      implicit none
 
      integer              :: ISPIN,NKPTS,NBANDS
      real(DP)             :: CELTOT(NBANDS,NKPTS,ISPIN),bkpts(3,NKPTS)
      real(DP)             :: proj_wt(NBANDS,NKPTS,ISPIN),bbmat(3,3)

!              wavecar file variables
      real(DP)             :: vkpt(3),rnv
      integer              :: nv

      real(DP)             :: ef

      logical,allocatable  :: l_good_kpt(:)

      integer              :: irec
      integer              :: i,k,l,n,nb,jspin,ikpt,ios,kpt,nkpts0

!     file name and handle for band data
      integer,   parameter :: bskfh = 73
      character(len=30)    :: bskfn = 'bs_projected.dat'

      character(len=256)   :: line

      character(len=2)     :: ADUM

      integer, allocatable :: segment_endpoint(:),epts(:)
      real(DP),allocatable :: x_kpt(:)
      real(DP)             :: dk(3),dk_old(3),x0
      integer              :: kpt_lines,nseg,k0,max_segments
      
      ! Fermi level remains unshifted, because it does not appear in EIGENVAL
      ef = 0.0_DP 

!--->    break up the k-points into lines, assuming that the order is correct
      allocate( l_good_kpt(NKPTS) )
      l_good_kpt = .true.
      max_segments = 80000
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

      write(*,'(3X,A)')"skip repeated k-points? (Y/N)"
      read(*,*)ADUM
      if(ADUM=="N".or.ADUM=="n".or.ADUM=="NO".or.ADUM=="no".or.&
     & ADUM=="No".or.ADUM=="nO")then
        l_good_kpt = .true.
      endif

      nkpts0 = count( l_good_kpt )
      if( nkpts0 /= NKPTS ) then
       write(*,'(/,5X,i0," repeated kpts; now total of ",i0/)') NKPTS-nkpts0,nkpts0
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
      write(bskfh,'(6f12.6)') &
     & x_kpt( segment_endpoint(1:kpt_lines) )/AUTOA
      
      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &             kpt,NBANDS,jspin,bkpts(:,ikpt),x_kpt(ikpt)/AUTOA
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                             &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo 

      close( bskfh )

      end subroutine out_bands

