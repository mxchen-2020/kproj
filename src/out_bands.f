      subroutine out_bands
      ! 
      implicit none


!              wavecar file variables
      real(DP),allocatable :: e(:),occ(:)
      real(DP)             :: vkpt(3),rnv
      integer              :: nv

!              k-projected information
      real(DP),allocatable :: proj_wt(:,:,:),eig(:,:,:),bkpts(:,:)

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
      bskfn ='bs_projected.dat'

      IF(phonon_unfold.or.lcao_unfold=.true.) ISPIN=1  

      allocate( e(NBANDS),occ(NBANDS))
      emin_ef = -1000
      emax_ef =  1000
      do jspin=1,ISPIN
        do ikpt=1,NKPTS
           e(1:NBANDS)   = CELTOT(1:NBANDS,ikpt,jspin)
           occ(1:NBANDS) = FERTOT(1:NBANDS,ikpt,jspin) 
!  get a better estimate of the fermi energy
         occ(1:nbands) = occ(1:nbands)/occ(1) ! assumes 1 occupied state
         emin_ef = max( maxval( e, mask = ( occ >= 0.5_fp ) ) , emin_ef )
         emax_ef = min( minval( e, mask = ( occ <  0.5_fp ) ) , emax_ef )
        enddo
      enddo
      ef = ( emin_ef + emax_ef )/2
      write(IUO,'(/,5X,"Fermi level obtained from bands = ",f0.5," eV")') ef
      write(IUO,'("    (may need to change in bs_projected.dat)")')
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
         if( all( abs( dk ) < 1.e-05_fp ) ) then ! vasp sometimes does
            l_good_kpt(ikpt) = .false.           ! k-points multiple times
            cycle
         endif
         dk = dk/maxval( abs( dk ) )
         if( any( abs( dk - dk_old ) > 1.e-04_fp ) ) then
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
            x_kpt(k)=x_kpt(k)/AUTOA
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
      write(bskfh,'(6f12.6)') x_kpt( segment_endpoint(1:kpt_lines) )
      
      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &             kpt,NBANDS,jspin,bkpts(:,ikpt),x_kpt(ikpt)
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                             &
     &             eig(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo 

      if(ISPIN==1)then
        open(66,file="band_proj_gnu.dat",status='unknown')
        do nb=1,NBANDS
          do ikpt=1,NKPTS
            if( .not. l_good_kpt(ikpt) ) cycle
            write(66,'(3f16.8)')x_kpt(ikpt),CELTOT(nb,ikpt,1),proj_wt(nb,ikpt,1)
          enddo
            write(66,*)
        enddo
        close(66)
      else
        open(66,file="band_proj_gnu_up.dat",status='unknown')
        do nb=1,NBANDS
          do ikpt=1,NKPTS
            if( .not. l_good_kpt(ikpt) ) cycle
            write(66,'(3f16.8)')x_kpt(ikpt),CELTOT(nb,ikpt,1),proj_wt(nb,ikpt,1)
          enddo
            write(66,*)
        enddo
        close(66)

        open(66,file="band_proj_gnu_dn.dat",status='unknown')
        do nb=1,NBANDS
          do ikpt=1,NKPTS
            if( .not. l_good_kpt(ikpt) ) cycle
            write(66,'(3f16.8)')x_kpt(ikpt),CELTOT(nb,ikpt,2),proj_wt(nb,ikpt,2)
          enddo
            write(66,*)
        enddo
        close(66)
      endif

      close( bskfh )

      end subroutine out_bands

