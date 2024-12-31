       program collect_eig
        ! put eigenvalues and weights in seperated bs_projected.dat together

        integer             :: mode

        integer, parameter  :: DP =SELECTED_real_KIND(10)
        integer             :: IUIN,neigen,NKPTS,NBANDS,NK0,ISPIN,NBDIM
        integer,allocatable :: NK(:)
        character(3)        :: efile
        double precision    :: ef,EFERMI
        logical             :: ex,LFLIP
        integer             :: KDUM,SPDUM,ig,ik,ikpt,jspin,ib,nb

        integer              :: nseg,nflip,NKMAX
        integer,allocatable  :: nseg0(:),nbnd(:)
        real(DP),allocatable :: x_kpt_re(:),x_lines(:)
        real(DP),allocatable :: CELTOT(:,:,:),proj_wt(:,:,:)
        real(DP),allocatable :: CELTOT_TMP(:,:,:,:),proj_wt_TMP(:,:,:,:)
        real(DP),allocatable :: CELTOT_INV(:,:,:),proj_wt_INV(:,:,:)
        real(DP),allocatable :: bkpts_INV(:,:),x_kpt_INV(:)
        real(DP),allocatable :: bkpts(:,:),x_kpt(:)
        real(DP),allocatable :: bkpts_TMP(:,:,:),x_kpt_TMP(:,:)

        real(DP),allocatable ::CELTOT_TOT(:,:,:),proj_wt_TOT(:,:,:)
        real(DP),allocatable ::bkpts_TOT(:,:),x_kpt_TOT(:)

        LFLIP=.false.
        write(*,'(5X,A)')"No. of eigenvalue files"
        read(*,*)neigen
        !if(neigen<=1)then
        !    write(*,'(5X,A)') "Maybe we have only one efile, Aborting!"
        !    stop
        !endif
 
        mode = 1
        write(*,'(5X,A)')"Mode 1: merge k-points"
        write(*,'(5X,A)') &
     & "Mode 2: merge bands (kpts for each bs file must be the same)"
        read(*,*)mode
        if(mode/=2)mode = 1

        nflip=0
        write(*,'(5X,A)')"Flip bands around the first kpt?"
        write(*,'(5X,A)')"0 for No and 1 for yes"
        read(*,*)nflip
        if(nflip==1)LFLIP=.true.

        IUIN = 1
        NKPTS = 0
        nseg = 0 
        NKMAX = -1


        ! get the total number of KPTS and bands
        allocate(nseg0(neigen),nbnd(neigen))
        do ig=1,neigen
          call open_eigen(IUIN,ig)
          read(IUIN,*) NK0,ISPIN,EFERMI
          read(IUIN,*) nseg0(ig)
          
          allocate(x_lines(nseg0(ig)))
          read(IUIN,*) (x_lines(i),i=1,nseg0(ig))
          read(IUIN,*)KDUM,nbnd(ig)
          deallocate(x_lines)
          close(IUIN)

          if(mode==1)then
            NKPTS = NKPTS + NK0
            nseg = nseg + nseg0(ig)
          else
            NKPTS = NK0
            nseg =  nseg0(ig)
          endif

          NKMAX = max(NKMAX,NK0)
        enddo

        if(mode==1)then
           NBANDS=minval(nbnd)
           NBDIM=maxval(nbnd)
        else
          do ig=1,neigen
             NBANDS = NBANDS + nbnd(ig)
          enddo
          NBDIM = nbnd(1)
        endif

        write(*,'(5X,A,3I7)')"ISPIN,NKPTS,NBANDS:",ISPIN,NKPTS,NBANDS

        allocate(bkpts(3,NKPTS),x_kpt(NKPTS))
        allocate(CELTOT(NBANDS,NKPTS,ISPIN))
        allocate(proj_wt(NBANDS,NKPTS,ISPIN))
        allocate(CELTOT_TMP(NBDIM,NKMAX,ISPIN,neigen))
        allocate(proj_wt_TMP(NBDIM,NKMAX,ISPIN,neigen))
        allocate(bkpts_TMP(3,NKMAX,neigen),x_kpt_TMP(NKMAX,neigen))
        allocate(x_kpt_re(nseg)) 
      
        allocate(NK(neigen))

        ! get eigenvalues, and weights 
        ef = 1000
        nseg=0
        
        do ig=1,neigen
          call open_eigen(IUIN,ig)
          read(IUIN,*) NK(ig),ISPIN,EFERMI
          ef = min(EFERMI,ef)
          read(IUIN,*) nseg0(ig)
          
          if(mode==1)then
            read(IUIN,*) (x_kpt_re(i),i=nseg+1,nseg+nseg0(ig))
            nseg = nseg + nseg0(ig)
          else
            read(IUIN,*) (x_kpt_re(i),i=1,nseg0(ig))
            nseg = nseg0(1)
          endif

          do jspin = 1,ISPIN
            do ikpt = 1,NK(ig)
              read(IUIN,*)KDUM,BDUM,SPDUM,bkpts_TMP(1:3,ikpt,ig), &
     & x_kpt_TMP(ikpt,ig)
              do nb = 1,nbnd(ig)
                 read(IUIN,*)CELTOT_TMP(nb,ikpt,jspin,ig), &
     & proj_wt_TMP(nb,ikpt,jspin,ig)
              enddo
            enddo
          enddo
          close(IUIN)
        enddo

        select case(mode)
        case(1)

        ! put all eigenvalues and weights into an array
        do jspin = 1,ISPIN
           ik=0
         do ig=1,neigen
          do ikpt=1,NK(ig)
            ik=ik+1
            x_kpt(ik)=x_kpt_TMP(ikpt,ig)
            bkpts(1:3,ik)=bkpts_TMP(1:3,ikpt,ig)
            do nb = 1,NBANDS
              CELTOT(nb,ik,jspin)=CELTOT_TMP(nb,ikpt,jspin,ig)
              proj_wt(nb,ik,jspin)=proj_wt_TMP(nb,ikpt,jspin,ig)
            enddo
          enddo
         enddo
       enddo

       case(2)

        ! put all eigenvalues and weights into an array
        x_kpt(:)=x_kpt_TMP(:,1)
        bkpts(:,:)=bkpts_TMP(:,:,1)

        do jspin = 1,ISPIN
           ib=0
         do ig=1,neigen
          do nb=1,nbnd(ig)
            ib=ib+1
            do ikpt = 1,NKPTS
              CELTOT(ib,ikpt,jspin)=CELTOT_TMP(nb,ikpt,jspin,ig)
              proj_wt(ib,ikpt,jspin)=proj_wt_TMP(nb,ikpt,jspin,ig)
            enddo
          enddo
         enddo
       enddo

       do jspin = 1,ISPIN
         do ikpt = 1, NKPTS
            call order_band(CELTOT(:,ikpt,jspin),proj_wt(:,ikpt,jspin),NBANDS)
         enddo
       enddo

       end select
    
        ! remove repeated kpoints and out bands
        call out_bands(ISPIN,NKPTS,NBANDS,ef,CELTOT,bkpts,proj_wt,x_kpt)
        !call out_bands(ISPIN,NKPTS,NBANDS,ef,CELTOT,bkpts,proj_wt)

        deallocate(x_kpt_re,CELTOT,bkpts,proj_wt,x_kpt)

       ! if flip, k_1 is the reference point 
       if(LFLIP)then
          open(1,file='bs_projected.dat',status='old')
          read(IUIN,*) NKPTS,ISPIN,ef
          read(IUIN,*) nseg
          allocate(x_kpt_re(nseg))
          read(IUIN,*) (x_kpt_re(i),i=1,nseg)

        allocate(bkpts(3,NKPTS),x_kpt(NKPTS))
        allocate(CELTOT(NBANDS,NKPTS,ISPIN),proj_wt(NBANDS,NKPTS,ISPIN))

          do jspin = 1,ISPIN
            do ikpt = 1,NKPTS
              ik = ik + 1
              read(IUIN,*)KDUM,BDUM,SPDUM,bkpts(1:3,ikpt),x_kpt(ikpt)
              do nb = 1,NBANDS
            read(IUIN,*)CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
              enddo
            enddo
          enddo
          
          allocate(CELTOT_INV(NBDIM,NKPTS-1,ISPIN))
          allocate(proj_wt_INV(NBDIM,NKPTS-1,ISPIN))
          allocate(bkpts_INV(3,NKPTS-1),x_kpt_INV(NKPTS-1))

          allocate(CELTOT_TOT(NBANDS,2*NKPTS-1,ISPIN))
          allocate(proj_wt_TOT(NBANDS,2*NKPTS-1,ISPIN))
          allocate(bkpts_TOT(3,2*NKPTS-1),x_kpt_TOT(2*NKPTS-1))

        do jspin = 1,ISPIN
         do ikpt = 1,NKPTS-1
           do nb = 1,NBANDS
              CELTOT_INV(nb,NKPTS-ikpt,jspin)=CELTOT(nb,ikpt+1,jspin)
              proj_wt_INV(nb,NKPTS-ikpt,jspin)=proj_wt(nb,ikpt+1,jspin)
              x_kpt_INV(NKPTS-ikpt) = -x_kpt(ikpt+1)
              bkpts_INV(:,NKPTS-ikpt) = bkpts(:,ikpt+1)-bkpts(:,1)
           enddo
         enddo
       enddo

       do jspin = 1,ISPIN
         do ikpt = 1,NKPTS-1
           do nb = 1,NBANDS
              CELTOT_TOT(nb,ikpt,jspin)=CELTOT_INV(nb,ikpt,jspin)
              proj_wt_TOT(nb,ikpt,jspin)=proj_wt_INV(nb,ikpt,jspin)
              x_kpt_TOT(ikpt) = x_kpt_INV(ikpt)
              bkpts_TOT(:,ikpt) = bkpts_INV(:,ikpt)
           enddo
         enddo

         do ikpt = 1,NKPTS
           do nb = 1,NBANDS
              CELTOT_TOT(nb,ikpt-1+NKPTS,jspin)=CELTOT(nb,ikpt,jspin)
              proj_wt_TOT(nb,ikpt-1+NKPTS,jspin)=proj_wt(nb,ikpt,jspin)
              x_kpt_TOT(ikpt-1+NKPTS) = x_kpt(ikpt)
              bkpts_TOT(:,ikpt-1+NKPTS) = bkpts(:,ikpt)
           enddo
         enddo
       enddo

       NKPTS=2*NKPTS-1

       do ikpt = 2,NKPTS
          x_kpt_TOT(ikpt)=x_kpt_TOT(ikpt)-x_kpt_TOT(1)
          !write(*,'(I5,f12.6)')ikpt,x_kpt_TOT(ikpt)
       enddo
       x_kpt_TOT(1)=0.0
       rewind(1)

       call out_bands(ISPIN,NKPTS,NBANDS,ef,CELTOT_TOT,                 &
     & bkpts_TOT,proj_wt_TOT,x_kpt_TOT)
       endif

       contains
       SUBROUTINE open_eigen(inf,ig)
        integer,intent(in):: inf,ig
        character(3)      :: efile
        logical           :: ex

        if(ig<10)then
          write(efile,'(i1)')ig
        else if(ig>=10.and.ig<100)then
          write(efile,'(i2)')ig
        else
          write(efile,'(i3)')ig
        endif
        inquire (file="bs_projected.dat_"//efile, exist=ex)
        if (ex.eq..false.) then
          write(*,'(5X,A,I,A)')                                         &
     & "Error: file = bs_projected.dat_",ig,                            &
     & " is not found, ABORTING!"
          stop
        endif
        open (unit=inf,FILE="bs_projected.dat_"//efile,status='old')
       END SUBROUTINE open_eigen

      subroutine order_band(enband,wt,NBANDS)

        integer              :: NBANDS
        real(DP)             :: enband(NBANDS),wt(NBANDS)
        real(DP)             :: energex
        integer              :: l,p

        do l=1,NBANDS-1
          do p=l+1,NBANDS
            if(enband(l)>enband(p))then
              energex = enband(p)
              enband(p) = enband(l)
              enband(l) = energex

              energex = wt(p)
              wt(p) = wt(l)
              wt(l) = energex
            else
            endif
          enddo
        enddo
      end subroutine

      subroutine out_bands(ISPIN,NKPTS,NBANDS,ef,CELTOT,bkpts,proj_wt,x_kpt)

      integer,intent(in)   :: ISPIN,NKPTS,NBANDS
      real(DP),intent(in)  :: CELTOT(NBANDS,NKPTS,ISPIN),bkpts(3,NKPTS)
      real(DP),intent(in)  :: proj_wt(NBANDS,NKPTS,ISPIN)
      real(DP),intent(in)  :: ef,x_kpt(NKPTS)

      real(DP)             :: vkpt(3),rnv
      integer              :: nv

      logical,allocatable  :: l_good_kpt(:)

      integer              :: irec
      integer              :: i,k,l,p,n,nb,jspin,ikpt,ios,kpt,nkpts0

!     file name and handle for band data
      integer,   parameter :: bskfh = 73
      character(len=30)    :: bskfn = 'bs_projected.dat'

      character(len=256)   :: line

     
      integer, allocatable :: segment_endpoint(:),epts(:)
      real(DP)             :: dk(3),dk_old(3),x0,dk_1(3),dk_2(3)
      integer              :: kpt_lines,nseg,k0,max_segments
      real(DP),allocatable :: seg(:),x_kpt_tmp(:)
 
      real(DP)             :: energex
      real(DP),allocatable :: enband_all_spin(:),wt_all_spin(:)

!--->    break up the k-points into lines, assuming that the order is correct
      allocate( l_good_kpt(NKPTS) )
      l_good_kpt = .true.
      max_segments = 2000
      allocate( segment_endpoint(max_segments) )
      nseg = 0
      dk_old = bkpts(:,2) - bkpts(:,1)
      dk_old = dk_old/maxval( abs( dk_old ) ) ! normalize the direction
      do ikpt = 2,NKPTS-1
         dk = bkpts(:,ikpt+1) - bkpts(:,ikpt)
         if( all( abs( dk ) < 1.e-04_DP ) ) then ! vasp sometimes does
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
       write(IUO,'(/,5X,i0," repeated kpts; now total of ",i0/)')       &
     & NKPTS-nkpts0,nkpts0
      endif

      ! get x_kpt_tmp for each k-point
      allocate(x_kpt_tmp(NKPTS))
      x_kpt_tmp = 0.0 
      do ikpt = 2,NKPTS
         dk = bkpts(:,ikpt) - bkpts(:,ikpt-1)

         if( all( abs( dk ) < 1.e-05_DP ) ) then
           ! the two neighboring k-points are the same
           x_kpt_tmp(ikpt) = x_kpt_tmp(ikpt-1)

           
         else
           ! the two k-points are different but may in the same bs file in the correct sequence
           if(abs(x_kpt(ikpt))>abs(x_kpt(ikpt-1)))then
             x0 = x_kpt(ikpt) - x_kpt(ikpt-1)
             x_kpt_tmp(ikpt) = x_kpt_tmp(ikpt-1) + x0
           else
 
             ! in this case the two k-points are different, too
             ! but ikpt is in bs_n and ikpt-1 is in bs_n-1, 
             ! becasue x_kpt(ikpt) stored in bs_n is 0 
             ! usually 2 < n <= NKPTS and
             ! I use dk to judge if it is along the direction in bs_n-1
             ! or the direction in bs_n
             dk_1 = bkpts(:,ikpt-1) - bkpts(:,ikpt-2)
             dk_2 = bkpts(:,ikpt+1) - bkpts(:,ikpt)
 
             dk_1 = dk_1 - dk
             dk_2 = dk_2 - dk
           
             if( all( abs( dk_1 ) < 1.e-05_DP )) then
               x0 = x_kpt(ikpt-1) - x_kpt(ikpt-2)
               x_kpt_tmp(ikpt) = x_kpt_tmp(ikpt-1) + x0

             elseif(all( abs( dk_2 ) < 1.e-05_DP ))then
               x0 = x_kpt(ikpt+1) - x_kpt(ikpt)
               x_kpt_tmp(ikpt) = x_kpt_tmp(ikpt-1) + x0
           
             else
               write(*,'(3X,A)')"ERROR in treating dk, stopping"
               stop
             endif

           endif

         endif  
        
      enddo

!--->    write to bs_projected.dat

      open(bskfh,file=bskfn,form='formatted',status='replace')
      write(bskfh,'(i0,2x,i0,f20.8,"  ! NKPTS, jspins, ef (eV)")')      &
     &       nkpts0,ISPIN,ef
      write(bskfh,'(i0,"  ! number of segments; positions:")') kpt_lines
      write(bskfh,'(6f12.6)') x_kpt_tmp( segment_endpoint(1:kpt_lines) ) 
      
      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &             kpt,NBANDS,jspin,bkpts(:,ikpt),x_kpt_tmp(ikpt)
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo
      close( bskfh )

      if(ISPIN>1)then

      open(33,file="bs_projected.dat_all_spin", &
     & form='formatted',status='unknown') 

      write(33,'(i0,2x,i0,f20.8,"  ! NKPTS, jspins, ef (eV)")')      &
     &       nkpts0,1,ef
      write(33,'(i0,"  ! number of segments; positions:")') kpt_lines
      write(33,'(6f12.6)') x_kpt_tmp( segment_endpoint(1:kpt_lines) ) 

        allocate(enband_all_spin(NBANDS*2),wt_all_spin(NBANDS*2))

        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(33,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &             kpt,NBANDS*2,1,bkpts(:,ikpt),x_kpt_tmp(ikpt)
 
          enband_all_spin = 0.0_DP
          wt_all_spin = 0.0_DP

          ib = 0
          do jspin=1,ISPIN
            do nb=1,NBANDS
               ib = ib + 1
               enband_all_spin(ib) = CELTOT(nb,ikpt,jspin)
               wt_all_spin(ib) = proj_wt(nb,ikpt,jspin)
            enddo
          enddo

          do l=1,NBANDS*2-1
            do p=l+1,NBANDS*2
              if(enband_all_spin(l)>enband_all_spin(p))then
                energex = enband_all_spin(p)
                enband_all_spin(p) = enband_all_spin(l)
                enband_all_spin(l) = energex

                energex = wt_all_spin(p)
                wt_all_spin(p) = wt_all_spin(l)
                wt_all_spin(l) = energex
              else
              endif
            enddo
          enddo

          do nb = 1,NBANDS*2
            write(33,'(f16.10,es22.5)')                              &
     &            enband_all_spin(nb),wt_all_spin(nb)
          enddo
        enddo

      close( 33 )
      endif 


      end subroutine out_bands

       end
