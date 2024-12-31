       program combine_bs
        ! sum over weights in seperated bs_projected.dat 
        !(summing over different spatial windows)

        integer, parameter  :: DP =SELECTED_real_KIND(10)
        integer             :: IUIN,neigen,NKPTS,ISPIN,NBANDS
        integer,allocatable :: NK(:)
        character(3)        :: efile
        real(DP),allocatable:: ef(:)
        logical             :: ex
        integer             :: KDUM,SPDUM,ig,ik,ikpt,jspin,nb

        integer,allocatable  :: nseg(:),nbnd(:)
        real(DP),allocatable :: x_kpt_re(:,:)
        real(DP),allocatable :: CELTOT(:,:,:,:),proj_wt(:,:,:,:)
        real(DP),allocatable :: bkpts(:,:,:),x_kpt(:,:)

        real(DP),allocatable :: proj_wt_TOT(:,:,:)

        write(*,'(5X,A)')"No. of bs_prpjected.dat files"
        write(*,'(5X,A)')"I assume that they are bs_prpjected.dat_1, &
     & bs_prpjected.dat_2, ..., bs_prpjected.dat_n"
        read(*,*)neigen
        if(neigen<=1)then
            write(*,'(5X,A)') "Maybe we have only one efile, Aborting!"
            stop
        endif

        IUIN = 1
        NKPTS = 0
        nseg = 0 
        NKPTS = -1

        ! get the total number of KPTS and bands
        allocate(NK(neigen),nseg(neigen),nbnd(neigen),ef(neigen))
        do ig=1,neigen
          call open_eigen(IUIN,ig)
          read(IUIN,*) NK(ig),ISPIN,ef(ig)
          read(IUIN,*) nseg(ig)
          read(IUIN,*) 
          read(IUIN,*)KDUM,nbnd(ig)
          close(IUIN)
        enddo

        do ig=1,neigen-1
           if(NK(ig+1)/=NK(ig))then
            write(*,'(5X,A)')"NK is not consistent! stoping"
            stop
           endif
           if(nseg(ig+1)/=nseg(ig))then
            write(*,'(5X,A)')"nseg is not consistent! stoping"
            stop
           endif
           if(nbnd(ig+1)/=nbnd(ig))then
            write(*,'(5X,A)')"nbnd is not consistent! stoping"
            stop
           endif
        enddo
        allocate(x_kpt_re(nseg(1),neigen))
        x_kpt_re = 0.0_DP

        NKPTS = NK(1)
        NBANDS = nbnd(1)

        allocate(CELTOT(NBANDS,NKPTS,ISPIN,neigen))
        allocate(proj_wt(NBANDS,NKPTS,ISPIN,neigen))
        allocate(proj_wt_TOT(NBANDS,NKPTS,ISPIN))
        allocate(bkpts(3,NKPTS,neigen),x_kpt(NKPTS,neigen))

        CELTOT = 0.0_DP
        proj_wt= 0.0_DP
        proj_wt_TOT= 0.0_DP
        bkpts  = 0.0_DP
        x_kpt  = 0.0_DP
      
        ! get eigenvalues, and weights 
        do ig=1,neigen
          call open_eigen(IUIN,ig)
          read(IUIN,*) NK(ig),ISPIN,ef(ig)
          read(IUIN,*) nseg(ig)
          read(IUIN,*) (x_kpt_re(i,ig),i=1,nseg(ig))

          do jspin = 1,ISPIN
            do ikpt = 1,NK(ig)
              read(IUIN,*)KDUM,BDUM,SPDUM,bkpts(1:3,ikpt,ig), &
     & x_kpt(ikpt,ig)
              do nb = 1,nbnd(ig)
                 read(IUIN,*)CELTOT(nb,ikpt,jspin,ig), &
     & proj_wt(nb,ikpt,jspin,ig)
              enddo
            enddo
          enddo
          close(IUIN)
        enddo

        ! check consistence
        do ig= 1,neigen-1
          do jspin = 1,ISPIN
            do ikpt = 1,NK(1)

                if(abs(x_kpt(ikpt,ig+1)-x_kpt(ikpt,ig))>1E-05)then
                  write(*,'(5X,A)')"Weights are inconsistent, stoping!" 
                  stop
                endif

                if(abs(bkpts(1,ikpt,ig+1)-bkpts(1,ikpt,ig))>1E-05.or.&
     & abs(bkpts(2,ikpt,ig+1)-bkpts(2,ikpt,ig))>1E-05.or. &
     & abs(bkpts(3,ikpt,ig+1)-bkpts(3,ikpt,ig))>1E-05)then
                  write(*,'(5X,A)')"K-points are inconsistent, stoping!"
                  stop
                endif

              do nb = 1,nbnd(1)
                if(abs(CELTOT(nb,ikpt,jspin,ig+1)- &
     & CELTOT(nb,ikpt,jspin,ig))>1E-06)then
                  write(*,'(5X,A)')"Bands are inconsistent, stoping!"
                  stop
                endif
              enddo
            enddo
          enddo
        enddo

        do jspin = 1,ISPIN
          do ikpt = 1,NK(1)
            do nb = 1,nbnd(1)
              do ig= 1,neigen
                proj_wt_tot(nb,ikpt,jspin) =  &
     & proj_wt_tot(nb,ikpt,jspin) + proj_wt(nb,ikpt,jspin,ig)
              enddo
            enddo
          enddo
        enddo

!--->    write to bs_projected.dat

      open(2,file="bs_projected.dat",status='unknown')
      write(2,'(i0,2x,i0,f20.8,"  ! NKPTS, jspins, ef (eV)")')      &
     &       NKPTS,ISPIN,ef(1)
      write(2,'(i0,"  ! number of segments; positions:")') nseg(1)
      write(2,'(6f12.6)') x_kpt_re(1:nseg(1),1)

      do jspin = 1,ISPIN
        do ikpt = 1,NKPTS

          write(2,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &             ikpt,NBANDS,jspin,bkpts(1:3,ikpt,1),x_kpt(ikpt,1)
          do nb = 1,NBANDS
            write(2,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin,1),proj_wt_tot(nb,ikpt,jspin)
          enddo
        enddo
      enddo

      close(2)    

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


       end
