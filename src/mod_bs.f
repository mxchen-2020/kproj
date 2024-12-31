       module mod_bs

        ! read eigeivalues and weights from bs_projected.dat

        USE mod_comp

        integer,public              :: ISPIN,NKPTS,NBANDS
        real(DP),allocatable,public :: CELTOT(:,:,:),proj_wt(:,:,:)
        real(DP),allocatable,public :: proj_wt_BS(:,:,:)
        real(DP),allocatable,public :: bkpts(:,:),x_kpt(:)
        real(DP),allocatable,public :: x_kpt_re(:)
        integer,public              :: nseg_tot

        real(DP),public   :: EFERMI

        contains

        subroutine read_bs
        logical         :: ex,LFLIP
        integer         :: KDUM,BDUM,SPDUM,ig,ik,ikpt,jspin,nb,isg,IUIN
        integer         :: nseg

!       integer              :: nseg
!       real(DP),allocatable :: x_kpt_re(:)
!       real(DP),allocatable ::CELTOT_TMP(:,:,:),proj_wt_TMP(:,:,:)
!       real(DP),allocatable ::bkpts(:,:),x_kpt(:)


!=============================================================================
!       read bs_projected.dat to 
!       get NKPTS,ISPIN,EFERMI, NBANDS, eigenvalues and k-projected weights

        IUIN = 4

        ! open bs_projected.dat 
        inquire (file="bs_projected.dat", exist=ex)
        if (ex.eq..false.) then
          write(*,'(5X,A,I,A)')                                         &
     & "Error: file = bs_projected.dat",ig,                            &
     & " is not found, ABORTING!"
          stop
        endif
        open (IUIN,FILE="bs_projected.dat",status='old')

        ! get the total number of KPTS and bands
        read(IUIN,*) NKPTS,ISPIN,EFERMI
        read(IUIN,*) nseg

        nseg_tot = nseg 
        allocate(x_kpt_re(nseg)) 
        read(IUIN,*)(x_kpt_re(isg),isg=1,nseg)

        allocate(bkpts(3,NKPTS),x_kpt(NKPTS))

        read(IUIN,*)KDUM,NBANDS,SPDUM,bkpts(1:3,1),x_kpt(1)
        write(*,'(3X,A,3I5)')"ISPIN,NKPTS,NBANDS:",ISPIN,NKPTS,NBANDS
!       write(*,'(3X,A,I5)')"nseg_tot:",nseg_tot
!       write(*,*)(x_kpt_re(isg),isg=1,nseg)

        allocate(CELTOT(NBANDS,NKPTS,ISPIN),proj_wt(NBANDS,NKPTS,ISPIN))
        allocate(proj_wt_BS(NBANDS,NKPTS,ISPIN))

        ! get eigenvalues, and weights 
!       rewind(IUIN)
!       read(IUIN,*) NKPTS,ISPIN,EFERMI
!       read(IUIN,*) nseg
!       read(IUIN,*)(x_kpt_re(isg),isg=1,nseg)

        backspace(IUIN)
         
        do jspin = 1,ISPIN
          do ik = 1,NKPTS
              read(IUIN,*)KDUM,BDUM,SPDUM,bkpts(1:3,ik),x_kpt(ik)
!             write(*,*)KDUM,BDUM,SPDUM,bkpts(1:3,ik),x_kpt(ik)
              do nb = 1,NBANDS
            read(IUIN,*)CELTOT(nb,ik,jspin),proj_wt_BS(nb,ik,jspin)
!           write(*,*)CELTOT_BS(nb,ik,jspin),proj_wt_BS(nb,ik,jspin)
              enddo
          enddo
        enddo
        close(IUIN)

        write(*,'(3X,A)')"bs_projected.dat OK!"
        end subroutine 

       end module mod_bs 
