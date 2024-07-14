      module mod_wave
        USE mod_comp
        USE mod_input
        USE read_information_from_qe
        USE read_information_from_ab
        integer,public             :: NRSPINORS
        integer,public             :: ISPIN,NKPTS,NBANDS
        real(DP),allocatable,public:: CELTOT(:,:,:),FERTOT(:,:,:)!E,OCC
        real(DP),allocatable,public:: KPTVEC(:,:)
        real(DP),allocatable,public:: NORM_FACTOR(:,:,:)
        logical,public             :: SINGLE_PREC
        real(DP),public            :: kmax,gmax,ENMAX
        integer                    :: NGX,NGY,NGZ,IRECLW
        logical,public ::lsda,noncolin,spinorbit
        logical,public ::exvasp,exqe,exab
        real(DP),allocatable,public  :: proj_wt(:,:,:)
      contains

      subroutine init_wav
      ! get NKPTS,NBANDS,ISPIN
      ! get the record length for opening WAVECAR

      IMPLICIT REAL(DP) (A-H,O-Z)

! local work arrays
      real(DP)                :: RKPTS,RBANDS,RDUM,RISPIN,RTAG
      integer                 :: IREC,IDUM
      
       inquire (file='WAVECAR', exist=exvasp) !VASP
       inquire (file='data-file-schema.xml', exist=exqe) !QE
       inquire (file=WFCAB,exist=exab)  !abinit
!-----------------------------------------------QE-init
       IF(exqe==.true.) THEN  
       write(IUO,'(/,3X,A)')"Read in parameters from QE_xml"
       call read_spin(lsda,noncolin,spinorbit)
       call read_cutoff(ENMAX)
       call read_cell(CELL)
       call read_nbnd(NBANDS)
       call read_nkp(NKPTS)
          IF (lsda==.true.) THEN
             ISPIN = 2
          ELSE
             ISPIN = 1
          ENDIF
       write(IUO,'(5X,A,I2,A,I4,A,I4)')"ISPIN=",ISPIN, &
     & "  NKPTS=",NKPTS,"  NBANDS=",NBANDS
       write(IUO,'(5X,A)')"Lattice vectors:(.a.u)"
        do i=1,3
          write(IUO,'(3X,3f12.6)')CELL(:,i)
        enddo
!-----------------------------------------------------------abinit-init
       ELSE IF(exab==.true.) THEN
        call read_information_from_abinit(ISPIN,NBANDS,NKPTS,CELL,ENMAX,NRSPINOR)
        write(IUO,'(5X,A,I2,A,I4,A,I4)')"ISPIN=",ISPIN, &
        & "  NKPTS=",NKPTS,"  NBANDS=",NBANDS
        write(IUO,'(5X,A)')"Lattice vectors:(.a.u)"
        do i=1,3
          write(IUO,'(3X,3f12.6)')CELL(:,i)
        enddo
!-----------------------------------------------------------------!VASP-init
       ELSE IF(exvasp==.true.) THEN
       write(IUO,'(/,3X,A)')"Read in parameters from WAVECAR"
       open(IUW,file="WAVECAR",status='old', &
     & recl=16,form='unformatted',access='direct',iostat=iost)

        NRSPINORS =1
        RTAG=0
        READ(IUW,REC=1,ERR=200) RDUM,RISPIN,RTAG
        close(IUW)
        write(IUO,'(5X,A,F20.8)')"RDUM:",RDUM
        write(IUO,'(5X,A,F20.8)')"RISPIN:",RISPIN
        write(IUO,'(5X,A,F20.8)')"RTAG:",RTAG

        close(IUW)

        IF (RTAG==45200) THEN
           SINGLE_PREC=.TRUE.
        ELSE
           WRITE(IUO,'(5X,A)') & 
     & "double precision WAVECAR encountered, converting it"
           SINGLE_PREC=.FALSE.
        ENDIF

        IRECLW=NINT(RDUM)
        write(IUO,'(5X,A,I10)')"IRECLW:",IDUM
        IF ((IRECLW<=0).OR.(IDUM>10000000)) IDUM=16
        ISPIN=NINT(RISPIN)

        OPEN(IUW,FILE='WAVECAR',ACCESS='DIRECT', &
     & FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IRECLW)

        RKPTS=0.0_DP; RBANDS=0.0_DP; ENMAX=0.0_DP
        READ(IUW,REC=2,ERR=200)RKPTS,RBANDS,ENMAX,((CELL(I,J),I=1,3),J=1,3)

        ! check if we get the correct record length 
        if(abs(RKPTS)<1E-6.or.abs(RBANDS)<1E-6.or.abs(ENMAX)<1E-6.or.&
     & abs(RKPTS)>1E6.or.abs(RBANDS)>1E6.or.abs(ENMAX)>5E3)then
          close(IUW)
          write(IUO,'(3X,A)')"Re-open WAVECAR"
          IRECLW=NINT(RDUM)/4
          OPEN(IUW,FILE='WAVECAR',ACCESS='DIRECT', &
     &        FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IRECLW)

          READ(IUW,REC=1,ERR=200) RDUM,RISPIN,RTAG
          write(IUO,'(5X,A,F20.8)')"RDUM:",RDUM
          write(IUO,'(5X,A,F20.8)')"RISPIN:",RISPIN
          write(IUO,'(5X,A,F20.8)')"RTAG:",RTAG

          READ(IUW,REC=2,ERR=200)RKPTS,RBANDS,ENMAX,((CELL(I,J),I=1,3),J=1,3)

          ! check record length again
          if(abs(RKPTS)<1E-6.or.abs(RBANDS)<1E-6.or.abs(ENMAX)<1E-6.or.&
     & abs(RKPTS)>1E6.or.abs(RBANDS)>1E6.or.abs(ENMAX)>5E3)then

            close(IUW)
            write(IUO,'(3X,A)')"Re-open WAVECAR again"
            IRECLW=NINT(RDUM)*4
            OPEN(IUW,FILE='WAVECAR',ACCESS='DIRECT', &
     &        FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IRECLW)

            READ(IUW,REC=1,ERR=200) RDUM,RISPIN,RTAG
            write(IUO,'(5X,A,F20.8)')"RDUM:",RDUM
            write(IUO,'(5X,A,F20.8)')"RISPIN:",RISPIN
            write(IUO,'(5X,A,F20.8)')"RTAG:",RTAG

            READ(IUW,REC=2,ERR=200)RKPTS,RBANDS,ENMAX,((CELL(I,J),I=1,3),J=1,3)

            ! check record length again and again
            if(abs(RKPTS)<1E-6.or.abs(RBANDS)<1E-6.or.abs(ENMAX)<1E-6.or.&
     & abs(RKPTS)>1E6.or.abs(RBANDS)>1E6.or.abs(ENMAX)>5E3)then
              write(IUO,'(3X,A)')"Error: can not find the correct RECL!"
              stop
            endif

          endif
        endif
 
        NKPTS=NINT(RKPTS)
        NBANDS=NINT(RBANDS)
        write(IUO,'(5X,A,I2,A,I4,A,I4)')"ISPIN=",ISPIN, &
     & "  NKPTS=",NKPTS,"  NBANDS=",NBANDS
        write(IUO,'(5X,A)')"Lattice vectors:(A)" 
        do i=1,3
           write(IUO,'(3X,3f12.6)')CELL(:,i)
        enddo

       ELSE
           CALL fstop(IUO,"erro init kproj!!")
       ENDIF!which file 
!---------------------------------------------------------------------
       allocate(CELTOT(NBANDS,NKPTS,ISPIN),FERTOT(NBANDS,NKPTS,ISPIN))
       allocate(NORM_FACTOR(NBANDS,NKPTS,ISPIN))
       allocate(KPTVEC(3,NKPTS))
       kmax = sqrt(2*ENMAX/HARTREE )
       gmax = 2*kmax
       write(IUO,'(5X,A,f12.6)')"kmax:",kmax

       return
  200 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading WAVECAR, header is corrupt'
      STOP
      end subroutine init_wav

      end module mod_wave

