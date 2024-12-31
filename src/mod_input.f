      module mod_input
       use mod_comp

       character(5),public:: task
       integer,public     :: MAT_P2S(3,3)   ! matrix that transforms Primitive Cell -> Supercell
       real(DP)           :: MAT_INVERT(3,3)
       integer,public     :: KSHIFT(3)      ! integer array that shifts the origin of the bands
                                            ! outside the first BZ of the supercell. this is
                                            ! useful if one has the supercell wave functions and
                                            ! wants the projected bands for one of the other folded
                                            ! back k-points. the values are in terms of the
                                            ! reciprocal lattice vectors of the supercell.

       real(DP),public    :: zlay1,zlay2    ! layer boundaries in direct coordinates for layer integration along z
       real(DP),public    :: ylay1,ylay2    ! layer boundaries in direct coordinates for layer integration along y
       real(DP),public    :: xlay1,xlay2    ! layer boundaries in direct coordinates for layer integration along y
       logical,public     :: LZLAYER        ! if do a layer integration along z
       logical,public     :: LYLAYER        ! if do a layer integration along y
       logical,public     :: LXLAYER        ! if do a layer integration along X
       logical,public     :: LSQPSI         ! if output \psi * conjg(\psi) on a FFT mesh
       integer,public     :: SPIN_INDEX,KPT_INDEX,BND_INDEX
       logical,public     :: LPOLA          ! if taking polarization effect into account
       real(DP),public    :: EPSL(3)        ! epsilon: polarization vector
       logical,public     :: LHSE           ! if the WAVECAR from a HSE calculation
       integer,public     :: NK_HSE         ! the number of KPOINTS, of which the weights are non-zero in the HSE calculation
       logical,public     :: LSORBIT        ! if the WAVECAR from SOC/noncollinear calculation
       real(DP),public    :: EINT(2)        ! energy window for bands to be selected
       logical,public     :: LPARD          ! if out PARCHG
       character(len=100),public   :: WFCAB         ! abinit wfcfile_name
       integer,public     :: NTHREAD          ! FFTW parallel
       logical,PUBLIC     :: spin_texture,phonon_unfold,lcao_unfold
!----------------------------------------------------other unfolding
       INTEGER,PUBLIC                   :: orbital_tot,det,ex ,NSPIN
       INTEGER,ALLOCATABLE,PUBLIC       :: proj_part(:)
       REAL(DP),ALLOCATABLE,PUBLIC     :: KPOINT_DIRECT(:,:),KPOINT_CARTESIAN(:,:)
       REAL(DP),PUBLIC                 :: CELL(3,3),FACTOR,CELL_P(3,3)
       INTEGER ,ALLOCATABLE,PUBLIC     :: POS_MAP_DIRECT(:,:,:,:)
       REAL(DP),ALLOCATABLE,PUBLIC     :: POS_MAP_CARTESIAN(:,:,:,:)
       INTEGER,ALLOCATABLE,PUBLIC      :: ions_type(:),tot_ions(:),ions_n(:,:),orbital_n(:,:)
       character(len=20),allocatable,public       :: ions_name(:,:)
       INTEGER,PUBLIC                  :: part

       contains
     
       subroutine input
         logical      :: LDUM, lopen
         character    :: CHARAC
         integer      :: IDUM, IERR, N, II
         real(DP)     :: RDUM,fdum
         complex(DPC) :: CDUM
         CHARACTER(len=100)    :: SZNAM
         !---------------------------------------------

         inquire (file='INKPROJ', exist=ex)
         if (ex.eq..false.)  & 
     & call fstop(IUO, 'Error: file = INKPROJ is not found, ABORTING!')

         open(IUI,file='INKPROJ',status='old')
         write(IUO,'(/,3x,A)')"Read in control paramters"

         lopen = .FALSE.
          call RDATAB(LOPEN,'INKPROJ',IUI,'WFCAB','=','#','//;','S', &
     &            IDUM,RDUM,CDUM,LDUM,WFCAB,N,100,IERR)
          IF ((IERR/=0).AND.(IERR/=3)) THEN
          IF (IUO>=0) &
          WRITE(IUO,*)'Error reading item ''wfc'' from file INCAR.'
          GOTO 150
          ENDIF

         
        inquire (file='INKPROJ', exist=ex)
        if (ex.eq..false.) then
        print*,"Error: file = INKPROJ is not found, ABORTING!"
        stop
        endif

        spin_texture=.false.
        call RDATAB(LOPEN,'INKPROJ',IUI,'spin_texture','=','#',';','L', &
    &            IDUM,RDUM,CDUM,spin_texture,CHARAC,N,1,IERR)
        if (((IERR/=0).AND.(IERR/=3)).OR. &
    &                    ((IERR==0).AND.(N<1))) THEN
           write(IUO,*) &
    & 'Error reading item ''spin_texture'' from file INKPROJ.'
        endif
           
        phonon_unfold=.false.
!        call RDATAB(LOPEN,'INKPROJ',IUI,'phonon_unfold','=','#',';','L', &
!    &            IDUM,RDUM,CDUM,phonon_unfold,CHARAC,N,1,IERR)
!        if (((IERR/=0).AND.(IERR/=3)).OR. &
!    &                    ((IERR==0).AND.(N<1))) THEN
!           write(IUO,*) &
!    & 'Error reading item ''phonon_unfold'' from file INKPROJ.'
!        endif

        lcao_unfold=.false.
!        call RDATAB(LOPEN,'INKPROJ',IUI,'lcao_unfold','=','#',';','L', &
!    &            IDUM,RDUM,CDUM,lcao_unfold,CHARAC,N,1,IERR)
!        if (((IERR/=0).AND.(IERR/=3)).OR. &
!    &                    ((IERR==0).AND.(N<1))) THEN
!           write(IUO,*) &
!    & 'Error reading item ''lcao_unfold'' from file INKPROJ.'
!        endif

       IF (lcao_unfold .and. phonon_unfold ==.true.) then
          print*,"erro!"
          stop
       ENDIF

         NSPIN = 1
         call rdatab(lopen,'INKPROJ',IUI,'NSPIN','=','#',';', 'I', &
     &     NSPIN, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
                     &  ((IERR==0).AND.(N<1))) THEN
           write(IUO,*) &
    & 'Error reading item ''NSPIN'' from file INKPROJ.'
        endif

       
       IF (NSPIN/=1 .and. NSPIN/=2) THEN
          PRINT*,"erro! NSPIN must equal 1 or 2"
          stop
       ENDIF        


         det = 1
         call rdatab(lopen,'INKPROJ',IUI,'det','=','#',';', 'I', &
     &     det, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
           write(IUO,'(5X,A)') &
     & 'error in reading item "det" from file INKPROJ.'
         endif


         ions_type = 1
         call rdatab(lopen,'INKPROJ',IUI,'ions_type','=','#',';', 'I', &
     &     ions_type, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
           write(IUO,'(5X,A)') &
     & 'error in reading item "ions_type" from file INKPROJ.'
         endif


!        ions_n=0
!        call rdatab(lopen,'INKPROJ',IUI,'ions_n', '=', '#', ';', 'I', &
!     &     ions_n(1), RDUM, CDUM, LDUM, CHARAC, N, ions_type, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!     &      ((IERR==0).and.(N<1))) then
!            write(IUO,'(5X,A)') &
!     & 'error in reading item "ions_n" from file INKPROJ.'
!         endif


!         allocate(proj_part(part))
!         proj_part=0
!        call rdatab(lopen,'INKPROJ',IUI,'proj_part', '=', '#', ';', 'I', &
!     &     proj_part(1), RDUM, CDUM, LDUM, CHARAC, N, part, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!     &      ((IERR==0).and.(N<1))) then
!            write(IUO,'(5X,A)') &
!     & 'error in reading item "proj_part" from file INKPROJ.'
!         endif
!         print*,proj_part
!
!        start_type=1
!        call rdatab(lopen,'INKPROJ',IUI,'start_type','=','#',';', 'I', &
!     &     start_type, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!                     &  ((IERR==0).AND.(N<1))) THEN
!           write(IUO,*) &
!    & 'Error reading item ''select_type'' from file INKPROJ.'
!        endif
!        end_type=ions_type
!        call rdatab(lopen,'INKPROJ',IUI,'end_type','=','#',';', 'I', &
!              &  end_type, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!                     &  ((IERR==0).AND.(N<1))) THEN
!           write(IUO,*) &
!    & 'Error reading item ''select_type'' from file INKPROJ.'
!        endif




!        ALLOCATE(orbital_n(ions_type))
!        
!        orbital_n=0
!        call rdatab(lopen,'INKPROJ',IUI,'orbital_n', '=', '#', ';', 'I', &
!     &     orbital_n(1), RDUM, CDUM, LDUM, CHARAC, N, ions_type, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!     &      ((IERR==0).and.(N<1))) then
!            write(IUO,'(5X,A)') &
!     & 'error in reading item "orbital_n" from file INKPROJ.'
!         endif



!         NKPOINT = 1
!         call rdatab(lopen,'INKPROJ',IUI,'NKPOINT','=','#',';', 'I', &
!     &     NKPOINT, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
!         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
!     &      ((IERR==0).and.(N<1))) then
!           write(IUO,'(5X,A)') &
!     & 'error in reading item "NKPOINT" from file INKPROJ.'
!         endif


    

!   IF (lcao_unfold .or. spin_texture ==.true.) then
!
! !  ALLOCATE(POS_MAP_DIRECT(ions_type,maxval(ions_n),3),POS_MAP_CARTESIAN(ions_type,maxval(ions_n),3))
!   inquire (file='pos_map.dat', exist=ex)
!   IF (ex == .false.) then
!      print*,"Error: file = pos_map.dat is not found, ABORTING!"
!      stop
!   ENDIF
!ENDIF
!     OPEN(unit=9,file="pos_map.dat")
!
!
!     CLOSE(9)
!      
     !ALLOCATE(KPOINT_DIRECT(NKPOINT,3),KPOINT_CARTESIAN(NKPOINT,3))
     !KPOINT_DIRECT=0


     !----------------------------------------------------------------------------
         LZLAYER = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LZLAYER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LZLAYER,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==LSORBIT).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LZLAYER'' from file INKPROJ.'
            GOTO 150
         endif
         if(LZLAYER)then
           write(IUO,'(5x,(A,L,2X,A))')"LZLAYER:",LZLAYER, &
     & "(Layer integration will be done in [zlay1, zlay2]"
         else
           write(IUO,'(5x,(A,L,2X,A))')"LZLAYER:",LZLAYER, &
     & "(Layer intergration will be done for the whole space along Z)"
         endif

         zlay1 = 0.0
         call RDATAB( lopen, 'INKPROJ', IUI, 'zlay1', '=', '#',';','F', &
     &        IDUM, zlay1, CDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''zlay1'' from file INKPROJ.'
            goto 150
         endif

         zlay2 = 0.5
         call RDATAB(lopen,'INKPROJ', IUI, 'zlay2', '=', '#', ';', 'F', &
     &        IDUM, zlay2, MAT_P2SCDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''zlay2'' from file INKPROJ.'
            goto 150
         endif
    
         if (zlay2 < zlay1) then
            write(IUO,&
     & '(/,''***************************************'',/)')
            write(IUO,*)' WARNING: zlay2 < zlay1'
            write(IUO,*)' setting zlay2 and zlay1 correctly'
            write(IUO,&
     & '(/,''***************************************'',/)')

            fdum=zlay1
            zlay1=zlay2
            zlay2=fdum   
            write(IUO,*)'We LSORBITexchange zlay1 and zlay2'
         end if
         if(LZLAYER)then
           write(IUO,'(5x,A,f8.2)')"zlay1:",zlay1
           write(IUO,'(5x,A,f8.2)')"zlay2:",zlay2
         endif

         LYLAYER = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LYLAYER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LYLAYER,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LYLAYER'' from file INKPROJ.'
            GOTO 150
         endif
         if(LYLAYER)then
           write(IUO,'(5x,(A,L,2X,A))')"LYLAYER:",LYLAYER, &
     & "(Layer integration will be done in [ylay1, ylay2]"
         else
           write(IUO,'(5x,(A,L,2X,A))')"LYLAYER:",LYLAYER, &
     & "(Layer intergration will be done for the whole space along y)"
         endif

         ylay1 = 0.0
         call RDATAB( lopen, 'INKPROJ', IUI, 'ylay1', '=', '#',';','F', &
     &        IDUM, ylay1, CDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''ylay1'' from file INKPROJ.'
            goto 150
         endif

         ylay2 = 0.5
         call RDATAB(lopen,'INKPROJ', IUI, 'ylay2', '=', '#', ';', 'F', &
     &        IDUM, ylay2, CDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''ylay2'' from file INKPROJ.'
            goto 150
         endif
    
         if (ylay2 < LSORBITylay1) then
            write(IUO,&
     & '(/,''***************************************'',/)')
            write(IUO,*)' WARNING: ylay2 < ylay1'
            write(IUO,*)' setting ylay2 and ylay1 correctly'
            write(IUO,&
     & '(/,''***************************************'',/)')

            fdum=ylay1
            ylay1=ylay2
            ylay2=fdum   
            write(IUO,*)'We exchange ylay1 and ylay2'
         end if
         if(LYLAYER)then
           write(IUO,'(5x,A,f8.2)')"ylay1:",ylay1
           write(IUO,'(5x,A,f8.2)')"ylay2:",ylay2
         endif

         LXLAYER = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LXLAYER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LXLAYER,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LXMAT_P2SLAYER'' from file INKPROJ.'
            GOTO 150
         endif
         if(LXLAYER)then
           write(IUO,'(5x,(A,L,2X,A))')"LXLAYER:",LXLAYER, &
     & "(Layer integration will be done in [xlay1, xlay2]"
         else
           write(IUO,'(5x,(A,L,2X,A))')"LXLAYER:",LXLAYER, &
     & "(Layer intergration will be done for the whole space along X)"
         endif

         xlay1 = 0.0
         call RDATAB( lopen, 'INKPROJ', IUI, 'xlay1', '=', '#',';','F', &
     &        IDUM, xlay1, CDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''xlay1'' from file INKPROJ.'
            goto 150
         endif

         xlay2 = 0.5
         call RDATAB(lopen,'INKPROJ', IUI, 'xlay2', '=', '#', ';', 'F', &
     &        IDUM, xlay2, CDUM, LDUM, CHARAC, N, 1, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''xlay2'' from file INKPROJ.'
            goto 150
         endif
    
         if (xlay2 < xlay1) then
            write(IUO,&
     & '(/,''***************************************'',/)')
            write(IUO,*)' WARNING: xlay2 < xlay1'
            write(IUO,*)' setting xlay2 and xlay1 correctly'
            write(IUO,&
     & '(/,''***************************************'',/)')

            fdum=xlay1
            xlay1=xlay2
            xlay2=fdum   
            write(IUO,*)'We exchange xlay1 and xlay2'
         end if
         if(LXLAYER)then
           write(IUO,'(5x,A,f8.2)')"xlay1:",xlay1
           write(IUO,'(5x,A,f8.2)')"xlay2:",xlay2
         endif

         NTHREAD = 2
         call rdatab(lopen,'INKPROJ',IUI,'NTHREAD', '=', '#', ';', 'I', &
     &     NTHREAD, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "NTHREAD" from file INKPROJ.'
            goto 150
         endif


         MAT_P2S(:,:) = 0
         do II = 1,3
            MAT_P2S(II,II) = 1        ! original matrix is diagonal
         enddo
         call rdatab(lopen,'INKPROJ',IUI,'MAT_P2S', '=', '#', ';', 'I', &
     &     MAT_P2S(1,1), RDUM, CDUM, LDUM, CHARAC, N, 9, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "MAT_P2S" from file INKPROJ.'
            goto 150
         endif

         write(IUO,'(5X,A)') &
     & "Transform matrix between primitive and super cells:"
         do II = 1,3
            write(IUO,'(5X,3I4)')MAT_P2S(1:3,II)
         enddo
 
        inquire (file="running_nscf.log",exist=ex)
        IF (ex == .true.) then
           lcao_unfold = .true.
           call  read_posmap
        ENDIF
        inquire (file="band.yaml",exist=ex)
        IF (ex == .true.) then
           phonon_unfold = .true.
           call  read_posmap
        ENDIF

         allocate(proj_part(part))
         proj_part=0
        call rdatab(lopen,'INKPROJ',IUI,'proj_part', '=', '#', ';', 'I', &
     &     proj_part(1), RDUM, CDUM, LDUM, CHARAC, N, part, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "proj_part" from file INKPROJ.'
         endif



         KSHIFT(:) = 0
         call rdatab(lopen,'INKPROJ',IUI,'KSHIFT', '=', '#', ';', 'I', &
     &     KSHIFT(1), RDUM, CDUM, LDUM, CHARAC, N, 3, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "KSHIFT" from file INKPROJ.'
            goto 150
         endif

         if( any( KSHIFT /= 0 ) ) then
            write(IUO,'(5x,(A,3(I3,A)))') &
     & "k-points shifted by ",KSHIFT(1),"b_1 +",KSHIFT(2),"b_2 +", &
     & KSHIFT(3),"b_3"
            write(IUO,'(7X,A)')"b_i are reciprocal lattice vectors"
         endif

         ! if there is a polarization vector
         LPOLA = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LPOLA','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LPOLA,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LPOLA'' from file INKPROJ.'
            GOTO 150
         endif

         if(LPOLA)then
           write(IUO,'(5x,(A,L,2X,A))')"LPOLA:",LPOLA, &
     & "(Taking polarization effect into account)"

         EPSL(1) = 1.0; EPSL(2) = 0.0; EPSL(3) = 0.0
         call RDATAB(lopen,'INKPROJ',IUI,'EPSL', '=', '#', ';', 'F', &
     &        IDUM, EPSL(1), CDUM, LDUM, CHARAC, N, 3, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''EPSL'' from file INKPROJ.'
            goto 150
         endif

         write(IUO,'(5X,A,3f6.2)')"Polarization vector:",EPSL(1:3)
         endif    ! LPOLA

         LSQPSI = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LSQPSI','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LSQPSI,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LSQPSI'' from file INKPROJ.'
            GOTO 150
         endif
         if(LSQPSI)then
           write(IUO,'(5x,(A,L,2X,A))')"LSQPSI:",LSQPSI, &
     & "Output square of wavefunctions on a FFT mesh"
           write(IUO,'(7X,A)')&
     & "Make sure that LPOLA must be .FALSE."
         endif

         if(LSQPSI)then
         SPIN_INDEX = 1
         call rdatab(lopen,'INKPROJ',IUI,'SPIN_INDEX','=','#',';', 'I', &
     &     SPIN_INDEX, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
           write(IUO,'(5X,A)') &
     & 'error in reading item "SPIN_INDEX" from file INKPROJ.'
            goto 150
         endif

         KPT_INDEX = 1
         call rdatab(lopen,'INKPROJ',IUI,'KPT_INDEX','=','#',';', 'I', &
     &     KPT_INDEX, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "KPT_INDEX" from file INKPROJ.'
            goto 150
         endif

         BND_INDEX = 1
         call rdatab(lopen,'INKPROJ',IUI,'BND_INDEX','=','#',';', 'I', &
     &     BND_INDEX, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "BND_INDEX" from file INKPROJ.'
            goto 150
         endif

         write(IUO,'(7X,A,1X,3(A,1X,I0,2X))')"Indices for |\psi|^2:",&
     & "spin",SPIN_INDEX,"kpoint",KPT_INDEX,"band",BND_INDEX
         endif

         LHSE = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LHSE','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LHSE,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LHSE'' from file INKPROJ.'
            GOTO 150
         endif
         if(LHSE)then
           write(IUO,'(5x,(A,L,2X,A))')"LHSE:",LHSE, &
     & "WAVECAR was obtained from a HSE calculation."
           write(IUO,'(7X,A)')&
     & "Please set NK_HSE correspondingly, & 
     & which is the number of k-points with non-zero weights"
         endif

         if(LHSE)then
         NK_HSE = 0
         call rdatab(lopen,'INKPROJ',IUI,'NK_HSE','=','#',';', 'I', &
     &     NK_HSE, RDUM, CDUM, LDUM, CHARAC, N, 1, IERR)
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,'(5X,A)') &
     & 'error in reading item "NK_HSE" from file INKPROJ.'
            goto 150
         endif

         if(NK_HSE<0)then
           NK_HSE = 0
           write(IUO,'(5x,A)')"NK_HSE can not be negtive, set it to 0"
         endif

         write(IUO,'(5x,A,I4)')"NK_HSE:",NK_HSE
         endif

         ! if there is a polarization vector
         LSORBIT = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LSORBIT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LSORBIT,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LSORBIT'' from file INKPROJ.'
            GOTO 150
         endif
          write(IUO,'(5x,(A,L,2X,A))')"LSORBIT:",LSORBIT, &
     & "If WAVECAR was obtained from a SOC/noncollinear calculation."

         ! if OUTPUT PARCHG
         LPARD = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LPARD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LPARD,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
            write(IUO,*) &
     & 'Error reading item ''LPARD'' from file INKPROJ.'
            GOTO 150
         endif

         if(LPARD)then
         EINT(1) = 0.0; EINT(2) = 0.0
         call RDATAB(lopen,'INKPROJ',IUI,'EINT', '=', '#', ';', 'F', &
     &        IDUM, EINT(1), CDUM, LDUM, CHARAC, N, 2, IERR )
         if ( ( (IERR/=0) .and. (IERR/=3) ).or. &
     &      ((IERR==0).and.(N<1))) then
            write(IUO,*)'error reading item ''EINT'' from file INKPROJ.'
            goto 150
         endif

         if(EINT(1)>EINT(2))then
           fdum = EINT(1)
           EINT(1) = EINT(2)
           EINT(2) = fdum 
         endif
          
         write(IUO,'(5X,A)')"LPARD = .TRUE.. OUTPUT PARCHG."
         write(IUO,'(5X,A,3f6.2)')"Energy window for PARCHG:",EINT(1:2)
         endif
         
         close(IUI)

         task = "kproj"         ! default: do kprojction
         if(LSQPSI)then
           task = "psi2"        ! output wavefunction
         elseif(LPARD)then
           task = "pchg"        ! output PARCHG
         endif

         RETURN
150      continue
         write(IUO,151) IERR, N

151      format(' Error code was IERR=',I1,' ... . Found N=',I5,' data.')
         STOP
      end subroutine input

        subroutine read_posmap
        integer::type_index,det_matp2s
        det_matp2s = MAT_P2S(1,1) * (MAT_P2S(2,2) * MAT_P2S(3,3) - MAT_P2S(3,2) * MAT_P2S(2,3)) - &
                     MAT_P2S(2,1) * (MAT_P2S(1,2) * MAT_P2S(3,3) - MAT_P2S(3,2) * MAT_P2S(1,3)) + &
                     MAT_P2S(3,1) * (MAT_P2S(1,2) * MAT_P2S(2,3) - MAT_P2S(2,2) * MAT_P2S(1,3))         
                     
        inquire (file='pos_map.dat', exist=ex)
        IF (ex == .false.) then
           print*,"Error: file = pos_map.dat is not found, ABORTING!"
           stop
        ENDIF

        open(unit=99,file="pos_map.dat")
        read(99,*)part
        allocate(ions_type(part),tot_ions(part))
        read(99,*)ions_type(:)
        read(99,*)tot_ions(:)
        allocate(orbital_n(part,maxval(ions_type)))
        allocate(ions_n(part,maxval(ions_type)))
        allocate(ions_name(part,maxval(ions_type)))
        allocate(POS_MAP_DIRECT(part,maxval(ions_type),maxval(tot_ions),3))
        allocate(POS_MAP_CARTESIAN(part,maxval(ions_type),maxval(tot_ions),3))
        ions_n=0
        POS_MAP_DIRECT=0
        POS_MAP_CARTESIAN=0
        type_index=1
        do i = 1, part
          do j =1,tot_ions(i)
             read(99,*)type_index,POS_MAP_DIRECT(i,type_index,ions_n(i,type_index)+1,:),ions_name(i,type_index)  ! part,ions_type,ions_n,shift
             ions_index=ions_n(i,type_index)+1
             ions_n(i,type_index)=ions_n(i,type_index)+1
          enddo
        enddo
     close(99)
     end subroutine


      end module mod_input
