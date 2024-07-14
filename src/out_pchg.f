      subroutine out_pchg

      use mod_comp
      use mod_input
      use mod_wave
      use mod_fft
      use mod_lattice
      use read_information_from_qe
!              fft mesh arrays
      complex(DP),allocatable :: q_pw(:,:,:,:)

!              WAVECAR related
      complex(DP),allocatable :: WFC(:,:)
      integer,allocatable     :: GVECT(:,:)
      complex(DP),allocatable :: CW2(:)
      complex(QS),allocatable :: CRD(:)
      real(DP)                :: VKPT(3)
      integer                 :: IREC,IDUM,spin_direction
      real(DP)                :: RNPL,NORM_TEST,IMEN
      real(DP)                :: b1_l(3),b2_l(3),b3_l(3),scalef
      integer                 :: GVECTOR(3)
      integer                 :: i,j,k,l,n,ib,isp,ik,ip,ic,ios,ikpt
      integer                 :: NSPINOR,NPLG
      CHARACTER(len=100)      :: wfcname,ikp,spin
      real(DP),allocatable    :: pchg(:,:,:),weight(:),eigenvalues(:,:)
      logical                 :: gamma_only
      write(IUO,'(/,3X,A)')"Write out PARCHG, &
     & which exactly has the same unit as that from VASP ( &
     &\rho * Volume )"

      write(IUO,'(/,47X,A,f8.4,A,f8.4,2X,A)')"[",EINT(1),",",EINT(2),"]"


      call fft_allocate( q_pw, nfft, NSPINOR ) 
      allocate(pchg(nfft(1),nfft(2),nfft(3)))
      pchg = 0.0_DP
!-------------------------------------------------------------    
      IF (exqe==.true.) THEN
        IF (noncolin == .true.) THEN
          NSPINOR = 2
        ELSE
          NSPINOR = 1
      ENDIF
      allocate(weight(NKPTS))
      call read_kpoint_energy(KPTVEC,eigenvalues,weight)
      do ib=1,NBANDS
        do isp = 1,ISPIN
        FERTOT(ib,:,isp)=weight(:)
        enddo
      enddo
      do isp=1,ISPIN !spin
        do ik = 1,NKPTS !kpoint
        !
        IF (ISPIN==1) spin =""
        IF (isp==1 .and. ISPIN==2 ) spin ="up"
        IF (isp==2 .and. ISPIN==2 ) spin ="dw"
        write(ikp,*) ik
        wfcname = "wfc"//trim(spin)//trim(adjustl(ikp))//".dat"
        OPEN ( UNIT = IUO, FILE = wfcname , &
              FORM='unformatted', STATUS = 'unknown' )
        !
        read(IUO) ikpt,(VKPT(i),i=1,3),spin_direction,gamma_only,scalef
        read(IUO) ngw,NPLG,NSPINOR,NBANDS   !当有自旋时候，这里NBANDS就是对应自旋的能带数
        ! 
        IF(spin_direction == 1) THEN
        CELTOT(:,ik,isp)=eigenvalues(1:NBANDS,ik)
        ELSE IF(spin_direction == 2 ) THEN
        CELTOT(:,ik,isp)=eigenvalues(NBANDS+1:2*NBANDS,ik)
        END IF
        read(IUO) b1_l,b2_l,b3_l
        !
        allocate (CW2(NPLG*NSPINOR))
        allocate (GVECT(3,NPLG))
        allocate (WFC(NPLG,NSPINOR))
        read(IUO) GVECT
        call set_nsphk(NPLG,GVECT)
       
        do ib = 1,NBANDS !band
           if(.not.(real(CELTOT(ib,ik,isp),DP)>EINT(1).and. &
           & real(CELTOT(ib,ik,isp),DP)<EINT(2)))cycle

           write(IUO,'(5X,A,I4,2X,A,3f8.4,1X,A,2X,A,I4,2X,A,f10.5)') &
           & "kpoint",ik,"(",KPTVEC(1:3,ik),")","band",ib,"E=", &
           & real(CELTOT(ib,ik,isp),DP)

          read(IUO) CW2    !CW2
          !  normalization     
          NORM_TEST = 0.0_DP
          DO j=1,NSPINOR
            DO i=1,NPLG
              WFC(i,j) = CW2(i+(j-1)*NPLG)
              NORM_TEST = NORM_TEST + WFC(i,j)*conjg(WFC(i,j))
            ENDDO
          ENDDO
            DO j=1,NSPINOR
               call wf_to_fftmesh(q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
               call fft( q_pw(:,:,:,j), 'k->r', nfft )
            ENDDO
               ! u* x u
               q_pw = q_pw*conjg(q_pw)
               if(NSPINOR>1)q_pw(:,:,:,1) = q_pw(:,:,:,1) +q_pw(:,:,:,2)

               do i=1,nfft(1)
                 do j=1,nfft(2)
                   do k=1,nfft(3)
                      pchg(i,j,k) = pchg(i,j,k) + &
        & real(q_pw(i,j,k,1)/NORM_TEST,DP)
                   enddo
                 enddo
               enddo

        enddo! band
        deallocate(CW2)
        deallocate(GVECT)
        deallocate(WFC)
        enddo! kpoint
      enddo! spin      
!-------------------------------------------------------------
      ELSE IF (exvasp==.true.) THEN
!-------------------------------------------------------------    
      NSPINOR = 1
      IF(LSORBIT) NSPINOR = 2
      IREC = 2
        do isp=1,ISPIN

          do ik=1,NKPTS 

        ! read WAVECAR, get the number of plane waves, k vector
        ! eigenvalue
        !
        IREC=IREC+1

        READ(IUW,REC=IREC,ERR=220) RNPL,VKPT, &
     & (CELTOT(J,ik,isp),IMEN,FERTOT(J,ik,isp),J=1,NBANDS)
        KPTVEC(1:3,ik)=VKPT(1:3)
        NPL=NINT(RNPL)

        ! the number of plane waves has been doubled for SOC
        NPLG = NPL/NSPINOR

        ! allocate array for G vectors and wave functions
        allocate(GVECT(3,NPLG),WFC(NPL,NSPINOR))

        ! generate G vectors
        call gen_gvect(VKPT,NPLG,NPL_TMP,GVECT)
        !
        ! check if we have generated correct G vectors
        !
        if(NPL_TMP/=NPLG) then
!         write(IUO,'(5X,A,2I6)') &
!    & "NPL_TMP and read NPLG:",NPL_TMP,NPLG
          call fstop(IUO,'Error in generating G vectors ')
        endif

        call set_nsphk(NPLG,GVECT)

        allocate(CW2(NPL),CRD(NPL))

        band:  do ib = 1,NBANDS
            
              ! get coefficients of wave functions from WAVECAR
              IREC=IREC+1

              if(.not.(real(CELTOT(ib,ik,isp),DP)>EINT(1).and. &
     & real(CELTOT(ib,ik,isp),DP)<EINT(2)))cycle

              write(IUO,'(5X,A,I4,2X,A,3f8.4,1X,A,2X,A,I4,2X,A,f10.5)') &
     & "kpoint",ik,"(",KPTVEC(1:3,ik),")","band",ib,"E=", &
     & real(CELTOT(ib,ik,isp),DP)

              IF (SINGLE_PREC) THEN
                READ(IUW,REC=IREC,ERR=240) (CRD(i),i=1,NPL)
                CW2(1:NPL)=CRD(1:NPL)
              ELSE
                READ(IUW,REC=IREC,ERR=240) (CW2(i),i=1,NPL)
              ENDIF

              NORM_TEST = 0.0_DP
              DO j=1,NSPINOR
              DO i=1,NPLG
                 WFC(i,j) = CW2(i+(j-1)*NPLG)
                 NORM_TEST = NORM_TEST + WFC(i,j)*conjg(WFC(i,j))
              ENDDO
              ENDDO

               ! using FFT to get the Bloch function u from G space to
               ! real space
               !
               DO j=1,NSPINOR
               call wf_to_fftmesh(q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
               call fft( q_pw(:,:,:,j), 'k->r', nfft )
               ENDDO
               ! u* x u
               q_pw = q_pw*conjg(q_pw)
               if(NSPINOR>1)q_pw(:,:,:,1) = q_pw(:,:,:,1) +q_pw(:,:,:,2)

               do i=1,nfft(1)
                 do j=1,nfft(2)
                   do k=1,nfft(3)
                      pchg(i,j,k) = pchg(i,j,k) + &
     & real(q_pw(i,j,k,1)/NORM_TEST,DP)
                   enddo
                 enddo
               enddo
              
        enddo band
        deallocate(GVECT,WFC)
        deallocate(CW2,CRD) 
      enddo! kpoint
      enddo! spin 
!---------------------------------------------------------
       ELSE IF(exab==.true.) THEN
       IF (NRSPINORS == 2) THEN
          NSPINOR = 2
        ELSE
          NSPINOR = 1
        ENDIF

       KPTVEC=kpt_ab

       do isp=1,ISPIN
        do ik=1,NKPTS
          READ(IUW) NPLG,NSPINOR,NBANDS

          ALLOCATE(GVECT(3,NPLG))
          ALLOCATE(CW2(NPLG*NSPINOR))
          ALLOCATE(WFC(NPLG,NSPINOR))

          READ(IUW) GVECT
          ! set weights for k-projected bands
          call set_nsphk(NPLG,GVECT)


          READ(IUW) CELTOT(:,ik,isp),FERTOT(:,ik,isp)
           CELTOT(:,ik,isp)=CELTOT(:,ik,isp)*27.211652_DP
            ! get coefficients of wave functions from _WFK
           do ib=1,NBANDS
            !
             READ(IUW) CW2
            !  normalization
             NORM_TEST = 0.0_DP
             DO j=1,NSPINOR
              DO i=1,NPLG
                WFC(i,j) = CW2(i+(j-1)*NPLG)
                NORM_TEST = NORM_TEST + WFC(i,j)*conjg(WFC(i,j))
              ENDDO
             ENDDO
             NORM_FACTOR(ib,ik,isp) = NORM_TEST

           ! using FFT to get the Bloch function u from G space to
               ! real space
               !
               DO j=1,NSPINOR
               call wf_to_fftmesh(q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
               call fft( q_pw(:,:,:,j), 'k->r', nfft )
               ENDDO
               ! u* x u
               q_pw = q_pw*conjg(q_pw)
               if(NSPINOR>1)q_pw(:,:,:,1) = q_pw(:,:,:,1) +q_pw(:,:,:,2)

               do i=1,nfft(1)
                 do j=1,nfft(2)
                   do k=1,nfft(3)
                      pchg(i,j,k) = pchg(i,j,k) + &
                      & real(q_pw(i,j,k,1)/NORM_TEST,DP)
                   enddo
                 enddo
               enddo

        enddo ! band
        deallocate(GVECT,WFC)
        deallocate(CW2)
      enddo! kpoint
      enddo! spin
      end if
!---------------------------------------------------------
      open(100,file="PARCHG_KPROJ",status='unknown')
      write(100,'(A)')"PARCHG from program kproj "
      write(100,'(f18.10)')1.0_DP
      write(100,'(3f18.10)')CELL(1:3,1)
      write(100,'(3f18.10)')CELL(1:3,2)
      write(100,'(3f18.10)')CELL(1:3,3)

      write(100,'(/,3I5)')nfft(1:3)

      write(100,'(5(1X,E18.11))') &
     & (((pchg(i,j,k), &
     & i=1,nfft(1)),j=1,nfft(2)),k=1,nfft(3))

      close(100)

      RETURN

  220 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
      STOP

  240 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP

      end subroutine out_pchg 
