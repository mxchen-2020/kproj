      subroutine out_psi2

      use mod_comp
      use mod_input
      use mod_wave
      use mod_fft
      use mod_lattice

!              fft mesh arrays
      complex(DP),allocatable :: q_pw(:,:,:,:)

!              WAVECAR related
      complex(DP),allocatable :: WFC(:,:)
      integer,allocatable     :: GVECT(:,:)
      complex(DP),allocatable :: CW2(:)
      complex(QS),allocatable :: CRD(:)
      real(DP)                :: VKPT(3)
      integer                 :: IREC,IDUM
      real(DP)                :: RNPL,NORM_TEST,IMEN

      integer                 :: GVECTOR(3)
      integer                 :: i,j,k,l,n,ib,isp,ik,ip,ic,ios
      integer                 :: NSPINOR,NPLG

      character(50)    :: filename

      write(IUO,'(/,3X,A)')"Write out |\psi|^2, &
     & which is in the form of \rho * Volume/(Nx,Ny,Nz)"

      NSPINOR = 1
      IF(LSORBIT) NSPINOR = 2

      call fft_allocate( q_pw, nfft, NSPINOR ) 
 
      IREC = 2
      spin: do isp=1,ISPIN

      kpoint:  do ik=1,NKPTS 

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
!              write(IUO,'(5X,A,2I6)') &
!    & "NPL_TMP and read NPLG:",NPL_TMP,NPLG
             call fstop(IUO,'Error in generating G vectors ')
          endif

          call set_nsphk(NPLG,GVECT)

          allocate(CW2(NPL),CRD(NPL))

          band: do ib = 1,NBANDS
            
            ! get coefficients of wave functions from WAVECAR
            IREC=IREC+1

            if(isp==SPIN_INDEX.and.ik==KPT_INDEX.and.ib==BND_INDEX)then 

              write(IUO,'(5X,A,I2,2X,A,I3,2X,A,I5)') &
     & "Spin:",isp,"kpoint:",ik,"band:",ib

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
!             NORM_TEST = 0.0_DP
!             DO i=1,NPL
!                NORM_TEST = NORM_TEST + CW2(i)*conjg(CW2(i))
!             ENDDO

!             write(IUO,'(/,3X,A,f12.6,/)')"Normalized factor:",NORM_TEST
 
              ! using FFT to get the Bloch function u from G space to real space 
              !
              DO j=1,NSPINOR
              call wf_to_fftmesh( q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
              call fft( q_pw(:,:,:,j), 'k->r', nfft )
              ENDDO
                
              ! u* x u
              q_pw = q_pw*conjg(q_pw)
              if(NSPINOR>1)then
                q_pw(:,:,:,1) = q_pw(:,:,:,1) + q_pw(:,:,:,2)
              endif
                
              open(100,file="PSI2.dat",status='unknown')
              write(100,'(A,I0,2X,A,I0,3f12.6,2X,A,I0)') &
     & "SPIN:",isp,"K-POINT:",ik,VKPT(1:3),"BAND:",ib

              write(100,'(3f18.10)')CELL(1:3,1)
              write(100,'(3f18.10)')CELL(1:3,2)
              write(100,'(3f18.10)')CELL(1:3,3)

              write(100,'(/,3I5)')nfft(1:3)

              write(100,'(5(1X,E18.11))') &
     & (((real(q_pw(i,j,k,1)/NORM_TEST/(nfft(1)*nfft(2)*nfft(3)),DP), &
     & i=1,nfft(1)),j=1,nfft(2)),k=1,nfft(3))

              close(100)
            endif

          enddo band
          deallocate(GVECT,WFC)
          deallocate(CW2,CRD) 
      enddo kpoint
      enddo spin

      return

  220 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
      STOP

  240 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP

      end subroutine out_psi2 
