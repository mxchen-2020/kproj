      module mod_kproj
        use mod_comp
        use mod_kproj_qe
        use mod_step_funciton
        use mod_kproj_abinit
      ! real(DP),allocatable,public  :: proj_wt(:,:,:)
      integer*8,public                :: plan,iret  !FFTW
      contains

      subroutine bands_wt_all

!     calculates the weight of states projected on to the unit cell over all space.
!     03/2013 (C.M.X.)

      use mod_input
      use mod_wave
      use mod_fft
      use mod_lattice

!              k-projected information
      real(DP),allocatable    :: wt_proj(:,:),wt_tmp(:)
      real(DP)                :: inv_mat(3,3),RMAT_P2S(3,3)
      real(DP)                :: eta(3)

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
      real(DP)                :: EPSL2,RT_EPSL2,KPG(3),KPG_CART(3)
      real(DP)                :: RE_CELL(3,3),KSHIFT_TMP(3)

      write(IUO,'(/,3X,A)')"Projection over the whole space"
!--------------------------------------------------------------------    
      IF(exqe==.true.) THEN
      call bands_wt_all_qe
      ELSE IF(exab==.true.) THEN
      call bands_wt_all_abinit
      else
!--------------------------------------------------------------------
      NSPINOR = 1
      IF(LSORBIT) NSPINOR = 2
      allocate(wt_tmp(NSPINOR))

      ! transform KSHIFT into cartessian coordinates
      if( any( KSHIFT /= 0 ) ) then
        KSHIFT_TMP(:) = matmul(KSHIFT(:),bmat)
        KSHIFT_TMP = KSHIFT_TMP / AUTOA
        write(IUO,'(/,3X,A,3f12.6,A)')"KSHIFT = ",KSHIFT_TMP(1:3), &
     &" (A^-1)"
      endif

      ! integer -> real (original MAT_P2S is integer)
      RMAT_P2S = real(MAT_P2S,DP)
      call invert_3x3(RMAT_P2S,inv_mat)
!     write(IUO,'(/,3X,A)') &
!    &"Transform matrix from supercell to primitive cell"
!     do i=1,3
!        write(IUO,'(3X,3f12.6)')inv_mat(1:3,i)
!     enddo

      ! reciprocal lattice vectors
      call LATTIC(CELL,RE_CELL)

      allocate( proj_wt(NBANDS,NKPTS,ISPIN))
      proj_wt(1:NBANDS,1:NKPTS,1:ISPIN) = 0.0_DP

      write(IUO,'(/,3X,A)')"Starting K-projection:"
      IREC = 2
      spin:  do isp=1,ISPIN

      kpoint: do ik=1,NKPTS
            ! read WAVECAR, get the number of plane waves, k vector
            ! eigenvalue
            !
            IREC=IREC+1

            READ(IUW,REC=IREC,ERR=320) RNPL,VKPT, &
     & (CELTOT(J,ik,isp),IMEN,FERTOT(J,ik,isp),J=1,NBANDS)
            KPTVEC(1:3,ik)=VKPT(1:3)
            NPL=NINT(RNPL)

            NPLG = NPL/NSPINOR
            ! allocate array for G vectors and wave functions
            allocate(GVECT(3,NPLG),wt_proj(NPLG,NSPINOR),WFC(NPLG,NSPINOR))

            wt_proj = 1.0_DP

            ! generate G vemod_kproj_qectors
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
            ! set weights for k-projected bands
            do ip = 1,NPLG

               ! find out G vectors contributing to the primitive cell
               ! this is the k-projection
               !
               GVECTOR = 0
               GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)
               eta(:)=matmul(GVECTOR(:),inv_mat)
               !
               if( any( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
                 wt_proj(ip,1:NSPINOR) = 0.0_DP
               endif
            enddo

            allocate(CW2(NPL),CRD(NPL))

            band: do ib = 1,NBANDS
            
              ! get coefficients of wave functions from WAVECAR
              IREC=IREC+1

              NORM_TEST = 0.0_DP
              IF (SINGLE_PREC) THEN
                READ(IUW,REC=IREC,ERR=340) (CRD(i),i=1,NPL)
                CW2(1:NPL)=CRD(1:NPL)
              ELSE
                READ(IUW,REC=IREC,ERR=340) (CW2(i),i=1,NPL)
              ENDIF
              ! normalization
              !
              DO j=1,NSPINOR
              DO i=1,NPLG
                WFC(i,j) = CW2(i+(j-1)*NPLG)
                NORM_TEST = NORM_TEST + WFC(i,j)*conjg(WFC(i,j))
              ENDDO
              ENDDO
!             write(IUO,'(7X,A,I4,A,f12.6)')"BAND:",ib, &
!    & "  Normalization factor:",NORM_TEST
              NORM_FACTOR(ib,ik,isp) = NORM_TEST
            
              ! weight coefficients of WFC by k-projection
              do ip = 1,NPLG
                !
                WFC(ip,1:NSPINOR)=WFC(ip,1:NSPINOR)*wt_proj(ip,1:NSPINOR)
                !
              enddo

              wt_tmp(1:NSPINOR) = 0.0_DP
              DO j=1,NSPINOR
              DO i=1,NPLG
                 wt_tmp(j) = wt_tmp(j) + WFC(i,j)*conjg(WFC(i,j))
              ENDDO
              ENDDO
              !
              DO j=1,NSPINOR
                proj_wt(ib,ik,isp) = proj_wt(ib,ik,isp)+wt_tmp(j)
              ENDDO
              proj_wt(ib,ik,isp) = proj_wt(ib,ik,isp)/NORM_FACTOR(ib,ik,isp)
             
            enddo band    
            !
            write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)
            !
            deallocate(GVECT,wt_proj,WFC)
            deallocate(CW2,CRD)
            !
      enddo kpoint
      enddo spin
      end if
      return

  320 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues fromWAVECAR'
      STOP

  340 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP
      end subroutine bands_wt_all

      subroutine bands_wt_layer

!     calculates the weight of states projected on to the unit cell for a selected spatial window.
!     2013 CMX 

      use mod_input
      use mod_wave
      use mod_fft
      use mod_lattice
      use omp_lib
      use, intrinsic :: iso_c_binding
      include "fftw3.f03"
!              fft mesh arrays
      complex(DP),allocatable :: q_tmp(:,:,:),q_pw(:,:,:,:)
      complex(DP),allocatable :: u_step_z(:),u_step_y(:),u_step_x(:)
      complex(DP),allocatable :: q_pw_x(:,:,:,:),q_pw_y(:,:,:,:),q_pw_z(:,:,:,:)
      complex(DP),allocatable :: wf_x(:,:),wf_y(:,:),wf_z(:,:)
      complex(DP)             :: proj_wt_x,proj_wt_y,proj_wt_z,proj_wt_t

!              k-projected information
      real(DP),allocatable    :: wt_proj(:,:)
      real(DP)                :: inv_mat(3,3),RMAT_P2S(3,3)
      real(DP)                :: eta(3)

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
      real(DP)                :: EPSL2,RT_EPSL2,KPG(3),KPG_CART(3)
      real(DP)                :: RE_CELL(3,3),KSHIFT_TMP(3)
      !------FFTW
      integer(C_INT) :: N_DIM1,N_DIM2,N_DIM3,erro
      type(C_PTR) :: plan1,plan2
      complex(C_DOUBLE_COMPLEX), pointer :: q_pw_fftw(:,:,:,:)
      type(C_PTR) :: fftw_allocate

      NSPINOR = 1
      IF(LSORBIT) NSPINOR = 2
      error =  fftw_init_threads()
      call fftw_plan_with_nthreads(omp_get_max_threads())
      call fft_allocate(q_pw, nfft, NSPINOR )  !only get nfft(:)

      allocate( q_tmp( nfft(1), nfft(2), nfft(3)) )
      q_tmp = cmplx( 0.0_DP, 0.0_DP, DP )

      N_DIM1=nfft(1)
      N_DIM2=nfft(2)
      N_DIM3=nfft(3)
      
      fftw_allocate = fftw_alloc_complex(int(N_DIM1*N_DIM2*N_DIM3*NSPINOR, C_SIZE_T))
      call c_f_pointer(fftw_allocate, q_pw_fftw,[N_DIM1,N_DIM2,N_DIM3,NSPINOR])

      plan1 = fftw_plan_dft_3d(N_DIM1,N_DIM2,N_DIM3,&
      q_pw_fftw(:,:,:,1),q_pw_fftw(:,:,:,1),1,FFTW_MEASURE)
      plan2 = fftw_plan_dft_3d(N_DIM1,N_DIM2,N_DIM3,&
      q_pw_fftw(:,:,:,1),q_pw_fftw(:,:,:,1),-1,FFTW_MEASURE)

!--------------------------------------------------------------------------------
      IF(exqe==.true.) THEN
      call bands_wt_layer_qe
      ELSE IF(exab==.true.) THEN
      call bands_wt_layer_abinit
      ELSE
!-------------------------------------------------------------------------------
      ! transform KSHIFT into cartessian coordinates
      if( any( KSHIFT /= 0 ) ) then
        KSHIFT_TMP(:) = matmul(KSHIFT(:),bmat)
        KSHIFT_TMP = KSHIFT_TMP / AUTOA
        write(IUO,'(/,3X,A,3f12.6,A)')"KSHIFT = ",KSHIFT_TMP(1:3), &
     &" (A^-1)"
      endif

      ! integer -> real (original MAT_P2S is integer)
      RMAT_P2S = real(MAT_P2S,DP)
      call invert_3x3(RMAT_P2S,inv_mat)
!     write(IUO,'(/,3X,A)') &
!    &"Transform matrix from supercell to primitive cell"
!     do i=1,3
!        write(IUO,'(3X,3f12.6)')inv_mat(1:3,i)
!     enddo

      ! reciprocal lattice vectors
      call LATTIC(CELL,RE_CELL)

      ! step function
      allocate( u_step_z(nfft(3)),u_step_y(nfft(2)),u_step_x(nfft(1))) 

      if(LXLAYER.and..not.LYLAYER.and..not.LZLAYER)then
        !
        call step_xlayer(u_step_x ,nfft, xlay1, xlay2)
        !
      elseif(.not.LXLAYER.and.LYLAYER.and..not.LZLAYER)then
        !
        call step_ylayer(u_step_y ,nfft, ylay1, ylay2)
        !
      elseif(.not.LXLAYER.and..not.LYLAYER.and.LZLAYER)then
        !
        call step_zlayer(u_step_z ,nfft, zlay1, zlay2)
        !
      elseif(LXLAYER.and.LYLAYER.and..not.LZLAYER)then
        !
        call step_xlayer(u_step_x ,nfft, xlay1, xlay2)
        !
        call step_ylayer(u_step_y ,nfft, ylay1, ylay2)
        !
      elseif(LXLAYER.and..not.LYLAYER.and.LZLAYER)then
        !
        call step_xlayer(u_step_x ,nfft, xlay1, xlay2)
        !
        call step_zlayer(u_step_z ,nfft, zlay1, zlay2)
        !
      elseif(.not.LXLAYER.and.LYLAYER.and.LZLAYER)then
        !
        call step_ylayer(u_step_y ,nfft, ylay1, ylay2)
        !
        call step_zlayer(u_step_z ,nfft, zlay1, zlay2)
        !
      elseif(LXLAYER.and.LYLAYER.and.LZLAYER)then
        !
        call step_xlayer(u_step_x ,nfft, xlay1, xlay2)
        !
        call step_ylayer(u_step_y ,nfft, ylay1, ylay2)
        !
        call step_zlayer(u_step_z ,nfft, zlay1, zlay2)
        !
      endif

!      call fft_allocate( q_pw, nfft, NSPINOR )
!      allocate( q_tmp( nfft(1), nfft(2), nfft(3)) )
!      q_tmp = cmplx( 0.0_DP, 0.0_DP, DP )

      ! for polarization effect 
      if(LPOLA)then
        call fft_allocate( q_pw_x, nfft, NSPINOR ) 
        call fft_allocate( q_pw_y, nfft, NSPINOR ) 
        call fft_allocate( q_pw_z, nfft, NSPINOR ) 
      endif

      allocate( proj_wt(NBANDS,NKPTS,ISPIN))

      write(IUO,'(/,3X,A)')"Starting K-projection:"

      IREC = 2
      spin:  do isp=1,ISPIN

      kpoint: do ik=1,NKPTS

            ! read WAVECAR, get the number of plane waves, k vector
            ! eigenvalue
            !
            IREC=IREC+1

            READ(IUW,REC=IREC,ERR=220) RNPL,VKPT, &
     & (CELTOT(J,ik,isp),IMEN,FERTOT(J,ik,isp),J=1,NBANDS)
            KPTVEC(1:3,ik)=VKPT(1:3)
            NPL=NINT(RNPL)

            NPLG = NPL/NSPINOR

            ! allocate array for G vectors and wave functions
            allocate(GVECT(3,NPLG),wt_proj(NPLG,NSPINOR),WFC(NPLG,NSPINOR))
            allocate(wf_x(NPLG,NSPINOR),wf_y(NPLG,NSPINOR),wf_z(NPLG,NSPINOR))

            wt_proj = 1.0_DP

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

            ! set weights for k-projected bands
            do ip = 1,NPLG

               ! find out G vectors contributing to the primitive cell
               ! this is the k-projection
               !
               GVECTOR = 0
               GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)
               eta(:)=matmul(GVECTOR(:),inv_mat)
               !
               if( any( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
                 wt_proj(ip,1:NSPINOR) = 0.0_DP
               endif
               !
            enddo

            allocate(CW2(NPL),CRD(NPL))

            band: do ib = 1,NBANDS
            
              ! get coefficients of wave functions from WAVECAR
              IREC=IREC+1

              NORM_TEST = 0.0_DP
              IF (SINGLE_PREC) THEN
                READ(IUW,REC=IREC,ERR=240) (CRD(i),i=1,NPL)
                CW2(1:NPL)=CRD(1:NPL)
              ELSE
                READ(IUW,REC=IREC,ERR=240) (CW2(i),i=1,NPL)
              ENDIF
              !
              ! normalization
              !
              DO j=1,NSPINOR
              DO i=1,NPLG
                 WFC(i,j) = CW2(i+(j-1)*NPLG)
                 NORM_TEST = NORM_TEST + WFC(i,j)*conjg(WFC(i,j))
              ENDDO
              ENDDO
!             write(IUO,'(7X,A,I4,A,f12.6)')"BAND:",ib, &
!    & "  Normalization factor:",NORM_TEST
              NORM_FACTOR(ib,ik,isp) = NORM_TEST
            
              ! weight coefficients of WFC by k-projection
              do ip = 1,NPLG
                !
                WFC(ip,1:NSPINOR)=WFC(ip,1:NSPINOR)*wt_proj(ip,1:NSPINOR)
                !
              enddo

               !
               ! do integrating over partial space (zlay1, zlay2) and (ylay1, ylay2) using FFT
               !
               ! using FFT to get the Bloch function u from G space to real space 
               !
               do j=1,NSPINOR
               call wf_to_fftmesh(q_pw_fftw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
               ! call fft( q_pw(:,:,:,j), 'k->r', nfft ) 2022-5-6 jxchen
               !
               call fftw_execute_dft(plan1,q_pw_fftw(:,:,:,j),q_pw_fftw(:,:,:,j)) 
               !
               ! u* x u
               q_tmp(:,:,:) = q_pw_fftw(:,:,:,j)
               q_tmp = q_tmp*conjg(q_tmp)/(nfft(1)*nfft(2)*nfft(3))
               q_pw_fftw(:,:,:,j) = q_tmp(:,:,:)
               !
               ! back FFT to put u* x u in G space
               !call fft( q_pw(:,:,:,j), 'r->g', nfft ) 2022-5-6 jxchen
               call fftw_execute_dft(plan2,q_pw_fftw(:,:,:,j),q_pw_fftw(:,:,:,j))

               ! q_pw(:,:,:,j)=q_pw(:,:,:,j)/(nfft(1)*nfft(2)*nfft(3))
               !
               enddo

               if(NSPINOR>1)q_pw_fftw(:,:,:,1) = q_pw_fftw(:,:,:,1) + q_pw_fftw(:,:,:,2)
 
               !
               if(LXLAYER.and..not.LYLAYER.and..not.LZLAYER)then
                !
                proj_wt(ib,ik,isp) = real(sum(q_pw_fftw(:,1,1,1)*u_step_x), DP)
                !
               elseif(.not.LXLAYER.and.LYLAYER.and..not.LZLAYER)then
                !
                proj_wt(ib,ik,isp) = real(sum(q_pw_fftw(1,:,1,1)*u_step_y), DP)
                !
               elseif(.not.LXLAYER.and..not.LYLAYER.and.LZLAYER)then
                !
                proj_wt(ib,ik,isp) = real(sum(q_pw_fftw(1,1,:,1)*u_step_z), DP)
                !
               elseif(LXLAYER.and.LYLAYER.and..not.LZLAYER)then
                !
                sum_stepf = 0.0_DP
                do i=1,nfft(1)
                  do j=1,nfft(2)
                 sum_stepf=sum_stepf+q_pw_fftw(i,j,1,1)*u_step_x(i)*u_step_y(j)
                  enddo
                enddo
                !
                proj_wt(ib,ik,isp) = real(sum_stepf, DP)
                !
               elseif(LXLAYER.and..not.LYLAYER.and.LZLAYER)then
                !
                sum_stepf = 0.0_DP
                do i=1,nfft(1)
                  do j=1,nfft(3)
                 sum_stepf=sum_stepf+q_pw_fftw(i,1,j,1)*u_step_x(i)*u_step_z(j)
                  enddo
                enddo
                !
                proj_wt(ib,ik,isp) = real(sum_stepf, DP)
                !
               elseif(.not.LXLAYER.and.LYLAYER.and.LZLAYER)then
                !
                sum_stepf = 0.0_DP
                do i=1,nfft(2)
                  do j=1,nfft(3)
                 sum_stepf=sum_stepf+q_pw_fftw(1,i,j,1)*u_step_y(i)*u_step_z(j)
                  enddo
                enddo
                !
                proj_wt(ib,ik,isp) = real(sum_stepf, DP)
                !
               elseif(LXLAYER.and.LYLAYER.and.LZLAYER)then
                !
                sum_stepf = 0.0_DP
                do i=1,nfft(1)
                  do j=1,nfft(2)
                   do k=1,nfft(3)
                      sum_stepf=sum_stepf+ &
     & q_pw_fftw(i,j,k,1)*u_step_x(i)*u_step_y(j)*u_step_z(k)
                   enddo
                  enddo
                enddo
                !
                proj_wt(ib,ik,isp) = real(sum_stepf, DP)
                !
               endif
               
               proj_wt(ib,ik,isp)=proj_wt(ib,ik,isp)/NORM_FACTOR(ib,ik,isp)
            enddo band    
            !
            write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)
            !
            deallocate(GVECT,wt_proj,WFC,wf_x,wf_y,wf_z)
            deallocate(CW2,CRD)
            !
      enddo kpoint
      enddo spin
      !----------------------------------------------
      end if
      return

  220 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
      STOP

  240 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP
      end subroutine bands_wt_layer

       end module mod_kproj
