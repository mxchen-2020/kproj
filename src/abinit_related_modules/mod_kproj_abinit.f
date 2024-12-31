MODULE mod_kproj_abinit
   
    !interface to abinit 
    !06/2021 jxchen
   
    USE mod_comp
    USE mod_input
    USE mod_fft
    USE mod_lattice
    USE mod_wave
    use mod_step_funciton
     contains

     SUBROUTINE bands_wt_all_abinit
     INTEGER :: NPLG,NSPINOR
     INTEGER,allocatable::GVECT(:,:)
     COMPLEX(DP),allocatable::CW2(:)
     integer :: isp,ik,ib,ip
     !kprojected information
     REAL(DP),allocatable       :: wt_tmp(:),wt_proj(:,:)
     REAL(DP)                   :: inv_mat(3,3),RMAT_P2S(3,3),eta(3),NORM_TEST
     REAL(DP)                   :: RE_CELL(3,3),KSHIFT_TMP(3)
     INTEGER                    :: GVECTOR(3)
     COMPLEX(DP),allocatable    :: WFC(:,:)
         
        if( any( KSHIFT /= 0 ) ) then
          KSHIFT_TMP(:) = matmul(KSHIFT(:),bmat)
          KSHIFT_TMP = KSHIFT_TMP / AUTOA
          write(IUO,'(/,3X,A,3f12.6,A)')"KSHIFT = ",KSHIFT_TMP(1:3), &
          &" (A^-1)"
        endif  
        RMAT_P2S = real(MAT_P2S,DP)
        call invert_3x3(RMAT_P2S,inv_mat)
        call LATTIC(CELL,RE_CELL)
        IF (NRSPINORS == 2) THEN
          NSPINOR = 2
        ELSE
          NSPINOR = 1
        ENDIF
        allocate(proj_wt(NBANDS,NKPTS,ISPIN))
        allocate(wt_tmp(NSPINOR))
        proj_wt = 0.0_DP
        write(IUO,'(/,3X,A)')"Starting K-projection:"
  
       KPTVEC=kpt_ab

       do isp=1,ISPIN
        do ik=1,NKPTS
          READ(IUW) NPLG,NSPINOR,NBANDS

          ALLOCATE(GVECT(3,NPLG))
          ALLOCATE(CW2(NPLG*NSPINOR))
          ALLOCATE (wt_proj(NPLG,NSPINOR))
          wt_proj = 1.0_DP
          ALLOCATE(WFC(NPLG,NSPINOR)) 

          READ(IUW) GVECT
          ! set weights for k-projected bands
          call set_nsphk(NPLG,GVECT)

          do ip = 1,NPLG
            ! find out G vectors contributing to the primitive cell
            ! this is the k-projection
            GVECTOR = 0
            GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)
            eta(:)=matmul(GVECTOR(:),inv_mat)
            !
            if(any ( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
               wt_proj(ip,1:NSPINOR) = 0.0_DP
            endif
          enddo
          
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

          ! weight coefficients of WFC by k-projection
          DO ip = 1,NPLG
          !         
            WFC(ip,1:NSPINOR)=WFC(ip,1:NSPINOR)*wt_proj(ip,1:NSPINOR)
          !
          ENDDO
          !
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

          enddo  !band
          !
          write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)

          DEALLOCATE(GVECT)
          DEALLOCATE(CW2)
          DEALLOCATE(wt_proj)
          DEALLOCATE(WFC)
        enddo !kpt
      enddo  !spin
      close(IUW)
      END SUBROUTINE bands_wt_all_abinit

      SUBROUTINE bands_wt_layer_abinit
  
      use, intrinsic :: iso_c_binding
      use omp_lib
      include "fftw3.f03"
 
     !            fft mesh arrays
      complex(DP),allocatable :: q_pw(:,:,:,:),q_tmp(:,:,:)
      complex(DP),allocatable :: u_step_z(:),u_step_y(:),u_step_x(:)
      complex(DP),allocatable :: q_pw_x(:,:,:,:),q_pw_y(:,:,:,:),q_pw_z(:,:,:,:)
      complex(DP),allocatable :: wf_x(:,:),wf_y(:,:),wf_z(:,:)
      complex(DP)             :: proj_wt_x,proj_wt_y,proj_wt_z,proj_wt_t

     !            k-projected information
      real(DP),allocatable    :: wt_proj(:,:),eigenvalues(:,:),weight(:)
      real(DP)                :: inv_mat(3,3),RMAT_P2S(3,3)
      real(DP)                :: eta(3)
     !           abinit information
      complex(DP),allocatable :: WFC(:,:)
      integer,allocatable     :: GVECT(:,:)
      complex(DP),allocatable :: CW2(:)
      real(DP)                :: NORM_TEST

      integer                 :: i,j,k,l,n,ib,isp,ik,ip,ic,ios
      integer                 :: NSPINOR,NPLG,GVECTOR(3)
      real(DP)                :: EPSL2,RT_EPSL2,KPG(3),KPG_CART(3)
      real(DP)                :: RE_CELL(3,3),KSHIFT_TMP(3)
!fftw
      integer(C_INT) :: N_DIM1,N_DIM2,N_DIM3,erro
      type(C_PTR) :: plan1,plan2

       write(IUO,'(/,3X,A)')"Projection over the selected spatial window"

       IF (NRSPINORS == 2) THEN
          NSPINOR = 2
       ELSE
          NSPINOR = 1
       ENDIF

      N_DIM1=nfft(1)
      N_DIM2=nfft(2)
      N_DIM3=nfft(3)

      error =  fftw_init_threads()
      call fftw_plan_with_nthreads(omp_get_max_threads())
      call fft_allocate(q_pw, nfft, NSPINOR )

      plan1 = fftw_plan_dft_3d(N_DIM1,N_DIM2,N_DIM3,&
      q_pw(:,:,:,1),q_pw(:,:,:,1),1,FFTW_MEASURE)
      plan2 = fftw_plan_dft_3d(N_DIM1,N_DIM2,N_DIM3,&
      q_pw(:,:,:,1),q_pw(:,:,:,1),-1,FFTW_MEASURE)



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

     ! call fft_allocate( q_pw, nfft, NSPINOR )
      allocate( q_tmp( nfft(1), nfft(2), nfft(3)) )
      q_tmp = cmplx( 0.0_DP, 0.0_DP, DP )

      ! for polarization effect
      if(LPOLA)then
        call fft_allocate( q_pw_x, nfft, NSPINOR )
        call fft_allocate( q_pw_y, nfft, NSPINOR )
        call fft_allocate( q_pw_z, nfft, NSPINOR )
      endif

      allocate( proj_wt(NBANDS,NKPTS,ISPIN))

      write(IUO,'(/,3X,A)')"Starting K-projection:"

       KPTVEC=kpt_ab

       do isp=1,ISPIN
        do ik=1,NKPTS
          READ(IUW) NPLG,NSPINOR,NBANDS

          ALLOCATE(GVECT(3,NPLG))
          ALLOCATE (wf_x(NPLG,NSPINOR),wf_y(NPLG,NSPINOR),wf_z(NPLG,NSPINOR))
          ALLOCATE(CW2(NPLG*NSPINOR))
          ALLOCATE (wt_proj(NPLG,NSPINOR))
          wt_proj = 1.0_DP
          ALLOCATE(WFC(NPLG,NSPINOR))

          READ(IUW) GVECT
          ! set weights for k-projected bands
          call set_nsphk(NPLG,GVECT)

          do ip = 1,NPLG
            ! find out G vectors contributing to the primitive cell
            ! this is the k-projection
            GVECTOR = 0
            GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)
            eta(:)=matmul(GVECTOR(:),inv_mat)
            !
            if(any ( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
               wt_proj(ip,1:NSPINOR) = 0.0_DP
            endif
          enddo

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

              ! weight coefficients of WFC by k-projection
          DO ip = 1,NPLG
          !
            WFC(ip,1:NSPINOR)=WFC(ip,1:NSPINOR)*wt_proj(ip,1:NSPINOR)
          !
          ENDDO
          !

         !
          ! do integrating over partial space (zlay1, zlay2) and (ylay1, ylay2) using FFT
          !
          ! using FFT to get the Bloch function u from G space to real space
          !
          DO j=1,NSPINOR
          call wf_to_fftmesh( q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
        !  call fft( q_pw(:,:,:,j), 'k->r', nfft )
        !  !
        !  ! u* x u
        !  q_tmp(:,:,:) = q_pw(:,:,:,j)
        !  q_tmp = q_tmp*conjg(q_tmp)
        !  q_pw(:,:,:,j) = q_tmp(:,:,:)
        !  !
        !  ! back FFT to put u* x u in G space
        !  call fft( q_pw(:,:,:,j), 'r->g', nfft )
          !
          call fftw_execute_dft(plan1,q_pw(:,:,:,j),q_pw(:,:,:,j)) 
           !
           ! u* x u
           q_tmp(:,:,:) = q_pw(:,:,:,j)
           q_tmp = q_tmp*conjg(q_tmp)/(nfft(1)*nfft(2)*nfft(3))
           q_pw(:,:,:,j) = q_tmp(:,:,:)
           !
           ! back FFT to put u* x u in G space
           call fftw_execute_dft(plan2,q_pw(:,:,:,j),q_pw(:,:,:,j))
          enddo

         if(NSPINOR>1)q_pw(:,:,:,1) = q_pw(:,:,:,1) + q_pw(:,:,:,2)

              !
              if(LXLAYER.and..not.LYLAYER.and..not.LZLAYER)then
               !
               proj_wt(ib,ik,isp) = real(sum(q_pw(:,1,1,1)*u_step_x), DP)
               !
              elseif(.not.LXLAYER.and.LYLAYER.and..not.LZLAYER)then
               !
               proj_wt(ib,ik,isp) = real(sum(q_pw(1,:,1,1)*u_step_y), DP)
               !
              elseif(.not.LXLAYER.and..not.LYLAYER.and.LZLAYER)then
               !
               proj_wt(ib,ik,isp) = real(sum(q_pw(1,1,:,1)*u_step_z), DP)
               !
              elseif(LXLAYER.and.LYLAYER.and..not.LZLAYER)then
               !
               sum_stepf = 0.0_DP
               do i=1,nfft(1)
                 do j=1,nfft(2)
                sum_stepf=sum_stepf+q_pw(i,j,1,1)*u_step_x(i)*u_step_y(j)
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
                sum_stepf=sum_stepf+q_pw(i,1,j,1)*u_step_x(i)*u_step_z(j)
                 enddo
               enddo
              proj_wt(ib,ik,isp) = real(sum_stepf, DP)
               !
              elseif(.not.LXLAYER.and.LYLAYER.and.LZLAYER)then
               !
               sum_stepf = 0.0_DP
               do i=1,nfft(2)
                 do j=1,nfft(3)
                sum_stepf=sum_stepf+q_pw(1,i,j,1)*u_step_y(i)*u_step_z(j)
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
          & q_pw(i,j,k,1)*u_step_x(i)*u_step_y(j)*u_step_z(k)
                  enddo
                 enddo
               enddo
               !
               proj_wt(ib,ik,isp) = real(sum_stepf, DP)
               !
              endif

              proj_wt(ib,ik,isp)=proj_wt(ib,ik,isp)/NORM_FACTOR(ib,ik,isp)
           enddo !band
           !
           write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)
           !
           deallocate(GVECT,wt_proj,WFC,wf_x,wf_y,wf_z)
           deallocate(CW2)
           !
      enddo! kpoint
      enddo! spin
      close(IUW)
      END SUBROUTINE bands_wt_layer_abinit
 
end module mod_kproj_abinit

      
      
      
      
      
      
      
      
      
      
      
























    
