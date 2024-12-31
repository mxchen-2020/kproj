  module mod_kproj_qe
      
      !interface to qe 
      !06/2021 jxchen

      USE mod_comp
      USE read_information_from_qe
      USE mod_input
      USE mod_step_funciton
      USE mod_fft
      USE mod_lattice
      USE mod_wave
     
      contains

      subroutine bands_wt_all_qe
      REAL(DP),allocatable       :: wt_tmp(:),wt_proj(:,:),eigenvalues(:,:),weight(:)
      REAL(DP)                   :: inv_mat(3,3),RMAT_P2S(3,3)
      
      INTEGER                    :: ik,spin_direction,i,j,ip,isp,ikpt
      INTEGER                    :: ngw,GVECTOR(3) 
      LOGICAL                    :: gamma_only
      INTEGER                    :: LBANDS
      REAL(DP)                   :: VKPT(3),eta(3),NORM_TEST
      REAL(DP)                   :: b1_l(3), b2_l(3), b3_l(3)  
      REAL(DP)                   :: scalef,RE_CELL(3,3),KSHIFT_TMP(3)
      CHARACTER(len=100)         :: wfcname,ikp,spin !The length of waven.dat cannot be greater than len
      COMPLEX(DP),allocatable    :: WFC(:,:),CW2(:)
      INTEGER,allocatable        :: GVECT(:,:)
       allocate(weight(NKPTS))
       call read_kpoint_energy(KPTVEC,eigenvalues,weight)
        do ib=1,NBANDS
          do isp = 1,ISPIN
            do ik =1,NKPTS
            FERTOT(ib,ik,isp)=weight(ik)
            enddo
          enddo
        enddo
!---------------------------------------------------------------------------------
        if( any( KSHIFT /= 0 ) ) then
          KSHIFT_TMP(:) = matmul(KSHIFT(:),bmat) 
          KSHIFT_TMP = KSHIFT_TMP / AUTOA  
          write(IUO,'(/,3X,A,3f12.6,A)')"KSHIFT = ",KSHIFT_TMP(1:3), &
          &" (A^-1)"
        endif  
        RMAT_P2S = real(MAT_P2S,DP)  
        call invert_3x3(RMAT_P2S,inv_mat) 
        call LATTIC(CELL,RE_CELL)
        IF (noncolin == .true.) THEN
          NSPINOR = 2
        ELSE
          NSPINOR = 1
        ENDIF
        allocate(proj_wt(NBANDS,NKPTS,ISPIN))
        allocate(wt_tmp(NSPINOR))
        proj_wt = 0.0_DP
        write(IUO,'(/,3X,A)')"Starting K-projection:"
!----------------------------------------------------------------------------------
    do isp=1,ISPIN !spin
       ! 
      do ik=1,NKPTS!kpoint
        !
        IF (ISPIN==1) spin =""
        IF (isp==1 .and. ISPIN==2 ) spin ="up"
        IF (isp==2 .and. ISPIN==2 ) spin ="dw"
        write(ikp,*) ik
        wfcname = "wfc"//trim(spin)//trim(adjustl(ikp))//".dat" 
        OPEN ( UNIT = IUW, FILE =wfcname, &
             & FORM='unformatted', STATUS = 'unknown' )  
        !
        read(IUW) ikpt,(VKPT(i),i=1,3),spin_direction,gamma_only,scalef 
        read(IUW) ngw,NPLG,NSPINOR,LNBANDS   !当有自旋时候，这里NBANDS就是对应自旋的能带数
        ! 
        IF(isp == 1) THEN
        CELTOT(:,ik,isp)=eigenvalues(1:NBANDS,ik)
        ELSE IF(isp == 2 ) THEN
        CELTOT(:,ik,isp)=eigenvalues(NBANDS+1:2*NBANDS,ik)
        END IF
        read(IUW) b1_l,b2_l,b3_l
        !
        allocate (CW2(NPLG*NSPINOR)) 
        allocate (GVECT(3,NPLG)) 
        allocate (wt_proj(NPLG,NSPINOR))
        wt_proj = 1.0_DP
        allocate (WFC(NPLG,NSPINOR))
        read(IUW) GVECT
!-----------------------------------------------------------------------! 
        ! set weights for k-projected bands
        call set_nsphk(NPLG,GVECT) 
        do ip = 1,NPLG

          ! find out G vectors contributing to the primitive cell
          ! this is the k-projection
          !
            GVECTOR = 0
            GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)  
            eta(:)=matmul(GVECTOR(:),inv_mat) 
           !
           if(any ( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
               wt_proj(ip,1:NSPINOR) = 0.0_DP   
           endif
        enddo
!-----------------------------------------------------------------------------!
        do ib  = 1,NBANDS !band
          read(IUW) CW2    !CW2
!-----------------------------------------------------------------------------!
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
!----------------------------------------------------------------------------------!  
        ENDDO! band
        !
        write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)
        !
        deallocate (CW2)
        deallocate (GVECT)
        deallocate (wt_proj)
        deallocate (WFC)
      
       end do !kpoiint
     end do  !spin
     close(IUW)

      end subroutine bands_wt_all_qe
      
      subroutine bands_wt_layer_qe

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
      complex(DP),allocatable :: q_pw(:,:,:,:),q_tmp(:,:,:)
      complex(DP),allocatable :: u_step_z(:),u_step_y(:),u_step_x(:)
      complex(DP),allocatable :: q_pw_x(:,:,:,:),q_pw_y(:,:,:,:),q_pw_z(:,:,:,:)
      complex(DP),allocatable :: wf_x(:,:),wf_y(:,:),wf_z(:,:)
      complex(DP)             :: proj_wt_x,proj_wt_y,proj_wt_z,proj_wt_t

!              k-projected information
      real(DP),allocatable    :: wt_proj(:,:),eigenvalues(:,:),weight(:)
      real(DP)                :: inv_mat(3,3),RMAT_P2S(3,3)
      real(DP)                :: eta(3)

!              wfc.dat related
      LOGICAL                 :: gamma_only
      INTEGER                 :: spin_direction,ikpt
      complex(DP),allocatable :: WFC(:,:)
      integer,allocatable     :: GVECT(:,:)
      complex(DP),allocatable :: CW2(:)
      real(DP)                :: VKPT(3)
      real(DP)                :: NORM_TEST
      REAL(DP)                :: scalef
      integer                 :: GVECTOR(3)
      integer                 :: i,j,k,l,n,ib,isp,ik,ip,ic,ios
      integer                 :: NSPINOR,NPLG,ngw,LNANDS
      real(DP)                :: EPSL2,RT_EPSL2,KPG(3),KPG_CART(3)
      real(DP)                :: RE_CELL(3,3),KSHIFT_TMP(3)
      REAL(DP)                :: b1_l(3), b2_l(3), b3_l(3)
      complex(DP)             :: sum_stepf,sum_stepf_x,sum_stepf_y,sum_stepf_z
      CHARACTER(len=100)      :: wfcname,ikp,spin
!FFTW
      integer(C_INT) :: N_DIM1,N_DIM2,N_DIM3,erro
      type(C_PTR) :: plan1,plan2


      write(IUO,'(/,3X,A)')"Projection over the selected spatial window"
  
       allocate(weight(NKPTS))
       call read_kpoint_energy(KPTVEC,eigenvalues,weight)
       do ib=1,NBANDS
         do isp = 1,ISPIN
         FERTOT(ib,:,isp)=weight(:)
         enddo
       enddo
     
      NSPINOR = 1
      IF (noncolin == .true.) NSPINOR = 2      
      
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
        write(IpolarisUO,'(/,3X,A,3f12.6,A)')"KSHIFT = ",KSHIFT_TMP(1:3), &
        &" (A^-1)"
      endif  !找kshift
      ! integer -> real (original MAT_P2S is integer)
      RMAT_P2S = real(MAT_P2S,DP) 
      call invert_3x3(RMAT_P2S,inv_mat)
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

      write(IUO,'(/,3X,A)')"Starting K-projection:"
!-----------------------------------------------------------
        IF (noncolin == .true.) THEN
          NSPINOR = 2
        ELSE
          NSPINOR = 1
        END IF
        allocate(proj_wt(NBANDS,NKPTS,ISPIN))
      !  allocate(wt_tmp(NSPINOR))
        write(IUO,'(/,3X,A)')"Starting K-projection:"
!----------------------------------------------------------------------------------
 ! lsda = .true.
        if (lsda ==.true.) THEN
          ISPIN=2
        ELSE
          ISPIN=1
        END IF
   do isp=1,ISPIN !spin
      do ik  = 1,NKPTS !kpoint
        !
        IF (ISPIN==1) spin =""
        IF (isp==1 .and. ISPIN==2 ) spin ="up"
        IF (isp==2 .and. ISPIN==2 ) spin ="dw"
        write(ikp,*) ik
        wfcname = "wfc"//trim(spin)//trim(adjustl(ikp))//".dat"
        OPEN ( UNIT = IUW, FILE = wfcname , &
           &   FORM='unformatted', STATUS = 'unknown' )
        !
        read(IUW) ikpt,(VKPT(i),i=1,3),spin_direction,gamma_only,scalef
        read(IUW) ngw,NPLG,NSPINOR,LNBANDS   !当有自旋时候，这里NBANDS就是对应自旋的能带数
        ! 
        IF(isp == 1) THEN
        CELTOT(:,ik,isp)=eigenvalues(1:NBANDS,ik)
        ELSE IF(isp == 2 ) THEN
        CELTOT(:,ik,isp)=eigenvalues(NBANDS+1:2*NBANDS,ik)
        END IF
        read(IUW) b1_l,b2_l,b3_l
        !
        allocate (wf_x(NPLG,NSPINOR),wf_y(NPLG,NSPINOR),wf_z(NPLG,NSPINOR))
        allocate (CW2(NPLG*NSPINOR))
        allocate (GVECT(3,NPLG))
        allocate (wt_proj(NPLG,NSPINOR))
        allocate (WFC(NPLG,NSPINOR))
        wt_proj = 1.0_DP
        read(IUW) GVECT                
!-----------------------------------------------------------------------! 
        ! set weights for k-projected bands
        call set_nsphk(NPLG,GVECT)
        do ip = 1,NPLG

          ! find out G vectors contributing to the primitive cell
          ! this is the k-projection
          !
            GVECTOR = 0
            GVECTOR(:)=GVECT(:,ip)+KSHIFT(:)  
            eta(:)=matmul(GVECTOR(:),inv_mat) 
          !
           if(any ( abs( eta-anint(eta) ) > 1.e-04_DP ) ) then
               wt_proj(ip,1:NSPINOR) = 0.0_DP   
           endif
        enddo
!-----------------------------------------------------------------------------!
        do ib = 1,NBANDS !band
          read(IUW) CW2    !CW2
!-----------------------------------------------------------------------------!
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
            WFC(ip,1:NSPINOR)=WFC(ip,1:NSPINOR)*wt_proj(ip,1:NSPINOR)
          ENDDO
          !
          ! do integrating over partial space (zlay1, zlay2) and (ylay1, ylay2) using FFT
          !
          ! using FFT to get the Bloch function u from G space to real space
          !
          DO j=1,NSPINOR
          call wf_to_fftmesh( q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
         ! call fft( q_pw(:,:,:,j), 'k->r', nfft )
         ! !
         ! ! u* x u
         ! q_tmp(:,:,:) = q_pw(:,:,:,j)
         ! q_tmp = q_tmp*conjg(q_tmp)
         ! q_pw(:,:,:,j) = q_tmp(:,:,:)
         ! !
         ! ! back FFT to put u* x u in G space
         ! call fft( q_pw(:,:,:,j), 'r->g', nfft )
         !
         call fftw_execute_dft(plan1,q_pw(:,:,:,j),q_pw(:,:,:,j))  
         !
           ! u* x u
           q_tmp(:,:,:) = q_pw(:,:,:,j)
           q_tmp = q_tmp*conjg(q_tmp)/(nfft(1)*nfft(2)*nfft(3))
           q_pw(:,:,:,j) = q_tmp(:,:,:)
           !
           ! back FFT to put u* x u in G space
           !call fft( q_pw(:,:,:,j), 'r->g', nfft )
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
      return

      220 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
      STOP

      240 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP

      end subroutine bands_wt_layer_qe

end module mod_kproj_qe







