      module mod_kproj
        use mod_comp
        real(DP),allocatable,public  :: proj_wt(:,:,:)

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

      return

  320 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
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

!              fft mesh arrays
      complex(DP),allocatable :: q_pw(:,:,:,:),q_tmp(:,:,:)
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

      complex(DP)      :: sum_stepf,sum_stepf_x,sum_stepf_y,sum_stepf_z

      write(IUO,'(/,3X,A)')"Projection over the selected spatial window"

      NSPINOR = 1
      IF(LSORBIT) NSPINOR = 2

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

      call fft_allocate( q_pw, nfft, NSPINOR ) 
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
               call wf_to_fftmesh( q_pw(:,:,:,j),nfft,NPLG,GVECT,WFC(:,j))
               call fft( q_pw(:,:,:,j), 'k->r', nfft )
               !
               ! u* x u
               q_tmp(:,:,:) = q_pw(:,:,:,j)
               q_tmp = q_tmp*conjg(q_tmp)
               q_pw(:,:,:,j) = q_tmp(:,:,:)
               !
               ! back FFT to put u* x u in G space
               call fft( q_pw(:,:,:,j), 'r->g', nfft )
               !
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
                !
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
            enddo band    
            !
            write(IUO,'(5X,A,I4,3f14.8)')"k-projection is done for kpt",ik,KPTVEC(1:3,ik)
            !
            deallocate(GVECT,wt_proj,WFC,wf_x,wf_y,wf_z)
            deallocate(CW2,CRD)
            !
      enddo kpoint
      enddo spin

      return

  220 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading eigenvalues from WAVECAR'
      STOP

  240 CONTINUE
      WRITE(IUO,'(5X,A)')'ERR: while reading coef. from WAVECAR'
      STOP

      end subroutine bands_wt_layer

      subroutine step_zlayer (f,nfft_tmp,zlay1_tmp,zlay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between zlay1 and zlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: zlay1_tmp,zlay2_tmp

      integer      :: k3,k3m,nf3
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k3m = min( k3d, (nfft_tmp(3)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k3=-k3m,k3m
         nf3 = modulo( k3+nfft_tmp(3), nfft_tmp(3) ) + 1

!--->    layer part: k1=0, k2=0
         if( k3 == 0 ) then
            ff = (zlay2_tmp - zlay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k3*zlay1_tmp
            arg2 = TWOPI*k3*zlay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k3*TWOPI)
         endif

         f(nf3) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine 

      subroutine step_ylayer (f,nfft_tmp,ylay1_tmp,ylay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between zlay1 and zlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: ylay1_tmp,ylay2_tmp

      integer      :: k2,k2m,nf2
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k2m = min( k2d, (nfft_tmp(2)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k2=-k2m,k2m
         nf2 = modulo( k2+nfft_tmp(2), nfft_tmp(2) ) + 1

!--->    layer part: k1=0, k2=0
         if( k2 == 0 ) then
            ff = (ylay2_tmp - ylay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k2*ylay1_tmp
            arg2 = TWOPI*k2*ylay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k2*TWOPI)
         endif

         f(nf2) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine 

      subroutine step_xlayer (f,nfft_tmp,xlay1_tmp,xlay2_tmp)

!     generates the complex values of the (complex conjugate of the)
!     step function for a layer between xlay1 and xlay2 on the full
!     reciprocal space fft mesh

      use mod_fft,      only : k1d,k2d,k3d

      implicit none

      integer,     intent(in)     :: nfft_tmp(3)
      complex(DP), intent(inout)  :: f(:)
      real(DP),    intent(in)     :: xlay1_tmp,xlay2_tmp

      integer      :: k1,k1m,nf1
      real(DP)     :: arg1,arg2
      complex(DP)  :: ff


!--->    fft mesh might be smaller than the total star list since,
!--->    for example, 2*kmax < gmax

      k1m = min( k1d, (nfft_tmp(1)-1)/2 )

      f = cmplx( 0.0_DP, 0.0_DP, DP )

      do k1=-k1m,k1m
         nf1 = modulo( k1+nfft_tmp(1), nfft_tmp(1) ) + 1

!--->    layer part: k1=0, k2=0
         if( k1 == 0 ) then
            ff = (xlay2_tmp - xlay1_tmp)  ! in lattice units
         else
            arg1 = TWOPI*k1*xlay1_tmp
            arg2 = TWOPI*k1*xlay2_tmp
            ff=cmplx(sin(arg2)-sin(arg1),cos(arg1)-cos(arg2),DP)/(k1*TWOPI)
         endif

         f(nf1) = ff   ! note: conjugate of step function

      enddo

      return
      end subroutine 

      subroutine invert_3x3(mat1,inv_mat1)
         !!! invert of matrix
         REAL(DP)        :: mat1(3,3),inv_mat1(3,3),det

         inv_mat1(1,1)=mat1(2,2)*mat1(3,3)-mat1(3,2)*mat1(2,3)
         inv_mat1(1,2)=mat1(3,2)*mat1(1,3)-mat1(1,2)*mat1(3,3)
         inv_mat1(1,3)=mat1(1,2)*mat1(2,3)-mat1(2,2)*mat1(1,3)
         inv_mat1(2,1)=mat1(2,3)*mat1(3,1)-mat1(3,3)*mat1(2,1)
         inv_mat1(2,2)=mat1(3,3)*mat1(1,1)-mat1(1,3)*mat1(3,1)
         inv_mat1(2,3)=mat1(1,3)*mat1(2,1)-mat1(2,3)*mat1(1,1)
         inv_mat1(3,1)=mat1(2,1)*mat1(3,2)-mat1(3,1)*mat1(2,2)
         inv_mat1(3,2)=mat1(3,1)*mat1(1,2)-mat1(1,1)*mat1(3,2)
         inv_mat1(3,3)=mat1(1,1)*mat1(2,2)-mat1(2,1)*mat1(1,2)
         call determ(mat1,det)
         inv_mat1(1:3,1:3)=inv_mat1(1:3,1:3)/det
       end subroutine invert_3x3

       subroutine determ(a,b)
       ! returnerer determinaten af a 3x3 matrix
       REAL(DP),intent(in)  :: a(3,3)
       real(DP),intent(out) :: b
         b=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+ &
     & a(1,3)*a(3,2)*a(2,1)-a(1,3)*a(2,2)*a(3,1)- &
     & a(2,3)*a(3,2)*a(1,1)-a(3,3)* &
         a(1,2)*a(2,1)
       END subroutine determ
        !
        !-------------------------------------------------------------
        ! reciprocal representation to cartessian representation
        !-------------------------------------------------------------
        !
        subroutine dire2car(A,B,re_cell)
          real(DP)     :: A(3), B(3), re_cell(3,3)
          B(1)=A(1)*re_cell(1,1)+A(2)*re_cell(1,2)+A(3)*re_cell(1,3)
          B(2)=A(1)*re_cell(2,1)+A(2)*re_cell(2,2)+A(3)*re_cell(2,3)
          B(3)=A(1)*re_cell(3,1)+A(2)*re_cell(3,2)+A(3)*re_cell(3,3)
        end subroutine dire2car

         SUBROUTINE LATTIC(A,B)
           real(DP)   ::  A(3,3),B(3,3),ANORM(3),BNORM(3),OMEGA

           CALL EXPRO(B(1,1),A(1,2),A(1,3))
           CALL EXPRO(B(1,2),A(1,3),A(1,1))
           CALL EXPRO(B(1,3),A(1,1),A(1,2))

           OMEGA =ABS(B(1,1)*A(1,1)+B(2,1)*A(2,1)+B(3,1)*A(3,1))

           do I=1,3
             do J=1,3
                B(I,J)=TWOPI*B(I,J)/OMEGA
             enddo
           enddo

         END SUBROUTINE LATTIC
         !
         !----------------- SUBROUTINE EXPRO   ---------------------------------
         ! EXPRO
         ! caclulates the x-product of two vectors
         !
         !----------------------------------------------------------------------
         !
         SUBROUTINE EXPRO(H,U1,U2)
           real(DP)    :: H(3),U1(3),U2(3)

           H(1)=U1(2)*U2(3)-U1(3)*U2(2)
           H(2)=U1(3)*U2(1)-U1(1)*U2(3)
           H(3)=U1(1)*U2(2)-U1(2)*U2(1)
         END SUBROUTINE EXPRO

      subroutine out_bands
      ! 
      use mod_input 
      use mod_wave
      use mod_lattice

      implicit none

!              wavecar file variables
      real(DP),allocatable :: e(:),occ(:)
      real(DP)             :: vkpt(3),rnv
      integer              :: nv

!              k-projected information
      real(DP),allocatable :: bkpts(:,:)

!              for HSE calculation
      real(DP),allocatable :: proj_wt_HSE(:,:,:),bkpts_HSE(:,:)
      real(DP),allocatable :: CELTOT_HSE(:,:,:),FERTOT_HSE(:,:,:)
      real(DP),allocatable :: KPTVEC_HSE(:,:)
      integer              :: NK_EFF ! the number of effctive kpoints

      logical,allocatable  :: l_good_kpt(:)

      integer              :: irec
      integer              :: i,k,l,n,nb,jspin,ikpt,ios,kpt,nkpts0
      real(DP)             :: emin_ef,emax_ef,ef

!     file name and handle for band data
      integer,   parameter :: bskfh = 73
      character(len=30)    :: bskfn = 'bs_projected.dat'

      character(len=256)   :: line

      integer, allocatable :: segment_endpoint(:),epts(:)
      real(DP),allocatable :: x_kpt(:)
      real(DP)             :: dk(3),dk_old(3),x0
      integer              :: kpt_lines,nseg,k0,max_segments

      if(LHSE)then
        NK_EFF = NKPTS - NK_HSE
        allocate(proj_wt_HSE(NBANDS,NK_EFF,ISPIN)) 
        allocate(CELTOT_HSE(NBANDS,NK_EFF,ISPIN)) 
        allocate(FERTOT_HSE(NBANDS,NK_EFF,ISPIN)) 
        allocate(KPTVEC_HSE(3,NK_EFF)) 

        do jspin=1,ISPIN
          do ikpt = 1,NK_EFF
            proj_wt_HSE(1:NBANDS,ikpt,jspin) = &
     &      proj_wt(1:NBANDS,ikpt+NK_HSE,jspin)

            CELTOT_HSE(1:NBANDS,ikpt,jspin) =  &
     &      CELTOT(1:NBANDS,ikpt+NK_HSE,jspin)

            FERTOT_HSE(1:NBANDS,ikpt,jspin) =  &
     &      FERTOT(1:NBANDS,ikpt+NK_HSE,jspin)

            KPTVEC_HSE(1:3,ikpt) = KPTVEC(1:3,ikpt+NK_HSE) 
          enddo
        enddo

        NKPTS = NK_EFF

        do jspin=1,ISPIN
          do ikpt = 1,NK_EFF
            proj_wt(1:NBANDS,ikpt,jspin) = &
     &      proj_wt_HSE(1:NBANDS,ikpt,jspin)

            CELTOT(1:NBANDS,ikpt,jspin) =  &
     &      CELTOT_HSE(1:NBANDS,ikpt,jspin)

            FERTOT(1:NBANDS,ikpt,jspin) =  &
     &      FERTOT_HSE(1:NBANDS,ikpt,jspin)

            KPTVEC(1:3,ikpt) = KPTVEC_HSE(1:3,ikpt)         
          enddo
        enddo
      endif
      
      allocate( e(NBANDS),occ(NBANDS))
      emin_ef = -1000
      emax_ef =  1000
      do jspin=1,ISPIN
        do ikpt=1,NKPTS
           e(1:NBANDS)   = CELTOT(1:NBANDS,ikpt,jspin)
           occ(1:NBANDS) = FERTOT(1:NBANDS,ikpt,jspin) 
!  get a better estimate of the fermi energy
         occ(1:nbands) = occ(1:nbands)/occ(1) ! assumes 1 occupied state
         emin_ef = max( maxval( e, mask = ( occ >= 0.5_DP ) ) , emin_ef )
         emax_ef = min( minval( e, mask = ( occ <  0.5_DP ) ) , emax_ef )
        enddo
      enddo
      ef = ( emin_ef + emax_ef )/2
      ef = 0.0_DP
      write(IUO,'(/,5X,"Fermi level obtained from bands = ",f0.5," eV")') ef
      write(IUO,'("    (need to change in bs_projected.dat)")')
        

!--->    break up the k-points into lines, assuming that the order is correct
      allocate( l_good_kpt(NKPTS) )
      allocate( bkpts(3,NKPTS) )
      bkpts = KPTVEC
      l_good_kpt = .true.
      max_segments = 8
      allocate( segment_endpoint(max_segments) )
      nseg = 0
      dk_old = bkpts(:,2) - bkpts(:,1)
      dk_old = dk_old/maxval( abs( dk_old ) ) ! normalize the direction
      do ikpt = 2,NKPTS-1
         dk = bkpts(:,ikpt+1) - bkpts(:,ikpt)
         if( all( abs( dk ) < 1.e-05_DP ) ) then ! vasp sometimes does
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
       write(IUO,'(/,5X,i0," repeated kpts; now total of ",i0/)') NKPTS-nkpts0,nkpts0
      endif

      allocate( x_kpt(NKPTS) )

      dk_old = bkpts(:,1)
      k0 = 1
      x0 = 0
      x_kpt(1) = x0
      kpt = 1
      do nseg = 1,kpt_lines
         do k = k0+1,segment_endpoint(nseg)
            if( .not. l_good_kpt(k) ) cycle
            kpt = kpt + 1
            dk = bkpts(:,k) - dk_old
            x_kpt(k) = x0 + sqrt( dot_product( dk, matmul( bbmat, dk ) ) )

            ! change to A^{-1}
            !x_kpt(k)=x_kpt(k)/AUTOA
         enddo
         x0 = x_kpt(segment_endpoint(nseg))
         k0 = segment_endpoint(nseg)
         dk_old = bkpts(:,segment_endpoint(nseg))
      enddo

!--->    write to bs_projected.dat

      open(bskfh,file=bskfn,form='formatted',status='replace')
      write(bskfh,'(i0,2x,i0,f20.8,"  ! NKPTS, jspins, ef (eV)")')      &
     &       nkpts0,ISPIN,ef
      write(bskfh,'(i0,"  ! number of segments; positions:")') kpt_lines
      write(bskfh,'(6f12.6)')                                           &
     & x_kpt( segment_endpoint(1:kpt_lines) )/AUTOA
      
      do jspin = 1,ISPIN
        kpt = 0
        do ikpt = 1,NKPTS

          if( .not. l_good_kpt(ikpt) ) cycle
          kpt = kpt + 1
          write(bskfh,'(i0,2x,i0,2x,i0,3f10.6,5x,f12.6)')               &
     &      kpt,NBANDS,jspin,bkpts(:,ikpt)+KSHIFT(:),x_kpt(ikpt)/AUTOA
          do nb = 1,NBANDS
            write(bskfh,'(f16.10,es22.5)')                              &
     &            CELTOT(nb,ikpt,jspin),proj_wt(nb,ikpt,jspin)
          enddo
        enddo
      enddo 

!     if(ISPIN==1)then
!       open(66,file="band_proj_gnu.dat",status='unknown')
!       do nb=1,NBANDS
!         do ikpt=1,NKPTS
!           if( .not. l_good_kpt(ikpt) ) cycle
!           write(66,'(3f16.8)')x_kpt(ikpt),&
!    & CELTOT(nb,ikpt,1),proj_wt(nb,ikpt,1)
!         enddo
!           write(66,*)
!       enddo
!       close(66)
!     else
!       open(66,file="band_proj_gnu_up.dat",status='unknown')
!       do nb=1,NBANDS
!         do ikpt=1,NKPTS
!           if( .not. l_good_kpt(ikpt) ) cycle
!           write(66,'(3f16.8)')x_kpt(ikpt),&
!    & CELTOT(nb,ikpt,1),proj_wt(nb,ikpt,1)
!         enddo
!           write(66,*)
!       enddo
!       close(66)

!       open(66,file="band_proj_gnu_dn.dat",status='unknown')
!       do nb=1,NBANDS
!         do ikpt=1,NKPTS
!           if( .not. l_good_kpt(ikpt) ) cycle
!           write(66,'(3f16.8)')x_kpt(ikpt),&
!    & CELTOT(nb,ikpt,2),proj_wt(nb,ikpt,2)
!         enddo
!           write(66,*)
!       enddo
!       close(66)
!     endif
      close( bskfh )

      end subroutine out_bands

       end module mod_kproj
