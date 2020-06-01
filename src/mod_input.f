      module mod_input
       use mod_comp

       character(5),public:: task
       integer,public     :: MAT_P2S(3,3)   ! matrix that transforms Primitive Cell -> Supercell
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

       contains
     
       subroutine input
         logical      :: LDUM, lopen,ex
         character    :: CHARAC
         integer      :: IDUM, IERR, N, II
         real(DP)     :: RDUM,fdum
         complex(DPC) :: CDUM

         inquire (file='INKPROJ', exist=ex)
         if (ex.eq..false.)  & 
     & call fstop(IUO, 'Error: file = INKPROJ is not found, ABORTING!')

         open(IUI,file='INKPROJ',status='old')
         write(IUO,'(/,3x,A)')"Read in control paramters"

         lopen = .FALSE.

         LZLAYER = .FALSE.
         call RDATAB(LOPEN,'INKPROJ',IUI,'LZLAYER','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LZLAYER,CHARAC,N,1,IERR)
         if (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
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
     &        IDUM, zlay2, CDUM, LDUM, CHARAC, N, 1, IERR )
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
            write(IUO,*)'We exchange zlay1 and zlay2'
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
    
         if (ylay2 < ylay1) then
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
     & 'Error reading item ''LXLAYER'' from file INKPROJ.'
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

      end module mod_input
