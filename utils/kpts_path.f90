      program kpts_path
        integer, parameter  :: DP =SELECTED_real_KIND(10)
        real(DP)            :: ki(3),kf(3),dk(3)
        integer             :: i,nkpts,ik
        real(DP),allocatable:: kvec(:,:)
        
        ki = 0.0_DP; kf = 0.0_DP

        write(*,'(5X,A)')"Please type ki(1),ki(2),ki(3) and &
     & kf(1),kf(2),kf(3)"
        read(*,*)ki(1:3),kf(1:3)

        write(*,'(5X,A)')"How many kpoints do you want along the path?"
        read(*,*)nkpts

        allocate(kvec(3,nkpts))
        kvec = 0.0_DP

        do i=1,3
          dk(i) = (kf(i) - ki(i)) / (nkpts-1)
        enddo

        open(1,file='kpts.dat',status='unknown')
        write(1,'(A)')"KPATH_generation"
        write(1,'(I8)')nkpts
        write(1,'(A)')"Direct"
        do ik = 1, nkpts
          kvec(1,ik) = ki(1) + (ik-1)*dk(1)
          kvec(2,ik) = ki(2) + (ik-1)*dk(2)
          kvec(3,ik) = ki(3) + (ik-1)*dk(3)
          write(1,'(3f20.14,I10)')kvec(1:3,ik),1
        enddo
        
        close(1) 
      end 
