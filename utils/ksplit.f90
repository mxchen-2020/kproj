	program main
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!                                           !!
	!!  splitting KPOINTS in IBZKPT into sevral  !!
	!!  files for EIGENVALUE calculations.       !!
	!!                                           !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        character*40         :: filename
	character*3 file_name
	double precision,dimension(50000,3) :: kp
	integer,dimension(50000) :: wei
	integer nkp,m,subnkp,nfile
	!subnkp=13

        write(*,'(3X,A)')"Name of the file containing k-points"
        read(*,*)filename
        
        write(*,'(5X,A)')"No. of kpts for each job"
        read(*,*)subnkp

        open(1,file=filename,status='old')
	read(1,*)
	read(1,*)nkp
	read(1,*)
	do i=1,nkp
	read(1,99)kp(i,1:3),wei(i)
	enddo
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!                                       !!
        !!  how many files for kpoint splitting  !!
	!!                                       !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nfile=nkp/subnkp
        nfile=int(nfile)
	if (nkp-nfile*subnkp>1) then
	nfile=nfile+1
	else
	endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!                                       !!
 	!!        kpoint splitting               !!
	!!                                       !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do m=1,nfile
	if(m<10)then
	write(file_name,'(i1)')m
	else if(m.gt.9.and.m.lt.100)then
	write(file_name,'(i2)')m
	else
	write(file_name,'(i3)')m
	endif
	OPEN (unit=2,FILE="KPOINTS_"//file_name)
	write(2,'(A,I0)')"KPOINT_file_",m
!       if(nkp-(m-1)*subnkp.le.subnkp)then
        if(m==nfile)then
	write(2,*)nkp-(m-1)*subnkp
	else
	write(2,*)subnkp
	endif
	write(2,'(A)')"Reciprocal"

        if(m<nfile)then
	do i=1,subnkp
	j=(m-1)*subnkp+i
	write(2,99)kp(j,1:3),wei(i)
        enddo
        else
        do i=1,nkp-(m-1)*subnkp
        j=(m-1)*subnkp+i
        write(2,99)kp(j,1:3),wei(i)
        enddo
        endif

	close(2)
	enddo
99	format(3f20.14,I14)
	close(1)
	close(3)
	end

