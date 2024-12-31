      program bsplt_k_projected

!     plots the projected bands, both xmgrace (2-D graph) and opendx
!     (a colored rubber-sheet). More than one bs_projected.dat file
!     can be combined by simply cat'ing them together.

      implicit none

      integer,       parameter    :: fp = selected_real_kind(13)
      real(kind=fp),parameter    ::  one_pi = 3.141592653589793238462643383279502884197_fp

!        file handles
      integer,       parameter    :: outfh=6

!     file name and handle for band data
      integer,          parameter :: bskfh = 73
      character(len=*), parameter :: bskfn ='bs_projected.dat'

      logical                    :: l_exist
      integer                    :: ik,ikpt,nkpts,kpt_lines,l,nps,ic
      integer                    :: js,jspin,jspins,ios
      integer                    :: nb,nbands,ip
      real(kind=fp)              :: ef,eig,wt,x_kpt,bkpt(3)
      real(kind=fp), allocatable :: x_lines(:)
      real(kind=fp)              :: emin,emax,de,dx,wtmin,dele,y2x,wtcut
      real(kind=fp)              :: emin0,emax0,omega,smearing
      integer                    :: i1,i2,it,nx

      real(kind=fp), allocatable :: e(:),ewts(:,:,:),br(:),etmp(:,:,:)
      real(kind=fp), allocatable :: x(:)

      character(len=132)         :: line,fname1,fname2

      real(fp)                   :: ewin_max,ewin_min,estep 

        real(fp),allocatable :: CELTOT(:,:,:),proj_wt(:,:,:)
        real(fp),allocatable :: x_kpt_t(:)

     
      write(outfh,'(/,"-------- k-projected bands ------",/)')
      inquire( file=bskfn, exist=l_exist )
      if( .not. l_exist ) then
         write(outfh,'("bsplt_k_project: file does not exist")')
         stop "bsplt_k_project: file does not exist"
      endif

      emin =  9999._fp
      emax = -9999._fp

      open(bskfh,file=bskfn,form='formatted')
      read(bskfh,*) nkpts,jspins,ef
      read(bskfh,*) kpt_lines

      if( .not. allocated( x_lines ) ) allocate( x_lines(kpt_lines) )
      read(bskfh,*) x_lines(1:kpt_lines )
      read(bskfh,*) ik,nbands

      allocate(x_kpt_t(nkpts),CELTOT(nbands,nkpts,jspins))
      allocate(proj_wt(nbands,nkpts,jspins))

!--->    read in the information and write to plot file

      rewind(bskfh)
      read(bskfh,*) nkpts,jspins,ef
      read(bskfh,*) kpt_lines
      read(bskfh,*) x_lines(1:kpt_lines )

      do jspin = 1,jspins

         do ikpt = 1,nkpts

            read(bskfh,*) ik,nbands,js,bkpt,x_kpt

            x_kpt_t(ikpt) = x_kpt

            do nb = 1,nbands
               read(bskfh,*) eig,wt

               CELTOT(nb,ikpt,jspin) = eig
               proj_wt(nb,ikpt,jspin)= wt

               emin = min( emin, eig )
               emax = max( emax, eig )
            enddo

         enddo

      enddo  ! spins

!===>    xmgrace plot
      open(30,file='bsplot.agr',                             &
     &                 status='unknown',form='formatted')

      write(30,'("@version 50120")')
      write(30,'("@page size 792, 612")')
      write(30,'("@with line")')
      write(30,'("@    line on")')
      write(30,'("@    line loctype world")')
      write(30,'("@    line g0")')
      write(30,'("@    line  0, 0,",f12.6,", 0")')             &
     &                        x_lines( kpt_lines )
      write(30,'("@    line linewidth 1.0")')
      write(30,'("@    line linestyle 1")')
      write(30,'("@    line color 1")')
      write(30,'("@    line arrow 0")')
      write(30,'("@    line arrow type 0")')
      write(30,'("@    line arrow length 1.0")')
      write(30,'("@    line arrow layout 1.0, 1.0")')
      write(30,'("@line def")')
      write(30,'("@with g0")')
      write(30,'("@    world 0, -20,",f12.6,", 10")')          &
     &                 x_lines(kpt_lines)
      write(30,'("@    view  0.25, 0.20, 1.00, 0.85")')
      write(30,'("@    xaxis  on")')
      write(30,'("@    xaxis  tick off")')
      write(30,'("@    xaxis  ticklabel off")')
      write(30,'("@    yaxis  on")')
      write(30,'("@    yaxis  tick on")')
      write(30,'("@    yaxis  tick major 5")')
      write(30,'("@    yaxis  tick minor ticks 1")')
      write(30,'("@    s0 type xysize")')
      write(30,'("@    s0 symbol 1")')
      write(30,'("@    s0 symbol size 0.5")')
      write(30,'("@    s0 symbol color 11")')
      write(30,'("@    s0 symbol linewidth 1.0")')
      write(30,'("@    s0 line type 0")')
      write(30,'("@    s0 line linestyle 1")')
      write(30,'("@    s0 line linewidth 1.0")')

      if(jspins>1)then
         write(30,'("@    s1 type xysize")')
         write(30,'("@    s1 symbol 1")')
         write(30,'("@    s1 symbol size 0.5")')
         write(30,'("@    s1 symbol color 4")')
         write(30,'("@    s1 symbol linewidth 1.0")')
         write(30,'("@    s1 line type 0")')
         write(30,'("@    s1 line linestyle 1")')
         write(30,'("@    s1 line linewidth 1.0")')
      endif

      write(30,'("@target g0.s0")')
      write(30,'("@type xysize")')

      do ikpt = 1,nkpts
        do nb = 1,nbands
          write(30,'(3f14.8)') x_kpt_t(ikpt),                       &
     & (CELTOT(nb,ikpt,1)-ef),proj_wt(nb,ikpt,1)
        enddo
      enddo
      write(30,'("&")')

      if(jspins>1)then
         write(30,'("@target g0.s1")')
         write(30,'("@type xysize")')

         do ikpt = 1,nkpts
           do nb = 1,nbands
              write(30,'(3f14.8)') x_kpt_t(ikpt),                       &
     & (CELTOT(nb,ikpt,2)-ef),proj_wt(nb,ikpt,2)
           enddo
         enddo
         write(30,'("&")')
      endif 

      write(30,'(A)')"#Print out to"
      write(30,'(A)')"@PRINT TO ""bsplot.eps"""
      write(30,'(A)')"@HARDCOPY DEVICE ""EPS"""
      write(30,'(A)')"@DEVICE ""EPS"" OP ""level2"""

!===>    opendx 2-d plot - uniform, 2-d mesh

!--->    set up an energy mesh (default = 0.05 eV)
      de = 0.05_fp
      emin =  emin -ef
      emax = emax -ef
      wtmin = 0.0_fp ! min. value for a point
      wtcut = 0.0_fp ! min. value for a point
      dele  = de
      y2x   = 1.0_fp ! ratio of plot: default 1 is for square
      emin0 = emin
      emax0 = emax
      emin  = floor( emin0 - 0.2_fp )
      emax  = min( real( ceiling( emax0 ), fp ) , 5.0_fp )
!--->    write out defaults and see whether to change:
      write(outfh,'("OpenDx plot defaults:")')
      write(outfh,'("   emin  = ",f12.6,"  (data: ",f12.6,")")')        &
     &emin,emin0
      write(outfh,'("   emax  = ",f12.6,"  (data: ",f12.6,")")')        & 
     &emax,emax0
      write(outfh,'("spacing in energy:")')
      write(outfh,'("   de    = ",f12.6)') de
      write(outfh,'("minimum weight for any state:")')
      write(outfh,'("   wtmin = ",f12.6)') wtmin
      write(outfh,'("energy broadening of state:")')
      write(outfh,'("   dele  = ",f12.6)') dele
      write(outfh,'("ratio of y to x of plot:")')
      write(outfh,'("   y2x   = ",f12.6)') y2x
      write(outfh,'("   wtcut = ",f12.6)') wtcut
      write(outfh,'("if you wish to change, then name=value")')
      write(outfh,'("to end, ""x"" on separate line")')
      loop_in: do
         write(outfh,'("> ")',advance='no')
         read (5,'(a)') line
         line = adjustl( line )
         do
            if( len_trim(line) == 0 ) exit
            if( line(1:1) == 'x' ) exit loop_in
            ic = index( line, '=' )
            if( ic == 0 ) then
               write(outfh,'("bad input. need name=value")')
               exit
            endif
            if( line(1:4) == 'dele' ) then
               read(line(ic+1:),*) dele
               write(outfh,'("dele  = ",f12.6)') dele
            elseif( line(1:2) == 'de' ) then
               read(line(ic+1:),*) de
               write(outfh,'("de    = ",f12.6)') de
            elseif(line(1:4) == 'emin' ) then
               read(line(ic+1:),*) emin
               write(outfh,'("emin  = ",f12.6)') emin
            elseif(line(1:4) == 'emax' ) then
               read(line(ic+1:),*) emax
               write(outfh,'("emax  = ",f12.6)') emax
            elseif(line(1:5) == 'wtmin' ) then
               read(line(ic+1:),*) wtmin
               write(outfh,'("wtmin = ",f12.6)') wtmin
            elseif(line(1:3) == 'y2x' ) then
               read(line(ic+1:),*) y2x
               write(outfh,'("y2x   = ",f12.6)') y2x
            elseif(line(1:5) == 'wtcut' ) then
               read(line(ic+1:),*) wtcut
               write(outfh,'("wtcut = ",f12.6)') wtcut
            else
               write(outfh,'("unrecognized input: ",a)') trim(line)
            endif
            it = scan( line(ic:) , 'edw' ) ! look for another value
            if( it == 0 ) then
               ic = len_trim(line)
               if( line(ic:ic) == 'x' ) exit loop_in
               exit
            endif
            line = line(ic+it-1: )
         enddo
      enddo loop_in

      ewin_max = emax
      ewin_min = emin
      estep    = 1.0

      nps = ( emax - emin )/de

      allocate(ewts(0:nps,nkpts,jspins),x(nkpts) )

      x(:) = x_kpt_t(:)
!--->    read in the information and write to plot file

      do jspin = 1,jspins
         do ikpt = 1,nkpts

            do ip = 0,nps
               omega = ip*de + emin + ef                

              do nb = 1,nbands
                smearing = &
     & dele/one_pi/((omega-CELTOT(nb,ikpt,jspin))**2+dele**2) 
                ewts(ip,ikpt,jspin) = ewts(ip,ikpt,jspin) + &
     & smearing * proj_wt(nb,ikpt,jspin)
              enddo

            enddo

         enddo

      enddo  ! spins

!        rescale x coordinates
      dx =  y2x*(emax-emin)/x(nkpts)
      x  = x*dx
!--->    write out the opendx files:
      dx = (emax-emin)/(y2x*(nkpts-1))  ! make the plot square as a default

      do js = 1,jspins
         write(fname1,'("plot_bands_",i0,".general")') js
         write(fname2,'("plot_bands_",i0,".dat")') js

         open(15,file=fname1,status='unknown',form='formatted')

         write(15,'("# OpenDx file")')
         write(15,'("# to run, do the following steps:")')
         write(15,'("# 1) --> Run Visual Programs")')
         write(15,'("# 2) --> select Band-plots.net &
     &(best to make a copy)")')
         write(15,'("# 3) --> modify as desired")')
!        write(15,'("# 1) --> Import Data")')
!        write(15,'("# 2) --> Grid or Scattered File")')
!        write(15,'("# 3) --> Describe Data")')
!        write(15,'("# 4) --> Data File (",a,")")') trim( fname1 )
!        write(15,'("# 5) close in File menu")')
!        write(15,'("# 6) --> Visualize Data")')
         write(15,'("#")')
         emax = maxval( ewts(:,:,js) )
         write(15,'("#")')
         write(15,'("#  maximum value:",es16.6)') emax
         write(outfh,'("spin ",i0,"  maximum value:",es16.6)') js,emax
         write(15,'("#")')
         write(15,'("file = ",a)') trim( fname2 )
         write(15,'("grid = ",i0," x ",i0)') nkpts,nps+1
         write(15,'("format = ascii")')
         write(15,'("majority = row")')
         write(15,'("field = field0")')
         write(15,'("structure = scalar")')
         write(15,'("type = float")')
         write(15,'("dependency = positions")')
!        write(15,'("positions = regular, regular, 0, ",f12.8,",",       &
!    &           f12.8,",",f12.8)') dx,emin,de
         write(15,'("positions = irregular, regular")',advance='no') 
         write(15,'(10000(", ",f0.6))',advance='no') x(1:nkpts)
         write(15,'(f0.6,", ",f0.6)') emin,de
!      write(15,'("positions = irregular, regular,",(1x,f0.6))',advance='no') &
!    &    x(1:nkpts),emin,de
!        write(15,'(1x,f0.6,",")',advance='no') x(1:nkpts),emin,de
         write(15,'(/,"end")')
         close(15)

         open(16,file=fname2,status='unknown',form='formatted')

         do i1=0,nps
           do i2 = 1,nkpts
              if(ewts(i1,i2,js)<=wtcut)then
                 ewts(i1,i2,js) = 0.0_fp
              endif
           enddo
         enddo

         write(16,'(4es20.12)') ( ( ewts(i1,i2,js), i1 = 0,nps ),        &
     &                                              i2 = 1,nkpts )
         close(16)

      enddo

      if(jspins==2)then
         ! all spins together
         fname1 = "plot_bands_total.general"
         fname2 = "plot_bands_total.dat"
 
         open(15,file=fname1,status='unknown',form='formatted')

         write(15,'("# OpenDx file")')
         write(15,'("# to run, do the following steps:")')
         write(15,'("# 1) --> Run Visual Programs")')
         write(15,'("# 2) --> select Band-plots.net &
     &(best to make a copy)")')
         write(15,'("# 3) --> modify as desired")')
       ! write(15,'("# 1) --> Import Data")')
       ! write(15,'("# 2) --> Grid or Scattered File")')
       ! write(15,'("# 3) --> Describe Data")')
       ! write(15,'("# 4) --> Data File (",a,")")') trim( fname1 )
       ! write(15,'("# 5) close in File menu")')
       ! write(15,'("# 6) --> Visualize Data")')
         write(15,'("#")')
         emax = maxval( ewts(:,:,:) )
         write(15,'("#")')
         write(15,'("#  maximum value:",es16.6)') emax
         write(outfh,'("spin ",i0,"  maximum value:",es16.6)') 1,emax
         write(15,'("#")')
         write(15,'("file = ",a)') trim( fname2 )
         write(15,'("grid = ",i0," x ",i0)') nkpts,nps+1
         write(15,'("format = ascii")')
         write(15,'("majority = row")')
         write(15,'("field = field0")')
         write(15,'("structure = scalar")')
         write(15,'("type = float")')
         write(15,'("dependency = positions")')
!!       write(15,'("positions = regular, regular, 0, ",f12.8,",",       &
!!   &           f12.8,",",f12.8)') dx,emin,de
       write(15,'("positions = irregular, regular")',advance='no') 
       write(15,'(10000(", ",f0.6))',advance='no') x(1:nkpts)
       write(15,'(f0.6,", ",f0.6)') emin,de
!!     write(15,'("positions = irregular, regular,",(1x,f0.6))',advance='no') &
!!   &    x(1:nkpts),emin,de
!!       write(15,'(1x,f0.6,",")',advance='no') x(1:nkpts),emin,de
         write(15,'(/,"end")')
         close(15)

         open(16,file=fname2,status='unknown',form='formatted')
         write(16,'(4es20.12)')((ewts(i1,i2,1)+ewts(i1,i2,2),i1 =0,nps),&
     &                                              i2 = 1,nkpts )
         close(16)

       endif


     ! for matplotlib and gnuplot
      do js=1,jspins
        write(fname2,'("band_kproj_",i0,".dat")') js
        open(1,file=fname2,status='unknown',form='formatted')
        write(1,'(A)')"# ik    E (eV)      weight"
 
        do ikpt = 1,nkpts
          do nb = 0,nps
            write(1,'(f12.6,2f12.6)')x_kpt_t(ikpt),nb*de+emin,ewts(nb,ikpt,js)
          enddo
          write(1,*)
        enddo
        close(1)
      enddo
      if(jspins>1)then
        open(1,file="band_kproj_tot.dat",status='unknown')
        write(1,'(A)')"# ik    E (eV)      weight"

        do ikpt = 1,nkpts
          do nb = 0,nps
            write(1,'(f12.6,2f12.6)')x_kpt_t(ikpt),nb*de+emin,ewts(nb,ikpt,1)+ewts(nb,ikpt,2)
          enddo
          write(1,*)
        enddo
        close(1)
      endif

        ! prepare the script for gnuplot plotting
        open(66,file="plotbnds.gnu",status='unknown')
        write(66,'("set term postscript enhanced color ")')
        write(66,'("#set term png truecolor transparent ")')
        write(66,'("set output ""plotbnds_up.eps""")')
        write(66,'("#set output ""plotbnds_up.png""")')
        write(66,'("set multiplot")')
        write(66,'("set pm3d map")')
        write(66,'("set palette rgb 21,22,23;")')
        write(66,'("set lmargin 13")')
        write(66,'("set rmargin 0.3")')
        write(66,'("#unset ytics")')
        write(66,'("set tmargin 0")')
        write(66,'("set bmargin 0")')
        write(66,'("#set tics out")')
        write(66,'("unset xtics")')
        write(66,'("unset xlabel")')
        write(66,'("set parametric")')
        write(66,'("const=0")')
        write(66,'("set trange [",f10.5," :",f10.5,"]")')ewin_min,ewin_max
        write(66,'("set xrange [0 : ",f12.6,"]")')x_lines(kpt_lines)
        write(66,'("set yrange [",f10.5," :",f10.5,"]")')ewin_min,ewin_max
        write(66,'("set cbrange [0:",f6.2,"]")')emax*0.5
        write(66,'(" ")')
        write(66,'("# tics pointing out")')
        write(66,'("set tics out")')
        write(66,'(" ")')
        write(66,'("# hide tics on x")')
        write(66,'("unset xtics")')
        write(66,'("set ytics font ""Helvetica,32"" nomirror #&
     & hide tics on the opposite side")')
        write(66,'("#set ytics ",3f8.3," font ""Helvetica,32"" nomirror #&
     & hide tics on the opposite side")')ewin_min,estep,ewin_max
        write(66,'(" ")')
        write(66,'("# set tics width")')
        write(66,'("set border lw 2 ")')
        write(66,'("set ylabel offset -5.5,0 ")')
        write(66,'("set ylabel ""E (eV)"" font ""Helvetica,32"" ")')
        write(66,'("set tmargin 1")')
        write(66,'("set bmargin 1")')
        write(66,'(" ")')

        write(66,'("#set label ""{/Symbol G}""  &
     & font ""Helvetica,32"" at 0.50,-15.3")')

        write(66,'("set label ""K"" font ""Helvetica,32"" at",f10.5,",",f10.5)')&
     & 0.0-x_lines(kpt_lines)/40,ewin_min-(ewin_max-ewin_min)/20

        do nx=1,kpt_lines
          write(66,'("set label ""K"" font ""Helvetica,32"" at",f10.5,",",f10.5)')&
     & x_lines(nx)-x_lines(kpt_lines)/40,ewin_min-(ewin_max-ewin_min)/20
        enddo
        
        write(66,'(" ")')
        write(66,'("splot  ""band_kproj_1.dat"" ")')

        write(66,'(" ")')

!     if(jspins>1)then
!       write(66,'("splot  ""band_kproj_2.dat"" ")')
!     endif
        close(66)


      end program bsplt_k_projected

