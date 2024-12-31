      program proj2mesh

!     put projected bands and projected weights on to a 2D mesh 
!     This will be used by matplotlib for plotting 
!     The y axe, band energy, will be in a uniform mesh

      implicit none

      integer,       parameter    :: fp = selected_real_kind(13)

!        file handles
      integer,       parameter    :: outfh=6, pltfh(2)= (/ 21,22 /)

!     file name and handle for band data
      integer,          parameter :: bskfh = 73
      character(len=*), parameter :: bskfn ='bs_projected.dat'

      logical                    :: l_exist
      integer                    :: ik,ikpt,nkpts,kpt_lines,l,nps,ic
      integer                    :: js,jspin,jspins,ios
      integer                    :: nb,nbands
      real(kind=fp)              :: ef,eig,wt,bkpt(3)
      real(kind=fp), allocatable :: x_lines(:),x_kpt(:)
      real(kind=fp)              :: emin,emax,de,dx,wtmin,dele,y2x
      integer                    :: i1,i2,it,nx

      real(kind=fp), allocatable :: e(:),ewts(:,:,:),br(:)

      character(len=132)         :: line,fname1,fname2
      
      inquire(file=bskfn, exist=l_exist )
      if( .not. l_exist ) then
         write(outfh,'("bsplt_k_project: file does not exist")')
         stop "bsplt_k_project: file does not exist"
      endif

      ! read bs_projected.dat
      open(bskfh,file=bskfn,form='formatted')

      read(bskfh,*) nkpts,jspins,ef
      read(bskfh,*) kpt_lines
      if( .not. allocated( x_lines ) ) allocate( x_lines(kpt_lines) )
      read(bskfh,*) x_lines(1:kpt_lines )
      allocate(x_kpt(nkpts))

      end 

