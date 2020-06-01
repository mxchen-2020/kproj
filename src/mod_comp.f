      module mod_comp
         !-------------------------
         ! paramters
         ! M. X. Chen -- 10/2010
         !-------------------------
         integer, parameter  :: I4B = selected_int_kind(9) 
         integer, parameter  :: I2B = selected_int_kind(4) 
         integer, parameter  :: I1B = selected_int_kind(2) 

         integer, parameter  :: SP = kind(1.0) 
         integer, parameter  :: DP =SELECTED_real_KIND(10) 

!        INTEGER, PARAMETER  :: fp=SELECTED_REAL_KIND(13)
         INTEGER, PARAMETER  :: fp=SELECTED_REAL_KIND(10)
         INTEGER, PARAMETER  :: QS=SELECTED_REAL_KIND(5)

         integer, parameter  :: SPC = kind( (1.0,1.0) ) 
         integer, parameter  :: DPC = kind( (1.0_DP,1.0_DP) ) 

         integer, parameter  :: LGT = kind( .true. ) 

         real(DP), parameter :: PI = 3.141592653589793238462643383279502884197_DP
         real(DP), parameter :: TWOPI = 2._DP*PI
         real(DP), parameter :: FOURPI = 4._DP*PI

         complex(DPC), parameter :: im = (0._DP,1._DP)
         complex(DPC), parameter :: unitary = (1._DP,0._DP)

         !  Some important Parameters, to convert to a.u.
         !  - AUTOA  = 1. a.u. in Angstroem
         !  - RYTOEV = 1 Ry in Ev
         !  - EVTOJ  = 1 eV in Joule
         !  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
         !  - BOLKEV = Boltzmanns constant in eV/K
         !  - BOLK   = Boltzmanns constant in Joule/K

         real(DP), parameter :: AUTOA=0.529177249_DP,RYTOEV=13.605826_DP
         real(DP), parameter :: EVTOJ=1.60217733E-19_DP,AMTOKG=1.6605402E-27_DP, &
     &           BOLKEV=8.6173857E-5_DP,BOLK=BOLKEV*EVTOJ, HPLANK=6.6262E-34_DP
         real(DP), parameter ::  HARTREE = 27.211652_DP

         real(DP), parameter ::  eps = 1.d-7
         real(DP), parameter ::  eps7  = 1e-7_DP
         real(DP), parameter ::  eps12 = 1e-12_DP

         !  extra large buffer length
         integer, parameter  :: xl_buffer=16384

         ! record time consuming
          real(DP)            :: time_start,time_end,consumetime

         integer,parameter,public  :: IUO=1  ! OUTKPROJ 
         integer,parameter,public  :: IUI=2  ! INKPROJ
         integer,parameter,public  :: IUW=3  ! WAVECAR
         integer,parameter,public  :: IUEG=4 ! eigen values
         
       contains

       subroutine init()
         character(8)         :: date
         character(10)        :: time
         character(5)         :: zone
         integer              :: values(8)

!        IUO  = 1   ! OUTKPROJ
!        IUI  = 2   ! INKPROJ
!        IUW  = 3   ! WAVECAR
!        IUEG = 4   ! energy band

         open(IUO,file="OUTKPROJ",status='unknown')
         write(IUO,'(3X,A)')"******** PROGRAM KPROJ ************"
         write(IUO,'(3X,A)')"****** MINGXING CHEN (HUNNU)*******"

         !
         ! running time
         !
         call date_and_time(date,time,zone,values)
         write(IUO,'(5X,A,I2,A,I2,A,I4,2X,I2,A,I2)')"Running Time: ",&
     & values(2),"/",values(3),"/",values(1),values(5),":",values(6)

       end subroutine init

       subroutine fstop(stout, errormessage )
       implicit none
 
       integer, intent(in)          :: stout
       character(len=*), intent(in) :: errormessage

        write(stout,'(5X,A)')&
      & "-----------------------------------------------"
       write(stout,'(5X,A)') errormessage
       write(stout,'(5X,A)')&
      & "-----------------------------------------------"
       stop

       return
       end subroutine fstop

       subroutine job_completed(consumetime)

       real(DP)   :: consumetime

       write(IUO,'(/,3X,A)')&
     & "=============== Ending KPROJ ================="
       write(IUO,'(/,3X,A,f12.3)') &
     & "Total CPU time used (sec): ", consumetime 

       end subroutine job_completed
    

       function job_time()
       !-------------------------------------
       ! calculate elapsed CPU time (seconds)
       !-------------------------------------
       implicit none

       real(DP) :: job_time

       ! t0 contains the time of the first call
       ! t1 contains the present time
       real(DP) :: time0, time1
       logical :: first=.true.
       save first, time0
       !
       call cpu_time(time1)
       !
       if (first) then
         time0 = time1
         job_time = 0.0_DP
         first = .false.
       else
         job_time = time1 - time0
       endif
       return
      end function job_time

      end module mod_comp

