    MODULE read_information_from_ab
    use mod_comp
    use mod_input
    Real(DP),public,allocatable :: kpt_ab(:,:)
    
       contains

    SUBROUTINE read_information_from_abinit(ISPIN_,NBANDS_,NKPTS_,CELL_,ENMAX_,NSPINOR_)
    ! https://docs.abinit.org/guide/abinit/#5-the-different-output-files
    CHARACTER(len=8) :: codvsn1 !for 9.x.x---? unknow
    CHARACTER(len=6) :: codvsn2 !for 8.x.x
    INTEGER  :: headform,fform
    INTEGER  :: bantot,cplex,date,icoulomb,intxc,ixc,kptopt,mband,natom,ngfft(3),&
    & nkpt,npsp,nshiftk_orig,nshiftk,nspden,nspinor,nsppol,nsym,ntypat,occopt,&
    & pawcpxocc,pertcase,usepaw,usewvl
    REAL(DP) :: acell(3),cellcharge,ecut,ecutdg,ecutsm,ecut_eff,etotal,&
    & fermie,nelect,qptn(3),residm,rprimd(3,3),shiftk(3),shiftk_orig(3),&
    & stmbias,tphysel,tsmear
    INTEGER,allocatable  :: istwfk(:),nband(:),&
    & npwarr(:),so_psp(:),symafm(:),symrel(:,:,:),typat(:)
    INTEGER  :: kptrlatt(3,3),kptrlatt_orig(3,3)
    REAL(DP),allocatable :: amu(:),occ(:),tnons(:,:),&
    & znucltypat(:),wtk(:),xred(:,:)
    CHARACTER(len=132):: title
    CHARACTER(len=32) :: md5_pseudos
    REAL(DP) :: znuclpsp,zionpsp
    INTEGER  :: pspso,pspdat,pspcod,pspxc,lmn_size
    integer  :: npw,nband_
    !  k-projected information
    INTEGER,INTENT(OUT)  ::ISPIN_,NBANDS_,NKPTS_
    INTEGER,INTENT(OUT),OPTIONAL  ::NSPINOR_
    REAL(DP),INTENT(OUT) ::CELL_(3,3),ENMAX_
      OPEN(unit=IUW,file=trim(adjustl(WFCAB)),form='unformatted',status='unknown')
      ! 
      READ(unit=IUW,iostat=ios) codvsn1,headform,fform
      if(ios/=0) then
        close(unit=IUW)
        OPEN(unit=IUW,file=trim(adjustl(WFCAB)),form='unformatted',status='unknown')
        read(unit=IUW,IOSTAT= ios) codvsn2,headform,fform
      end if
      !
      READ(unit=IUW) bantot,date,intxc,ixc,natom,ngfft(1:3),&
      & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase,usepaw,&
      & ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,&
      & tphysel,tsmear,usewvl,nshiftk_orig,nshiftk,mband
      !
      ALLOCATE(istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt))
      ALLOCATE(so_psp(npsp),symafm(nsym),symrel(3,3,nsym))
      ALLOCATE(typat(natom),kpt_ab(3,nkpt),occ(bantot),tnons(3,nsym))
      ALLOCATE(znucltypat(ntypat),wtk(nkpt))
      ALLOCATE(xred(3,natom),amu(ntypat))
      ! 
      READ(unit=IUW) istwfk(1:nkpt),nband(1:nkpt*nsppol),&
      & npwarr(1:nkpt),so_psp(1:npsp),symafm(1:nsym),symrel(1:3,1:3,1:nsym),&
      & typat(1:natom),kpt_ab(1:3,1:nkpt),occ(1:bantot),tnons(1:3,1:nsym),&
      & znucltypat(1:ntypat),wtk(1:nkpt)
      ! 
      READ(unit=IUW) residm,xred(1:3,1:natom),etotal,fermie,amu(1:ntypat)
      !
      READ(unit=IUW) kptopt,pawcpxocc,nelect,cellcharge,icoulomb,&
      & kptrlatt(3,3),kptrlatt_orig(3,3),shiftk_orig(3),shiftk(3)
      ! 
      do ipsp=1,npsp
      ! (npsp lines, 1 for each pseudo; npsp=ntypat, except if alchemical pseudo-atoms)
      READ(unit=IUW) title,znuclpsp,zionpsp,pspso,&
      &  pspdat,pspcod,pspxc,lmn_size,md5_pseudos
      enddo
      !(in case of usepaw==1, there are some additional records)
      if (usepaw==1)then ! it is useless to kproj
        read(unit=IUW)
        read(unit=IUW)
      endif
      !
      ISPIN_=nsppol
      ENMAX_=ecut*27.211652_DP !Hartree
      CELL_=rprimd !.a.u
      NBANDS_=nband(1)
      NKPTS_=nkpt
      NSPINOR_=nspinor

     END SUBROUTINE read_information_from_abinit
     end module
