Module read_information_from_qe
     USE qes_types_module
     USE Fox_dom
     USE mod_comp 
     
     contains

     SUBROUTINE read_spin(lsda,noncolin,spinorbit)
     IMPLICIT NONE
     TYPE(Node), POINTER :: xml_node,mydoc,out_node
     TYPE(spin_type)     :: obj
     INTEGER             :: ierr
     TYPE(Node), POINTER :: child_node
     TYPE(NodeList), POINTER :: childList
     INTEGER :: childList_size, index, iostati
     LOGICAL,INTENT(out) :: lsda,noncolin,spinorbit
     myDoc => parseFile("data-file-schema.xml")  !分析文件名
     xml_node => getDocumentElement(myDoc)  ! 元素
     obj%tagname = getTagName(xml_node)
     !
     childList => getElementsByTagname(xml_node,"output")
     out_node => item(childList, 0)
     !
     childList => getElementsByTagname(out_node, "lsda")
     child_node => item(childList, 0)
     IF (ASSOCIATED(child_node))&
       CALL extractDataContent(child_node, obj%lsda )
     !
      childList => getElementsByTagname(out_node, "noncolin")
      child_node => item(childlist, 0)
      IF (ASSOCIATED(child_node))&
        CALL extractDataContent(child_node, obj%noncolin )
     !
      childList => getElementsByTagname(out_node, "spinorbit")
      child_node => item(childList, 0)
      IF (ASSOCIATED(child_node))&
        CALL extractDataContent(child_node, obj%spinorbit )
     ! obj%lwrite = .TRUE.
     lsda=obj%lsda
     noncolin=obj%noncolin
     spinorbit=obj%spinorbit 
     !print*,obj%lsda,obj%noncolin,obj%spinorbit
     END SUBROUTINE read_spin

    SUBROUTINE read_cell(CELL) 
    IMPLICIT NONE
    !
    TYPE(Node), POINTER          :: xml_node,myDoc,out_node
    TYPE(cell_type)              :: obj
    INTEGER                      :: ierr
    REAL(DP),INTENT(OUT) ::CELL(3,3)
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    myDoc => parseFile("data-file-schema.xml")  
    xml_node => getDocumentElement(myDoc)  !
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node,"output")
    out_node => item(tmp_node_list, 0)
    !
    tmp_node_list => getElementsByTagname(out_node, "a1")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a1, IOSTAT = iostat_ )
    !
    tmp_node_list => getElementsByTagname(out_node, "a2")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a2, IOSTAT = iostat_ )
    tmp_node_list => getElementsByTagname(out_node, "a3")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a3, IOSTAT = iostat_ )
   ! print*,obj%a1,obj%a2,obj%a3
    CELL(1:3,1) = obj%a1
    CELL(1:3,2) = obj%a2
    CELL(1:3,3) = obj%a3
    END SUBROUTINE read_cell

    SUBROUTINE read_reciprocal_lattice(b)
    IMPLICIT NONE
    !
    TYPE(Node), POINTER          :: xml_node,myDoc,out_node
    TYPE(reciprocal_lattice_type)              :: obj
    INTEGER                      :: ierr
    REAL(DP),INTENT(OUT)::b(3,3)
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    myDoc => parseFile("data-file-schema.xml")  
    xml_node => getDocumentElement(myDoc)  
    obj%tagname = getTagName(xml_node)
    !
    tmp_node_list => getElementsByTagname(xml_node,"output")
    out_node => item(tmp_node_list, 0)
    !
    tmp_node_list => getElementsByTagname(out_node, "b1")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b1, IOSTAT = iostat_ )
    !
    tmp_node_list => getElementsByTagname(out_node, "b2")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b2, IOSTAT = iostat_ )
    tmp_node_list => getElementsByTagname(out_node, "b3")
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b3, IOSTAT = iostat_ )
    b(1:3,1) = obj%b1
    b(1:3,2) = obj%b2
    b(1:3,3) = obj%b3
   ! print*,b
    END SUBROUTINE

    SUBROUTINE read_nbnd(NBANDS)
    IMPLICIT NONE
    !
    TYPE(Node), POINTER    :: xml_node,Mydoc,out_node
    TYPE(band_structure_type)       :: obj
    INTEGER                :: ierr,NBANDS
    !
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    myDoc => parseFile("data-file-schema.xml")  
    xml_node => getDocumentElement(myDoc)  
    !
    !
    tmp_node_list => getElementsByTagname(xml_node,"output")
    tmp_node_list_size = getLength(tmp_node_list)
    out_node => item(tmp_node_list, 0)
    !
    !
    tmp_node_list => getElementsByTagname(out_node, "nbnd")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size>0) THEN
    obj%nbnd_ispresent = .TRUE.
    tmp_node => item(tmp_node_list, 0)
    CALL extractDataContent(tmp_node, obj%nbnd , IOSTAT = iostat_)
    NBANDS=obj%nbnd
    END IF
    
    tmp_node_list => getElementsByTagname(out_node, "nbnd_up")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size>0) THEN
    obj%nbnd_up_ispresent = .TRUE.
    tmp_node => item(tmp_node_list, 0)
    CALL extractDataContent(tmp_node, obj%nbnd_up , IOSTAT = iostat_)
    NBANDS=obj%nbnd_up
    END IF
    
    tmp_node_list => getElementsByTagname(out_node, "nbnd_dw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size>0) THEN
    obj%nbnd_dw_ispresent = .TRUE.
    tmp_node => item(tmp_node_list, 0)
    CALL extractDataContent(tmp_node, obj%nbnd_dw , IOSTAT = iostat_)
    END IF
    END SUBROUTINE read_nbnd

     SUBROUTINE read_nkp(NKPTS)
     TYPE(Node) ,POINTER      :: xml_node,myDoc,out_node
     TYPE(k_points_IBZ_type)  :: obj
     INTEGER                  :: ierr
     INTEGER,INTENT(OUT)      :: NKPTS
     !
     TYPE(Node), POINTER :: tmp_node
     TYPE(NodeList), POINTER :: tmp_node_list
     INTEGER :: tmp_node_list_size, index, iostat_
 
     myDoc => parseFile("data-file-schema.xml")  
     xml_node => getDocumentElement(myDoc)  
     obj%tagname = getTagName(xml_node)
     !
     tmp_node_list => getElementsByTagname(xml_node,"output")
     out_node => item(tmp_node_list, 0)
     !
     tmp_node_list => getElementsByTagname(out_node, "nk")
     tmp_node => item(tmp_node_list, 0)
     CALL extractDataContent(tmp_node, obj%nk , IOSTAT = iostat_)
     NKPTS = obj%nk
    ! print*,NKPTS
     END SUBROUTINE read_nkp

     SUBROUTINE read_cutoff(ENMAX)
     IMPLICIT NONE
     !
     TYPE(Node),  POINTER                 :: xml_node,myDoc,out_node
     TYPE(basis_type)                     :: obj
     INTEGER                              :: ierr
     REAL(DP),INTENT(OUT)                 :: ENMAX
     !
     TYPE(Node), POINTER :: tmp_node
     TYPE(NodeList), POINTER :: tmp_node_list
     INTEGER :: tmp_node_list_size, index, iostat_
     !
     myDoc => parseFile("data-file-schema.xml")  
     xml_node => getDocumentElement(myDoc)  
     obj%tagname = getTagName(xml_node)
     !
     tmp_node_list => getElementsByTagname(xml_node,"output")
     out_node => item(tmp_node_list, 0)
     !
     tmp_node_list => getElementsByTagname(out_node, "ecutwfc")
     tmp_node => item(tmp_node_list, 0)
     IF (ASSOCIATED(tmp_node))&
     CALL extractDataContent(tmp_node, obj%ecutwfc, IOSTAT = iostat_ )
     !print*,obj%ecutwfc
     ENMAX = 2*obj%ecutwfc*13.605826_DP
     END SUBROUTINE read_cutoff
     
     SUBROUTINE read_kpoint_energy(KPTVEC,eigenvalues,weight_)
     IMPLICIT NONE
     TYPE(Node),  POINTER                 :: xml_node,myDoc
     TYPE(k_point_type)                   :: obj
     TYPE(vector_type)                    :: obj2
     INTEGER                              :: ierr,i,j
     real(DP)  ::b(3,3),b_inverse(3,3)
     TYPE(Node), POINTER :: tmp_node,tmp_node1
     TYPE(NodeList), POINTER :: tmp_node_list,tmp_node_list1
     INTEGER :: tmp_node_list_size, index, iostat_,tmp_node_list_size1
     REAL(DP),allocatable,INTENT(OUT)  :: KPTVEC(:,:)
     REAL(DP),allocatable,INTENT(OUT)  :: eigenvalues(:,:),weight_(:)
     myDoc => parseFile("data-file-schema.xml")  
     xml_node => getDocumentElement(myDoc)  
     !进入根目录节点
 
     tmp_node_list => getElementsByTagname(xml_node,"output")
     tmp_node => item(tmp_node_list, 0)
     !进入根目录节点下的output节点
     call read_reciprocal_lattice(b)
     call inverse(b,b_inverse)
 
     tmp_node_list1 => getElementsByTagname(tmp_node, "ks_energies")
     tmp_node_list_size1 = getLength(tmp_node_list1)
     allocate(KPTVEC(3,tmp_node_list_size1))
     allocate(weight_(tmp_node_list_size1))
     do index = 1,tmp_node_list_size1
     tmp_node1 => item(tmp_node_list1, index-1)
 !------------------------------------------------------------------!进入k_point节点
     tmp_node_list => getElementsByTagname(tmp_node1, "k_point")
     tmp_node_list_size = getLength(tmp_node_list)
     tmp_node => item( tmp_node_list, 0 )
     obj%tagname = getTagName(tmp_node)
      
     IF (hasAttribute(tmp_node, "weight")) &
     & CALL extractDataAttribute(tmp_node, "weight", obj%weight)
     weight_(index)=obj%weight
     IF (hasAttribute(tmp_node, "label")) &
     & CALL extractDataAttribute(tmp_node, "label", obj%label)
       CALL extractDataContent(tmp_node, obj%k_point )
     KPTVEC(1:3,index)=matmul(obj%k_point(1:3),b_inverse)
     END DO
     DO i=1,3
        DO j=1,tmp_node_list_size1
         if (abs(KPTVEC(i,j))<1E-10) KPTVEC(i,j)=0
        ENDDO
     ENDDO
!------------------------------------------------------------------!进入eigenvalues节点
     tmp_node_list => getElementsByTagname(tmp_node1, "eigenvalues")
   
     tmp_node_list_size = getLength(tmp_node_list)
     tmp_node => item( tmp_node_list, 0 )
     obj2%tagname = getTagName(tmp_node)
       
     IF (hasAttribute(tmp_node, "size")) &
     & CALL extractDataAttribute(tmp_node, "size", obj2%size)
     allocate(eigenvalues(obj2%size,tmp_node_list_size1))

     do index = 1,tmp_node_list_size1
     tmp_node1 => item(tmp_node_list1, index-1)

     tmp_node_list => getElementsByTagname(tmp_node1, "eigenvalues")
     tmp_node => item( tmp_node_list, 0 )

     allocate(obj2%vector(obj2%size))
     CALL extractDataContent(tmp_node, obj2%vector) !get eigenvalues
     eigenvalues(:,index)=obj2%vector*27.211652_DP
     deallocate(obj2%vector)
     ENDDO
!---------------------------------------------------------------------!    
     END SUBROUTINE read_kpoint_energy
     
     SUBROUTINE inverse(b_,b_inverse)
     real(DP) :: omtil,b1_(3),b2_(3),b3_(3)
     real(DP),intent(out)::b_inverse(3,3)
     real(DP) ::b_(3,3)
     integer :: i,j
     b1_=b_(1:3,1)
     b2_=b_(1:3,2)
     b3_=b_(1:3,3)
     det_b_ = b1_(1)*b2_(2)*b3_(3) + b2_(1)*b3_(2)*b1_(3) +          &
     &         b3_(1)*b1_(2)*b2_(3) - b1_(3)*b2_(2)*b3_(1) -         &
     &         b2_(3)*b3_(2)*b1_(1) - b3_(3)*b1_(2)*b2_(1)
     ! b_=transpose(b_)
     b_inverse(1,1)=  b_(2,2)*b_(3,3)-b_(2,3)*b_(3,2)
     b_inverse(1,2)=-(b_(2,1)*b_(3,3)-b_(2,3)*b_(3,1))
     b_inverse(1,3)=  b_(2,1)*b_(3,2)-b_(2,2)*b_(3,1)
     b_inverse(2,1)=-(b_(1,2)*b_(3,3)-b_(1,3)*b_(3,2))
     b_inverse(2,2)=  b_(1,1)*b_(3,3)-b_(1,3)*b_(3,1)
     b_inverse(2,3)=-(b_(1,1)*b_(3,2)-b_(1,2)*b_(3,1))
     b_inverse(3,1)=  b_(1,2)*b_(2,3)-b_(1,3)*b_(2,2)
     b_inverse(3,2)=-(b_(1,1)*b_(2,3)-b_(1,3)*b_(2,1))
     b_inverse(3,3)=  b_(1,1)*b_(2,2)-b_(2,1)*b_(1,2)
     b_inverse=b_inverse/det_b_
     do i=1,3
       do j=1,3
       if(abs(b_inverse(i,j))<1E-7) b_inverse(i,j)=0
       end do
     end do
     END SUBROUTINE inverse

 END MODULE
