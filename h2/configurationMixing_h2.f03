 program coeffsMixing

    !   This program computes configuration mixing between various orbitals in a
    !   series of related geometries

    !****x* Main/matDiffNOCalc
    !   A. M. Kinyua, 2023
    !*  NAME
    !*      Natural Orbitals Calculator 
        
        use mqc_gaussian
        use iso_fortran_env, only: int32, int64, real64

    !
    ! Variables ----------------------------------------------------
        implicit none 
        type(mqc_gaussian_unformatted_matrix_file)::fileInfo
        character(len=:),allocatable::command,fileName,help_path
        character(len=256),dimension(:),allocatable::fileList
        integer(kind=int64)::iOut=6,iPrint=1,numFile=0,i,j,k,l,stat_num,iUnit
        type(mqc_scf_integral)::overlapIntegral,mo_listIntegral_1,mo_listIntegral_i,overlapQuantifier
        type(mqc_scf_integral)::homo,homo_1,homo_2,homo_3,homo_4,homo_5,homo_6,homo_7
        type(mqc_gaussian_unformatted_matrix_file)::temp_file
        type(mqc_scalar)::normHomo,normHomo_1,normHomo_2,c2c0,c2c1,c1c0
        !integer(kind=int64)::mo1=48,mo2=49,mo3=50,mo4=51,mo5=52,mo6=53,mo7=54,mo8=55 
        integer(kind=int64)::mo1=1,mo2=2,mo3=3,mo4=4
        type(mqc_matrix)::tmpmatrix
        type(mqc_vector)::tmpvector
    !
    !     Parse input Options ---------------------------------------
    !   

        j=1
        do i=1,command_argument_count()
            if(i.ne.j) cycle
            call mqc_get_command_argument(i,command)
           if(command.eq.'-f') then
    !
    !*      -f text file                   Input text file with matrix files.
    !*                                     The first line contains the number of matrix files in the
    !*                                     input, and then on each line is a separate matrix file.
    !*
               call mqc_get_command_argument(i+1,fileName)
               j = i+2
            endIf
       endDo


    !
    !*  Parse input file ---------------------------------------------
    !
        if(.not.allocated(fileName)) call mqc_error('No input file provided', iOut)
        open(newunit=iUnit, file=fileName,status='old',iostat=stat_num)
        if(stat_num/=0) call mqc_error_a('Error opening file',iOut,'fileName',fileName)
        read(unit=iUnit,fmt='(i20)',iostat=stat_num) numFile
        if(stat_num/=0) call mqc_error('Error reading file number',iOut)
        allocate(fileList(numFile))
        do i = 1, numFile
            read(unit=iUnit,fmt='(A)',iostat=stat_num) fileList(i)
            if((stat_num<0).and.(i<=numFile)) call mqc_error('File EOF reached early',iOut)
        endDo
        close(unit=iUnit)

    
    ! getting overlap integral and printing it 
        call temp_file%getESTObj('overlap',est_integral=overlapIntegral,fileName=fileList(1))
        call overlapIntegral%print(6,'Overlap Integral')       
        call temp_file%getESTObj('mo coefficients',est_integral=mo_listIntegral_1,filename=fileList(1))
        call mo_listIntegral_1%print(6,'mo coeffs in absence of pertubation')
    
    ! loop over all matrix files and get mo coeffs, especially alpha mo coeffs
        do i=1,numFile
            call temp_file%getESTObj('mo coefficients',est_integral=mo_listIntegral_1,filename=fileList(1))
            call temp_file%getESTObj('mo coefficients',est_integral=mo_listIntegral_i,filename=fileList(i))
           
    ! <C_o|S|C_t>
            write(*,*)'    '
            overlapQuantifier = matmul(dagger(mo_listIntegral_1),matmul(overlapIntegral,mo_listIntegral_i))
            call overlapQuantifier%print(6,'<C_1|S|C_'//trim(num2char(i))//'>')
    
    ! MO1 ----> this is homo for h2        
            homo=overlapQuantifier%orbitals('useStrings',betaOrbsIn=[mo1],axis=2)
            tmpmatrix = homo%getBlock('beta')
            tmpvector = tmpmatrix%vat([0],[1])
            call tmpvector%print(6,'HOMO vector at '//trim(num2char(i)))
            call mqc_print(tmpvector%norm('F'),6,'HOMO vector norm at '//trim(num2char(i)))
    
    ! MO2 ------> this is lumo for h2
            homo_1=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo2],betaOrbsIn=[mo2],axis=2)
            tmpmatrix = homo_1%getBlock('beta')
            tmpvector = tmpmatrix%vat([0],[1])
            call tmpvector%print(6,'LUMO vector at '//trim(num2char(i)))
            call mqc_print(tmpvector%norm('F'),6,'LUMO vector norm at '//trim(num2char(i)))
    
    ! MO3 -------> this is lumo+1 for h2
            homo_2=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo3],betaOrbsIn=[mo3],axis=2)
            tmpmatrix = homo_2%getBlock('beta')
            tmpvector = tmpmatrix%vat([0],[1])
            call tmpvector%print(6,'LUMO+1 vector at '//trim(num2char(i)))
            call mqc_print(tmpvector%norm('F'),6,'LUMO+1 vector norm at '//trim(num2char(i)))
    
    ! MO4 --------> this is lumo+2 for h2
            homo_3=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo4],betaOrbsIn=[mo4],axis=2)
            tmpmatrix = homo_3%getBlock('beta')
            tmpvector = tmpmatrix%vat([0],[1])
            call tmpvector%print(6,'LUMO+2 vector at '//trim(num2char(i)))
            call mqc_print(tmpvector%norm('F'),6,'LUMO+2 vector norm at '//trim(num2char(i)))
    
    ! MO5 --------> this is lumo+2 for h2
           ! homo_4=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo5],betaOrbsIn=[mo5],axis=2)
            !tmpmatrix = homo_4%getBlock('beta')
            !tmpvector = tmpmatrix%vat([0],[1])
            !call tmpvector%print(6,'LUMO vector at '//trim(num2char(i)))
            !call mqc_print(tmpvector%norm('F'),6,'LUMO vector norm at '//trim(num2char(i)))
    
    ! MO6 --------> this is lumo+2 for h2
            !homo_5=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo6],betaOrbsIn=[mo6],axis=2)
            !tmpmatrix = homo_5%getBlock('beta')
            !tmpvector = tmpmatrix%vat([0],[1])
            !call tmpvector%print(6,'LUMO+1 vector at '//trim(num2char(i)))
            !call mqc_print(tmpvector%norm('F'),6,'LUMO+1 vector norm at '//trim(num2char(i)))
    
    ! MO7 --------> this is lumo+2 for h2
            !homo_6=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo7],betaOrbsIn=[mo7],axis=2)
            !tmpmatrix = homo_6%getBlock('beta')
            !tmpvector = tmpmatrix%vat([0],[1])
            !call tmpvector%print(6,'LUMO+2 vector at '//trim(num2char(i)))
            !call mqc_print(tmpvector%norm('F'),6,'LUMO+2 vector norm at '//trim(num2char(i)))
    
    ! MO8 --------> this is lumo+2 for h2
            !homo_7=overlapQuantifier%orbitals('useStrings',alphaOrbsIn=[mo8],betaOrbsIn=[mo8],axis=2)
            !tmpmatrix = homo_7%getBlock('beta')
            !tmpvector = tmpmatrix%vat([0],[1])
            !call tmpvector%print(6,'LUMO+3 vector at '//trim(num2char(i)))
            !call mqc_print(tmpvector%norm('F'),6,'LUMO+3 vector norm at '//trim(num2char(i)))
    
    ! Configuration mixing
            !c2c0 = normHomo_2/normHomo
            !c2c1 = normHomo_2/normHomo_1
            !c1c0 = normHomo_1/normHomo
            !call c2c0%print(6,'|C2|^2/|C0|^2 at '//trim(num2char(i)))
            !call c2c1%print(6,'|C2|^2/|C1|^2 at '//trim(num2char(i)))
            !call c1c0%print(6,'|C1|^2/|C0|^2 at '//trim(num2char(i)))
       endDo
!


contains


 function mqc_integral_normalizer(integral,methodIn) result(norm)
!
       implicit none
       class(mqc_scf_integral),intent(in)::integral
       character(len=1),optional,intent(in)::methodIn
       character(len=1)::method
       type(mqc_scalar)::norm
       type(mqc_matrix)::ret_block
!
       if(Present(methodIn)) then
         method = methodIn
       else
         method = 'F'
       endIf
!      
      ret_block = integral%getBlock('full')
      call ret_block%print(6,'ret block')
      norm = mqc_matrix_norm(ret_block,method)
!
  end function mqc_integral_normalizer


end program coeffsMixing

