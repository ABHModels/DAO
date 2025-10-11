module redistribution
    use ieee_arithmetic
    implicit none

    integer,parameter::file_ntemp=70
    integer::file_nener
    real(8),allocatable::file_temper(:),file_kernel(:,:,:),file_ener(:),file_skn(:,:)

    integer::ioskn,iokernel,iotheta,ioener
    integer::debug_unit

    contains

    subroutine  read_compton(lun11,debug,dirname,model_nener)
        character(256) ftemp,fkernel,fener,fskn
        character(256)::out_file_name,temf
        character(80),intent(in)::dirname
        integer,intent(in)::debug
        integer,intent(in)::lun11
        integer,intent(in)::model_nener
        real(8) test_temper(3),test_theta,test_ener(3)
        integer::index_temper,index_ener
        integer m,i,j,k,io_status
        real(8) tstart,tend
        real(8),parameter::mec2 = 511.d0,boltz = 8.617333262145d-8
        real(8),allocatable::tk_slice(:,:)

        file_nener = model_nener
        allocate(file_temper(file_ntemp),&
        &   file_ener(file_nener),file_skn(file_nener,file_ntemp),tk_slice(file_nener,file_nener))

        call cpu_time(tstart)
        ftemp = trim(dirname)//'/theta.dat'
        call test_exist(ftemp)
        fkernel = trim(dirname)//'/kernel.dat'
        call test_exist(fkernel)
        fener = trim(dirname)//'/ener.dat'
        call test_exist(fener)
        fskn = trim(dirname)//'/skn.dat'
        call test_exist(fskn)

        call read_file(ftemp,fkernel,fener,fskn)

        call cpu_time(tend)


        write(lun11,*) 'Compton redistribution file read successfully load'
        write(lun11,*) '------Temperature file: ',trim(ftemp)
        write(lun11,*) '------Kernel file: ',trim(fkernel)
        write(lun11,*) '------Energy file: ',trim(fener)
        write(lun11,*) '------Cross section file: ',trim(fskn)
        write(lun11,*) '------Number of temperature points: ',file_ntemp
        write(lun11,*) '------Number of energy points: ',file_nener
        write(lun11,*) '------Time taken to read compton redistribution file: ',tend-tstart
        flush(lun11)

        ! WRITE TO TEST 
        if (debug.eq.1) then
            test_temper = (/1e7,5e7,5e8/)
            test_ener = (/1.e3,100.e3,500.e3/)
            do m =1,3
            !/// LOOP OVER TEMPERATURES

            ! Get Kernel function at specific temperature
            test_theta = test_temper(m)*boltz/mec2
            call get_index(test_theta,file_ntemp,file_temper,index_temper)
            read(iokernel,rec=index_temper) tk_slice


            do k=1,3
            !/// LOOP OVER ENERGIES

            call get_index(test_ener(k),file_nener,file_ener,index_ener)
            
            write(out_file_name,'(a,F0.3,a,F0.3,a,a)') 'debug/k_',test_ener(k)/1.e3,'keV_',&
            &   log10(test_temper(m)),trim(dirname),'.dat'
            open(newunit=debug_unit,file=trim(out_file_name),status='unknown',action='write',iostat=io_status)
                if (io_status/=0) then
                    write(*,*) 'Error opening kernel file in debug'
                    stop
                endif
            write(debug_unit,8) (file_ener(j),tk_slice(j,index_ener),j=1,file_nener)
            close(debug_unit)

            !/// END LOOP OVER ENERGIES
            enddo

            !/// END LOOP OVER TEMPERATURES
            enddo
        endif
        8 format(2E16.8)
        
        deallocate(tk_slice)
    end subroutine read_compton

    subroutine test_exist(filename)
        character(256),intent(in)::filename
        logical::exist

        inquire(file=filename,exist=exist)
        if(.not.exist) then
            write(*,*) 'Error: file ',trim(filename),' does not exist'
            stop
        end if
    end subroutine test_exist

    subroutine read_file(ftemp,fkernel,fener,fskn)
        character(256),intent(in)::ftemp,fkernel,fener,fskn
        integer::io_status,unit
        integer::rec_leng

        open(newunit=unit,file=ftemp,status='old',action='read',iostat=io_status,form='unformatted')
        if (io_status/=0) then
            write(*,*) 'Error opening temperature file'
            stop
        endif
        read(unit) file_temper
        close(unit)

        open(newunit=unit,file=fener,status='old',action='read',iostat=io_status,form='unformatted')
        if (io_status/=0) then
            write(*,*) 'Error opening energy file'
            stop
        endif
        read(unit) file_ener
        close(unit)

        open(newunit=unit,file=fskn,status='old',action='read',iostat=io_status,form='unformatted')
        if (io_status/=0) then
            write(*,*) 'Error opening skn file'
            stop
        endif
        read(unit) file_skn
        close(unit)

        rec_leng = 8*file_nener*file_nener
        open(unit=iokernel,file=fkernel,status='old',action='read',iostat=io_status,form='unformatted',&
        &   access='direct', recl=rec_leng)
        if (io_status/=0) then
            write(*,*) 'Error opening kernel file'
            stop
        endif

    end subroutine read_file

    subroutine get_cross_section(mo_nener,mo_ener,theta,crosssection)
        ! This routine read the exactly compton cross section from the stored array
        implicit none
        integer,intent(in)::mo_nener
        real(8),intent(in)::mo_ener(mo_nener),theta
        real(8),intent(out)::crosssection(mo_nener)

        integer::index_temper,i
        real(8)::file_cs(file_nener)

        crosssection = 0.0
        call get_index(theta,file_ntemp,file_temper,index_temper)

        do i=1,file_nener
            file_cs(i) = file_skn(i,index_temper)
        enddo

        call interpolate_withex(crosssection,mo_nener,mo_ener,file_cs,file_nener,file_ener)
        
        if(any(crosssection.EQ.0)) then
            print*,'zeros after interpolate'
            stop
        endif

        return
    end subroutine 
    
    subroutine get_ener_pointer(mo_nener,mo_ener,mo_ener_pointer)
        implicit none
        ! This routine get the relationship bewtween model array and file array
        integer,intent(in)::mo_nener
        real(8),intent(in)::mo_ener(mo_nener)
        integer,intent(out)::mo_ener_pointer(mo_nener)

        integer::i

        do i=1,mo_nener
            call get_index(mo_ener(i),file_nener,file_ener,mo_ener_pointer(i))
        enddo
        return
    end subroutine 

    subroutine get_temp_pointer(mo_ntau,mo_temp,mo_temp_pointer)
        implicit none
        integer,intent(in)::mo_ntau
        real(8),intent(in)::mo_temp(mo_ntau)
        integer,intent(out)::mo_temp_pointer(mo_ntau)

        real(8),parameter::mec2 = 511.d0,boltz = 8.617333262145d-8
        real(8) theta
        integer i

        do i=1,mo_ntau
            theta = mo_temp(i)*boltz/mec2
            call get_index(theta,mo_ntau,file_temper,mo_temp_pointer(i))
        enddo
        return
    end subroutine


end module redistribution