
subroutine readinput(lun11)
    use Reflection_Var
    use constants
    implicit none
    integer,intent(in)::lun11
    character(80)::filein
    integer::io_status,newunit
    namelist /input_params/&
    &   rfnmaxmain,rfnmaxrte,&
    &   ndrt,nfrt,nmurt,&
    &   rfppemin,rfppemax,rfppemax2,&
    &   rftaumin,rftaumaxrt,&
    &   rfmumin,rfmumax,ity,&
    &   rfdebug,    &

    &   rfit_temp,rftemp_gas,rftemp_gas_unit,&
    &   rfgamma,rfhcut,rfinmu,rfstinci,rfinci_file,&
    &   rfsbot,rfktbb,rfFx_frac,&

    &   rfzeta,rfnh,&
    
    &   rfop_spec,rfop_temp,rfop_inte,rfop_log,rfcomp_file,rfop_emis,&
    &   rfdataenv,rfmodel,rfabdfits,&
    &   rfh,rfhe,rfli,rfbe,rfba,rfc,rfn,rfo,rff,rfne,rfna,rfmg,&
    &   rfal,rfsi,rfp,rfs,rfcl,rfar,rfk,rfca,rfsca,rfti,rfva,&
    &   rfcr,rfmn,rffe,rfco,rfni,rfcu,rfzn


    CALL get_command_argument(1, filein)
    print*, ' infile: ', filein

    ! read input file
    open(newunit=newunit, file=filein, status='old', action='read', iostat=io_status)
    if (io_status /= 0) then
        print *, 'Error: Could not open file ', trim(filein)
        stop
    end if

    read(newunit, NML=input_params, IOSTAT=io_status)
    close(newunit)

    call getlunx(lun11) 
    open(unit=lun11,file=rfop_log,status='unknown') 

    
    write(lun11,'(A)') '--- Input Parameters Summary ---'
    write(lun11,'(A, T75,I6)') "Number of mu points: ",nmurt
    write(lun11,'(A, T75,2I6)') 'Number of iterations:', rfnmaxmain,rfnmaxrte
    write(lun11,'(A, T75,I6)') "Debug: ",rfdebug
    write(lun11,'(A, T75,I6)') "Thermal equilibrium iteration steps: ",rfit_temp
    write(lun11,'(A, T75,I6)') "Angle Grids Typr [1: GS point 2: Linear]: ",ity
    write(lun11,'(A, T75,A)') "Top illumination type: ",rfstinci
    if (rfstinci.ne.'file') then
        write(lun11,'(A,T75,F3.1,1X,E7.1,1X,F4.2)') 'Gamma | Ecut[eV] | angle[cos]',rfgamma,rfhcut,rfinmu
    else
        write(lun11,'(A,T75,A)') "Top illumination file: ",rfinci_file
    endif
    write(lun11,'(A,T75,I6)') "Bottom illumination: ",rfsbot
    write(lun11,'(A,T75,F6.2)') "Bottom illumination temperature [eV]: ",rfktbb

    write(lun11,'(A,T75,F6.2)') "F(top) / F(total)",rfFx_frac

    write(lun11,'(A, T75, F3.1)')  'Ionization parameter:', rfzeta
    write(lun11,'(A, T75, E7.1)')  'Hydrogen number density (cm^-3):', rfnh

    write(lun11,'(A, T75, A)') "Compton redistribution file: ",trim(rfcomp_file)

    write(lun11,'(A, T75, A)') "Atomic Database path: ",rfdataenv

    write(lun11,'(A, T75, 3A)') "Final output: ",rfop_spec,rfop_temp,rfop_inte

    write(lun11,'(A, T75, I6)') "Model [0: reflected model    1: pure scattering model]: ",rfmodel 

    if(rfmodel.eq.1) then
        write(lun11,'(A, T75,E7.1)') "Gas temperature: ",rftemp_gas
        write(lun11,'(A, T75,A)') "Unit of gas temperature: ",rftemp_gas_unit
    endif
        
    

    flush(lun11)
    ! convert initial temperature to K
    if (rftemp_gas_unit.eq.'K'.or.rftemp_gas_unit.eq.'k') then
        rftemp_gas = rftemp_gas
    else if (rftemp_gas_unit.eq.'ev') then
        rftemp_gas = rftemp_gas/boltz
    else if (rftemp_gas_unit.eq.'kev') then
        rftemp_gas = rftemp_gas*1e3/boltz
    else if (rftemp_gas_unit.eq.'erg') then
        rftemp_gas = rftemp_gas/ergsev/boltz
    endif

end subroutine readinput
    