
subroutine readinput(lun11)
    use Reflection_Var
    use constants
    implicit none
    integer,intent(in)::lun11
    character(80)::filein
    character(len=80) :: separator
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
    &   rfsbot,rfktbb,rfFx_frac,ktb_nthcomp,&

    &   rfzeta,rfnh,&
    
    &   rfop_spec,rfop_temp,rfop_inte,rfop_log,rfcomp_file,rfop_emis,&
    &   rfdataenv,rfmodel,rfabdfits,&
    &   rfh,rfhe,rfli,rfbe,rfba,rfc,rfn,rfo,rff,rfne,rfna,rfmg,&
    &   rfal,rfsi,rfp,rfs,rfcl,rfar,rfk,rfca,rfsca,rfti,rfva,&
    &   rfcr,rfmn,rffe,rfco,rfni,rfcu,rfzn

    call default_input

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

    
    separator = '----------------------------------------------------------------------'
    write(lun11,'(A)') separator
    write(lun11,'(A)') '                     SIMULATION PARAMETERS REPORT'
    write(lun11,'(A)') separator
    write(lun11,*) ! 空行

    ! --- INPUT PARAMETERS SUMMARY ---
    write(lun11,'(A)') '--- [Input Parameters Summary] ---'
    write(lun11,'(A55, T60, ": ", I8)') "Debug Mode", rfdebug
    write(lun11,'(A55, T60, ": ", I8)') "Model [0: Reflected, 1: Pure Scattering]", rfmodel 
    write(lun11,*) 

    ! --- GRIDS and ITERATIONS ---
    write(lun11,'(A)') '--- [Grids and Iterations] ---'
    write(lun11,'(A55, T60, ": ", I8)') "Max Main Iterations", rfnmaxmain
    write(lun11,'(A55, T60, ": ", I8)') "Max RTE Iterations", rfnmaxrte
    write(lun11,'(A55, T60, ": ", I8)') "Number of mu points", nmurt
    write(lun11,'(A55, T60, ": ", I8)') "Number of energy points", nfrt
    write(lun11,'(A55, T60, ": ", I8)') "Number of angle points", nmurt
    write(lun11,'(A55, T60, ": ", I8)') "Angle Grids Type [1: GS point, 2: Linear]", ity
    write(lun11,*)

    ! --- PLASMA PHYSICS ---
    write(lun11,'(A)') '--- [Plasma Physics] ---'
    write(lun11,'(A55, T60, ": ", I8)')     "Thermal equilibrium steps", rfit_temp
    write(lun11,'(A55, T60, ": ", ES12.4)') "Initial Gas Temperature [K]", rftemp_gas
    write(lun11,'(A55, T60, ": ", ES12.4)') "Hydrogen Number Density [cm^-3]", rfnh
    write(lun11,'(A55, T60, ": ", F12.4)')  "Ionization Parameter (log10)", rfzeta
    write(lun11,*)

    ! --- ILLUMINATION ---
    write(lun11,'(A)') '--- [Illumination Configuration] ---'
    write(lun11,'(A55, T60, ": ", A)')      "Top Illumination Type", trim(rfstinci)
    if (rfstinci.eq."file") then
        write(lun11,'(A)') "Custm Incident Path:"
        write(lun11,'(T5, A)') trim(rfinci_file)
    else if ((rfstinci.eq.'powlaw').or.(rfstinci.eq.'cutoff')) then
        write(lun11,'(A55, T60, ": ", F12.4)')  "Gamma", rfgamma
        write(lun11,'(A55, T60, ": ", F12.4)')  "Ecut[eV]", rfhcut
    else if (rfstinci.eq.'blackbody') then
        write(lun11,'(A55, T60, ": ", F12.4)')  "Blackbody Temperature", rfgamma
    else if (rfstinci.eq."nthcomp") then
        write(lun11,'(A55, T60, ": ", F12.4)')  "Gamma", rfgamma
        write(lun11,'(A55, T60, ": ", F12.4)')  "kT_E[eV]", rfhcut
        write(lun11,'(A55, T60, ": ", F12.4)')  "kT_bb[eV]",ktb_nthcomp
    endif
    write(lun11,'(A55, T60, ": ", F12.4)')  "Incidence Angle (cosine)", rfinmu
    write(lun11,'(A55, T60, ": ", I8)')     "Bottom Illumination Switch", rfsbot
    write(lun11,'(A55, T60, ": ", F12.4)')  "F(top) / F(total) Ratio", rfFx_frac   
    write(lun11,'(A55, T60, ": ", F12.4)')  "Bottom Illumination Temp [eV]", rfktbb
    write(lun11,*)

    ! --- FILE PATHS ---
    write(lun11,'(A)') '--- [File Paths & Environment] ---'
    write(lun11,'(A)') "Compton Redistribution File:"
    write(lun11,'(T5, A)') trim(rfcomp_file)
    
    write(lun11,'(A)') "Atomic Database Path:"
    write(lun11,'(T5, A)') trim(rfdataenv)


    write(lun11,'(A)') "Output Files (Spec, Temp, Inte, Fits):"
    write(lun11,'(T5, A)') trim(rfop_spec)
    write(lun11,'(T5, A)') trim(rfop_temp)
    write(lun11,'(T5, A)') trim(rfop_inte)
    write(lun11,'(T5, A)') trim(rfabdfits)
    
    write(lun11,'(A)') separator
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
    