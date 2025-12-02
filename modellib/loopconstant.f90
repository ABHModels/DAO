     subroutine loopconstant(trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,          &
        & numrec,npass,nlimd,rmax,xpxcol,xi,zeta,lfix,                     &
        & lun11,abel,cfrac,emult,taumax,xeemin,spectype,specfile,specunit, &
        & kmodelname,nloopctl,critf,vturbi,eptmp,zrtmp,numcon2,ncn2,radexp)
        ! read depth invariant parameters
        use Reflection_Var
        use globaldata
         implicit none                                                           
         real(8) eptmp(ncn),zrtmp(ncn),abel(nl),abel2(30) 
         character(8) stringst,kblnk8 
         character(80) specfile,spectype,stringsl,kblnk80,stringst2 
         character(30) kmodelname 
         integer nloopctl,specunit,ierr,ll,lcdd2,lun13,nenergy,ncn2 
         integer lwri,lpri,lcdd,numrec,npass,nlimd,lfix,lun11,numcon2,mm 
         real(8) trad,xlum,r,t,xpx,p,rmax,xpxcol,xi,zeta,cfrac,emult,taumax,&
        &     xeemin,critf,vturbi,ccc,xlum2,xpcol,r19,radexp                                                       
     
!       density                                                     
        xpx = rfnh

        numcon2 = nfrt
                                                                                                     
!       number of steps                                             
        numrec = ndrt

!       number of iterations                                        
        nlimd = rfit_temp
                                                      
!       abundances                                                  
        do ll=1,nl
        abel2(ll)=rfabel(ll)
        enddo
        do mm=1,nl 
            abel(mm)=abel2(mm)
        enddo 

        ncn2 = numcon2
        print*, " Xstar constant parameters done"
        return 
    end subroutine loopconstant