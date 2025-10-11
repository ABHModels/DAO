subroutine ener(epi,ncn2,df) 
!                                                                       
!     Name: ener.f90  
!     Description:  
!     This routine sets up the energy grid
!     Grid is logarithmic with two subranges:  0.1eV-40 keV, 40keV-1MeV.
!     This structure of epi is key to operation of various other routines
!     author: T. Kallman                                                
!     List of Parameters:
!           Output:
!           epi(ncn)=energy grid (ev)
!           ncn2=length of epi
!     Dependencies:  none
!     Called by:  xstar
!                                                                       
      implicit none 
!                                                                       
      real(8) epi(ncn2) ,df(ncn2)
      integer numcon,numcon2,numcon3,ncn2,ll,ll2 
      real(8) ebnd1,ebnd2,ebnd2o,dele 
!                                                                       
      numcon = ncn2 
    !   if (numcon.gt.ncn) stop 'ncn2 too large for arry dimension' 
    !   if (numcon.lt.4) write (6,*) 'in ener: numcon error' 
      numcon2=max(2,ncn2/50) 
      numcon3=numcon-numcon2 
      ebnd1=0.1 
!     nb changed energy grid for H only                                 
      ebnd2=4.e+5 
!      ebnd2=4.e+1                                                      
      ebnd2o=ebnd2 
      dele=(ebnd2/ebnd1)**(1./float(numcon3-1)) 
      epi(1)=ebnd1 
!      write (lun11,*)'in ener',ncn2,numcon,numcon2,numcon3                 
      do ll=2,numcon3 
        epi(ll)=epi(ll-1)*dele 
        enddo 
      ebnd2=1.e+6 
      ebnd1=ebnd2o 
      dele=(ebnd2/ebnd1)**(1./float(numcon2-1)) 
      do ll2=1,numcon2 
        ll=ll2+numcon3 
        epi(ll)=epi(ll-1)*dele 
        enddo 

    df(1)=0.d0
    do ll=2,ncn2
        df(ll) = epi(ll) - epi(ll-1)
    enddo
    df(1)=df(2)  !! Need to check this!
!                                                                       
      return 
      end  