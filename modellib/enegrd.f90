subroutine enegrd
  use Reflection_Var
  implicit none
  integer::nmaxp
  real(8),parameter::ppemin=0.1, ppemax=9.9d5, ppemax2=1d6
  integer::numcon, numcon2, numcon3, ll, ll2
  real(8):: ebnd1, ebnd2, ebnd2o, dele

  nmaxp = nfrt
  numcon = nfrt
  if (numcon.lt.4) stop 'in ener: numcon error'
  numcon2=max(2,nmaxp/50)
  numcon3=numcon-numcon2
  ebnd1=ppemin
  ebnd2=ppemax
  ebnd2o=ebnd2
  dele=(ebnd2/ebnd1)**(1./dfloat(numcon3-1))
  rfener(1)=ebnd1
  do ll=2,numcon3
    rfener(ll)=rfener(ll-1)*dele
  enddo
  ebnd2=ppemax2
  ebnd1=ebnd2o
  dele=(ebnd2/ebnd1)**(1./dfloat(numcon2-1))
  do ll2=1,numcon2
    ll=ll2+numcon3
    rfener(ll)=rfener(ll-1)*dele
  enddo

  wrfener(1)=0.d0
  do ll=2,nmaxp
    wrfener(ll) = rfener(ll) - rfener(ll-1)
  enddo
  wrfener(1)=wrfener(2)  !! Need to check this!
end subroutine


subroutine enegrd2
    ! This is XSTAR energy grids 
    use Reflection_Var
    implicit none           
    real(8) epi(nfrt),df(nfrt)                                                           
    integer numcon,numcon2,numcon3,ncn2,ll,ll2 
    real(8) ebnd1,ebnd2,ebnd2o,dele 
          
    ncn2 = nfrt                                                      
    numcon = ncn2 

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

    do ll=1,ncn2
      rfener(ll) = epi(ll)
      wrfener(ll)= df(ll)
    enddo
    return 

end subroutine