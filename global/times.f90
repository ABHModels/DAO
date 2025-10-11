 module times
      integer ntp
      parameter (ntp=99)
      real tread,tloop,tfunc,trates1,thcor,trates2,tucalc(ntp),theat
      integer ncall(ntp)
end module times