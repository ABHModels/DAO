module constants
!
!     This module contains definitions of natural constants
!     and unit conversion routines. Base unit system is cgs. 
!
!     Source:    
!     NIST 2018 CODATA recommended values of physical constants
!     https://pml.nist.gov/cuu/Constants/index.html      
!        
      implicit none 

      real(8), parameter :: sigth = 6.6524587321e-25   ! Thomson cross section (cm^2)
      real(8), parameter :: ergsev = 1.602176634e-12   ! erg per eV 
      real(8), parameter :: hev = 4.1356655385381E-15  ! ev per hz  
      real(8), parameter :: herg = 6.62606957030463E-27! erg per hz

      real(8), parameter :: emmh = 1.66053906660e-24   ! Atomic mass constant (g)
      real(8), parameter :: bk = 1.380649e-16          ! Boltzmann constant (erg K^-1)
      real(8), parameter :: ccc = 2.99792458e10        ! Speed of light (cm s^-1)
      real(8), parameter :: pi = 3.141592653589793     !
      real(8), parameter :: fourpi = 4*pi              !
      real(8), parameter :: Ry = 13.605693122994       ! Rydberg constant times hc (eV)
      real(8), parameter :: hh =  6.62607015e-27        ! Planck constant (erg s)
      real(8), parameter :: me = 9.1093837015e-28      ! Electron rest mass (g)
      real(8), parameter :: spol = 2.99792458e10       ! speed of light in cm/s
      real(8), parameter :: boltz = 8.61733326e-5     ! Boltzmann constant (eV K^-1)
      real(8), parameter :: mec2 = 511.d3             ! Electron rest mass (eV)
!    
      !data ergsev/1.602176634-12/
      !data sigth/6.6524587321e-25/                                                 
      !data emmh/1.66053906660-24/                                                 
      !data bk/1.380649e-16/ 
      !data pi/3.14159265358979323846/ 
      !data ccc/2.99792458e10/
!
!
      contains
!        
      function A2keV(x)
!       Convert wavelenght in Angstrom to energy in keV 
            implicit none
            real(8) A2keV
            real(8) :: x
            A2keV=1e5*hh*ccc/ergsev/x                        
      end function A2keV
      
      function keV2A(x)
!       Convert energy in keV to wavelenght in Angstrom
            real(8) :: keV2A
            real(8) :: x
            keV2A=A2keV(x)
      end function keV2A
      
      function A2eV(x)
!       Convert wavelenght in Angstrom to energy in eV
            implicit none
            real(8):: A2eV
            real(8) :: x
            A2eV=1e3*A2keV(x)
      end function A2eV
      
      function eV2A(x)
!       Convert energy in eV to wavelenght in Angstrom
            implicit none
            real(8) :: eV2A
            real(8) :: x
            eV2A=keV2A(1e-3*x)
      end function eV2A

      function erg2eV(x)
!       Convert energies from erg to eV
            implicit none
            real(8) :: erg2eV
            real(8) :: x
            erg2eV=ergsev*x
      end function erg2eV

      function eV2erg(x)
!       Convert energies from eV to erg
            implicit none
            real(8) :: eV2erg
            real(8) :: x
            eV2erg=x/ergsev
      end function eV2erg
                                                      
end module constants