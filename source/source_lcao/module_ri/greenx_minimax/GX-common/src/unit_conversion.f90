!> Unit Conversion, defined according to [CODATA 2018])(http://physics.nist.gov/constants) 
!> 
!> For all units, please:
!> a) Define using CODATA 2018
!> b) Define using double precision (dp). For example, elec_mass = 9.10938370e-31_dp;
!> c) Use descriptive naming
!>
!> See the 20th May 2019 redefinition of the SI units and CODATA 2018.
!> Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2019)
!> "The 2018 CODATA Recommended Values of the Fundamental Physical Constants"
!> (Web Version 8.0). Database developed by J. Baker, M. Douma, and S. Kotochigova.
!> Available at http://physics.nist.gov/constants,
!> National Institute of Standards and Technology, Gaithersburg, MD 20899.
!
module unit_conversion
    use kinds, only: dp 
    implicit none
    private

    !> Bohr to angstrom                            
    real(kind=dp), parameter, public :: bohr_to_angstrom = 0.529177210903_dp

end module
