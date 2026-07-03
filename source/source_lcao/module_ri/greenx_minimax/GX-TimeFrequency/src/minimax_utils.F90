! ***************************************************************************************************
!  Copyright (C) 2020-2022 Green-X library                                                          
!  This file is distributed under the terms of the APACHE2 License.                                 
!                                                                                                   
! ***************************************************************************************************
!> \brief This module contains auxiliary procedures and data structures for the main minimax routines
! ***************************************************************************************************
module minimax_utils
#include "gx_common.h"
  use kinds, only: dp
  implicit none

  private

  type :: er_aw_aux
     ! Sorted array of the energy ranges
     real(kind=dp), dimension(:), allocatable :: energy_range
     ! Matrices with coefficients and weights per energy region
     real(kind=dp), dimension(:, :), allocatable :: aw_erange_matrix
  end type er_aw_aux

  !> Transformation types
  integer, parameter, public :: cosine_tw = 1
  integer, parameter, public :: cosine_wt = 2
  integer, parameter, public :: sine_tw = 3

  public :: find_erange, er_aw_aux

contains

  !> \brief Find first element in unsorted array that is strictly greater than a given value
  !>        This algorithm is O(n), difficult to do better with unsorted arrays
  !! @param[in] lenght - lenght of sorted array
  !! @param[in] einter - sorted array of the energy intervals
  !! @param[in] eval - the energy value
  function find_erange(length, einter, eval) result(idx)
    integer, intent(in)                     :: length
    real(dp), dimension(length), intent(in) :: einter
    real(dp), intent(in)                    :: eval
    integer                                 :: idx

    ! Auxiliary variables
    integer                                 :: jdx
    real(dp)                                :: tmp_min_max

    ! Begin work
    tmp_min_max = huge(0.0_dp)
    idx = length + 1

    do jdx = 1, length
       if (eval < einter(jdx) .and. einter(jdx) < tmp_min_max) then
          idx = jdx
          tmp_min_max = einter(jdx)
       end if
    end do

  end function find_erange

end module minimax_utils
