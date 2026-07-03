module gx_minimax_wrp
  use iso_c_binding, only: c_int, c_double, c_ptr, c_f_pointer
  use gx_minimax, only: gx_minimax_grid_frequency
  use kinds, only: dp
  implicit none

  private
  public :: gx_minimax_grid_frequency_wrp

contains

  subroutine gx_minimax_grid_frequency_wrp(num_points, e_min, e_max, omega_points, omega_weights, ierr) &
      bind(C, name="gx_minimax_grid_frequency_wrp")
    integer(kind=c_int), intent(in) :: num_points
    real(kind=c_double), intent(in) :: e_min, e_max
    type(c_ptr), value :: omega_points
    type(c_ptr), value :: omega_weights
    integer(kind=c_int), intent(out) :: ierr

    real(dp), allocatable :: omega_points_f(:)
    real(dp), allocatable :: omega_weights_f(:)
    real(kind=c_double), pointer :: omega_points_fptr(:)
    real(kind=c_double), pointer :: omega_weights_fptr(:)
    integer :: ig

    call c_f_pointer(omega_points, omega_points_fptr, [num_points])
    call c_f_pointer(omega_weights, omega_weights_fptr, [num_points])

    call gx_minimax_grid_frequency(num_points, e_min, e_max, omega_points_f, omega_weights_f, ierr)
    if (ierr == 0) then
      do ig = 1, num_points
        omega_points_fptr(ig) = omega_points_f(ig)
        omega_weights_fptr(ig) = omega_weights_f(ig)
      end do
    end if

    if (allocated(omega_points_f)) deallocate(omega_points_f)
    if (allocated(omega_weights_f)) deallocate(omega_weights_f)
  end subroutine gx_minimax_grid_frequency_wrp

end module gx_minimax_wrp
