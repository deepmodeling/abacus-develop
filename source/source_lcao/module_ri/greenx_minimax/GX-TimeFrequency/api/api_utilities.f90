module api_utilites
    use error_handling, only: error_message__
    use minimax_omega,    only : omega_npoints_supported
    use minimax_tau, only : tau_npoints_supported
    implicit none
    private
    public :: gx_check_ntau, gx_get_error_message

contains

   !> \brief Check whether ntau is amongst the list of supported mesh-sizes.
   !> Return exit status in ierr and error message in msg (if any)
   subroutine gx_check_ntau(ntau, msg, ierr)

      integer, intent(in) :: ntau
      character(len=*), intent(out) :: msg
      integer, intent(out) :: ierr

      integer :: i
      character(len=1024) :: my_sbuf

      ierr = 1
      msg = "No error detected"

      if (size(omega_npoints_supported) /= size(tau_npoints_supported)) then
         msg = "Internal incosistency: size(omega_npoints_supported) /= size(tau_npoints_supported))"
         return
      end if

      do i = 1, size(omega_npoints_supported)
         if (omega_npoints_supported(i) /= tau_npoints_supported(i)) then
            msg = "Internal incosistency: omega_npoints_supported(ii) /= tau_npoints_supported(ii))"
         end if
         if (ntau == tau_npoints_supported(i)) then
            ierr = 0
            msg = "Input ntau found in internal list"
            exit
         end if
      end do

      if (ierr /= 0) then
         write (my_sbuf, "(a,i0,3a,*(i0,1x))") &
            "Cannot find ntau: ", ntau, " in the internal list.", char(10), &
            "Please choose among:", tau_npoints_supported
         msg = my_sbuf(1:min(len_trim(my_sbuf), len(msg)))
      end if

   end subroutine gx_check_ntau

   !> \brief Return error_message in the msg buffer provided by caller
   !>  This routine should be invoked to access the internal error message
   !>  if one of the GreenX APIs returns ierr /= 0.
   subroutine gx_get_error_message(msg)
      character(len=*), intent(out) :: msg
      msg = error_message__(1:min(len_trim(error_message__), len(msg)))
   end subroutine gx_get_error_message

end module
