! **************************************************************************************************
!  Copyright (C) 2020-2022 Green-X library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!> for correlation methods as well as weights for the inhomogeneous cosine/sine transform.
!>
!> NB: When dealing with runtime exceptions, we set ierr to a non-zero value and return immediately
!  to the caller so we don't need to goto to a cleanup section at the end of the procedure.
!  Assume -std=f2008: i.e. allocatable arrays are automatically deallocated when going out of scope.
! **************************************************************************************************

module minimax_grids
#include "gx_common.h"
  use kinds,             only: dp
  use error_handling,    only: register_exc
  use constants,         only: pi
  use minimax_tau,       only: get_points_weights_tau
  use minimax_omega,     only: get_points_weights_omega
  use minimax_utils,     only: cosine_wt, cosine_tw, sine_tw

  implicit none

  private

  !> Main entry point for client code.
  public :: gx_minimax_grid, gx_minimax_grid_frequency, gx_minimax_grid_rpa

contains

  !> \brief Compute minimax grid for GW calculation on imaginary time/frequency domain.
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] e_min: Minimum transition energy. Arbitrary units as we only need e_max/e_min
  !! @param[in] e_max: Maximum transition energy.
  !! @param[out] tau_points: imaginary time grid points
  !! @param[out] tau_weights: weights for imaginary time grid points weights
  !! @param[out] omega_points: imaginary frequency grid points
  !! @param[out] omega_weights: weights for imaginary frequency grid points
  !! @param[out] cosft_wt: weights for tau -> omega cosine transform. cos(w*t) factor is included.
  !! @param[out] cosft_tw: weights for omega -> tau cosine transform. cos(w*t) factor is included.
  !! @param[out] sinft_wt: weights for tau -> omega sine transform. sin(w*t) factor is included.
  !! @param[out] max_errors: Max error for the three kind of transforms (same order as previous args)
  !! @param[out] cosft_duality_error. Max_{ij} |AB - I| where A and B are the cosft_wt and cosft_tw matrices.
  !! @param[out] ierr: Exit status
  subroutine gx_minimax_grid(num_points, e_min, e_max, &
       tau_points, tau_weights, omega_points, omega_weights, &
       cosft_wt, cosft_tw, sinft_wt, &
       max_errors, cosft_duality_error, ierr, bare_cossine_weights)

    integer, intent(in)                               :: num_points
    real(kind=dp), intent(in)                         :: e_min, e_max
    real(kind=dp), allocatable, dimension(:), &  
         intent(out)                                  :: tau_points, tau_weights
    real(kind=dp), allocatable, dimension(:), &
         intent(out)                                  :: omega_points(:), omega_weights(:)
    real(kind=dp), allocatable, dimension(:, :), &
         intent(out)                                  :: cosft_wt(:, :), cosft_tw(:, :), sinft_wt(:, :)
    real(kind=dp), intent(out)                        :: max_errors(3), cosft_duality_error
    logical, intent(in), optional                     :: bare_cossine_weights
    integer, intent(out)                              :: ierr

    ! Internal variables
    logical                                           :: my_bare_cossine_weights
    integer, parameter                                :: cos_t_to_cos_w = 1
    integer, parameter                                :: cos_w_to_cos_t = 2
    integer, parameter                                :: sin_t_to_sin_w = 3
    integer                                           :: i_point, j_point
    real(kind=dp)                                     :: e_range, scaling
    real(kind=dp), dimension(:), allocatable          :: x_tw
    real(kind=dp), dimension(:, :), allocatable       :: mat

    my_bare_cossine_weights = .false.
    if (present(bare_cossine_weights)) then
      if(bare_cossine_weights) my_bare_cossine_weights = .true.
    endif

    ! Begin work
    e_range = e_max/e_min   
    ierr = 0

    ! Allocations
    allocate (x_tw(2*num_points))
    if (.not. allocated(omega_points)) then
       allocate (omega_points(num_points))
    end if
    if (.not. allocated(omega_weights)) then
       allocate (omega_weights(num_points))
    end if
    if (.not. allocated(tau_points)) then
       allocate (tau_points(num_points))
    end if
    if (.not. allocated(tau_weights)) then
       allocate (tau_weights(num_points))
    end if

    ! Get the frequency grid points and weights 
    call get_points_weights_omega(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! Scale the frequency grid points and weights from [1,R] to [e_min,e_max]
    ! Note: the frequency grid points and weights include a factor of two
    scaling = e_min
    omega_points = x_tw(1: num_points) *scaling
    omega_weights = x_tw(num_points+1: 2* num_points) *scaling   

    ! Get the time grid points and weights
    call get_points_weights_tau(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! Scale the time grid points and weights from [1,R] to [e_min,e_max]
    scaling = 2.0_dp *e_min
    tau_points = x_tw(1: num_points) /scaling
    tau_weights = x_tw(num_points+1: 2* num_points) /scaling

    if (.not. allocated(cosft_wt)) allocate (cosft_wt(num_points, num_points))
    if (.not. allocated(cosft_tw)) allocate (cosft_tw(num_points, num_points))
    if (.not. allocated(sinft_wt)) allocate (sinft_wt(num_points, num_points))

    ! get the weights for the cosine transform W^c(it) -> W^c(iw)
    call get_transformation_weights(num_points, tau_points, omega_points, cosft_wt, e_min, e_max, &
         max_errors(1), cos_t_to_cos_w, ierr)
    if (ierr /= 0) return

    ! get the weights for the cosine transform W^c(iw) -> W^c(it)
    call get_transformation_weights(num_points, tau_points, omega_points, cosft_tw, e_min, e_max, &
         max_errors(2), cos_w_to_cos_t, ierr)
    if (ierr /= 0) return

    ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
    call get_transformation_weights(num_points, tau_points, omega_points, sinft_wt, e_min, e_max, &
         max_errors(3), sin_t_to_sin_w, ierr)
    if (ierr /= 0) return

    ! Compute the actual weights used for the inhomogeneous cosine/ FT and check whether
    ! the two matrices for the forward/backward transform are the inverse of each other.
    if(.not.my_bare_cossine_weights) then
      do j_point = 1, num_points
         do i_point = 1, num_points
            cosft_wt(j_point, i_point) = cosft_wt(j_point, i_point)*cos(tau_points(i_point)*omega_points(j_point))
            cosft_tw(i_point, j_point) = cosft_tw(i_point, j_point)*cos(tau_points(i_point)*omega_points(j_point))
            sinft_wt(j_point, i_point) = sinft_wt(j_point, i_point)*sin(tau_points(i_point)*omega_points(j_point))
         end do
      end do
    end if

    allocate (mat(num_points, num_points))
    mat = matmul(cosft_wt, cosft_tw)
    do i_point = 1, num_points
       mat(i_point, i_point) = mat(i_point, i_point) - 1.0_dp
    end do
    cosft_duality_error = maxval(abs(mat))

    deallocate (mat)
    deallocate (x_tw)

  end subroutine gx_minimax_grid

  !> \brief Compute minimax grid for RPA correlation energy on imaginary time/frequency domain.
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] e_min: Minimum transition energy. Arbitrary units as we only need e_max/e_min
  !! @param[in] e_max: Maximum transition energy.
  !! @param[out] tau_points: imaginary time grid points
  !! @param[out] tau_weights: weights for imaginary time grid points weights
  !! @param[out] omega_points: imaginary frequency grid points
  !! @param[out] omega_weights: weights for imaginary frequency grid points
  !! @param[out] cosft_wt: weights for tau -> omega cosine transform. cos(w*t) factor is included.
  !! @param[out] max_errors: Max error for transforms (same order as previous args)
  !! @param[out] ierr: Exit status  
  subroutine gx_minimax_grid_rpa(num_points, e_min, e_max, &
       tau_points, tau_weights, omega_points, omega_weights, &
       cosft_wt,max_errors,ierr)

    integer, intent(in)                               :: num_points
    real(kind=dp), intent(in)                         :: e_min, e_max
    real(kind=dp), allocatable, dimension(:), &
         intent(out)                                  :: tau_points, tau_weights
    real(kind=dp), allocatable, dimension(:), &
         intent(out)                                  :: omega_points(:), omega_weights(:)
    real(kind=dp), allocatable, dimension(:, :), &
         intent(out)                                  :: cosft_wt(:, :)
    real(kind=dp), intent(out)                        :: max_errors
    integer, intent(out)                              :: ierr

    ! Internal variables
    integer, parameter                                :: cos_t_to_cos_w = 1
    ! integer                                           :: i_point, j_point
    real(kind=dp)                                     :: e_range, scaling
    real(kind=dp), dimension(:), allocatable          :: x_tw

    ! Begin work
    e_range = e_max/e_min
    ierr = 0

    ! Allocations
    allocate (x_tw(2*num_points))
    if (.not. allocated(omega_points)) then
       allocate (omega_points(num_points))
    end if
    if (.not. allocated(omega_weights)) then
       allocate (omega_weights(num_points))
    end if
    if (.not. allocated(tau_points)) then
       allocate (tau_points(num_points))
    end if
    if (.not. allocated(tau_weights)) then
       allocate (tau_weights(num_points))
    end if

    ! Get the frequency grid points and weights 
    call get_points_weights_omega(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! Scale the frequency grid points and weights from [1,R] to [e_min,e_max]
    ! Note: the frequency grid points and weights include a factor of two
    scaling = e_min
    omega_points = x_tw(1: num_points) *scaling
    omega_weights = x_tw(num_points+1: 2* num_points) *scaling

    ! Get the time grid points and weights
    call get_points_weights_tau(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! Scale the time grid points and weights from [1,R] to [e_min,e_max]
    scaling = 2.0_dp *e_min
    tau_points = x_tw(1: num_points) /scaling
    tau_weights = x_tw(num_points+1: 2* num_points) /scaling

    allocate (cosft_wt(num_points, num_points))

    ! get the weights for the cosine transform W^c(it) -> W^c(iw)
    call get_transformation_weights(num_points, tau_points, omega_points, cosft_wt, e_min, e_max, &
         max_errors, cos_t_to_cos_w, ierr)
    if (ierr /= 0) return

    deallocate (x_tw)

  end subroutine gx_minimax_grid_rpa

  !> \brief Retrieves the frequency grid for a canonical GW/RPA calculation
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] e_min: Minimum transition energy. Arbitrary units as we only need e_max/e_min
  !! @param[in] e_max: Maximum transition energy.
  !! @param[out] omega_points: imaginary frequency grid points
  !! @param[out] omega_weights: weights for imaginary frequency grid points
  !! @param[out] ierr: Exit status
  subroutine gx_minimax_grid_frequency (num_points, e_min, e_max, omega_points, omega_weights, ierr)
    integer, intent(in)                               :: num_points
    real(kind=dp), intent(in)                         :: e_min, e_max
    real(kind=dp), allocatable, dimension(:), &
         intent(out)                                  :: omega_points(:), omega_weights(:)
    integer, intent(out)                              :: ierr

    ! Internal variables
    real(kind=dp)                                     :: e_range, scaling
    real(kind=dp), dimension(:), allocatable          :: x_tw
    
    ! Begin work
    e_range = e_max/e_min   
    ierr = 0

    ! Allocations
    allocate (x_tw(2*num_points))
    if (.not. allocated(omega_points)) then
       allocate (omega_points(num_points))
    end if
    if (.not. allocated(omega_weights)) then
       allocate (omega_weights(num_points))
    end if

    ! Get the frequency grid points and weights 
    call get_points_weights_omega(num_points, e_range, x_tw, ierr)
    if (ierr /= 0) return

    ! Scale the frequency grid points and weights from [1,R] to [e_min,e_max]
    ! Note: the frequency grid points and weights include a factor of two
    scaling = e_min
    omega_points = x_tw(1: num_points) *scaling
    omega_weights = x_tw(num_points+1: 2* num_points) *scaling

    deallocate (x_tw)
    
  end subroutine gx_minimax_grid_frequency


  !> \brief Get the weights eiter for the cosine/sin transformation for tau to omega or viceversa
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] tau_points: imaginary time grid points
  !! @param[in] omega_points: imaginary frequency grid points
  !! @param[inout] weights: corresponding tranformation weights
  !! @param[in] e_min: Minimum transition energy.
  !! @param[in] e_max: Maximum transition energy.
  !! @param[inout] max_errors: Max error for the three kind of transforms
  !! @param[in] transformation type : 1 the cosine transform cos(it) -> cos(iw)
  !!                                : 2 the cosine transform cos(iw) -> cos(it)
  !!                                : 3 the sine transform   sin(it) -> sin(iw)
  !! @param[in] ierr: exit status
  subroutine get_transformation_weights(num_points, tau_points, omega_points, weights, e_min, e_max, &
       max_error, transformation_type, ierr)

    integer, intent(in)                                :: num_points
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: tau_points
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: omega_points
    real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                 :: weights
    real(kind=dp), intent(in)                          :: e_min, e_max
    real(kind=dp), intent(inout)                       :: max_error
    integer, intent(in)                                :: transformation_type
    integer, intent(out)                               :: ierr

    ! Internal variables
    integer                                            :: i_node, i_point, j_point, k_point, &
         num_x_nodes
    integer, parameter                                 :: nodes_factor = 200
    real(kind=dp)                                      :: current_point, x_factor
    real(kind=dp), allocatable, dimension(:)           :: weights_work, x_mu, psi
    real(kind=dp), allocatable, dimension(:, :)        :: mat_A

    integer                                            :: lwork
    integer, allocatable, dimension(:)                 :: iwork      
    real(kind=dp), allocatable, dimension(:)           :: vec_S, vec_UT_psi, work     
    real(kind=dp), allocatable, dimension(:, :)        :: mat_U, mat_VT, mat_VT_s

    external :: dgemm
    external :: dgesdd

    ! Begin work
    ierr = 0

    allocate (weights_work(num_points), source=0.0_dp)

    ! compute the number of x nodes per magnitude points per 10-interval
    num_x_nodes = (int(log10(e_max/e_min)) + 1)*nodes_factor

    ! make sure that the number of x nodes are at least as many integration points
    num_x_nodes = max(num_x_nodes, num_points)

    allocate (x_mu(num_x_nodes), source=0.0_dp)
    allocate (psi(num_x_nodes), source=0.0_dp)

    ! Allocations for the BLAS routines
    ! double the value nessary for 'A' to achieve good performance
    lwork = 8*num_points*num_points + 12*num_points + 2*num_x_nodes
    allocate (iwork(8*num_points), source=0)
    allocate (work(lwork), source=0.0_dp)

    allocate (mat_A(num_x_nodes, num_points), source=0.0_dp)
    allocate (mat_U(num_x_nodes, num_x_nodes), source=0.0_dp)
    allocate (mat_VT(num_x_nodes, num_points), source=0.0_dp)
    allocate (mat_VT_s(num_points, num_x_nodes), source=0.0_dp)
    allocate (vec_S(num_points), source=0.0_dp)
    allocate (vec_UT_psi(num_x_nodes), source=0.0_dp)

    ! set the x-mu logarithmically in the interval [e_min,e_max]
    x_factor = (e_max/e_min)**(1.0_dp/(real(num_x_nodes, kind=dp) - 1.0_dp))
    do i_node = 1, num_x_nodes
       x_mu(i_node) = e_min*x_factor**(i_node - 1)
    end do

    current_point = 0.0_dp
    max_error = 0.0_dp

    ! loop over all grid points
    do i_point = 1, num_points
       ! calculate psi and mat_A
       call calculate_psi_and_mat_A(num_points, tau_points, omega_points, num_x_nodes, x_mu, psi, &
            mat_A, i_point, current_point, transformation_type)

       ! Singular value decomposition of mat_A = U*Sigma*V^T
       call dgesdd('A', num_x_nodes, num_points, mat_A, num_x_nodes, vec_S, mat_U, num_x_nodes, &
            mat_VT, num_x_nodes, work, lwork, iwork, ierr)

       if (ierr /= 0) then
          _REGISTER_EXC("DGESDD returned ierr != 0")
          return
       end if

       ! integration weights = (V Sigma^-1 U^T)*psi

       ! 1) V * Sigma^-1
       do j_point = 1, num_points
          do k_point = 1, num_points
             mat_VT_s(k_point, j_point) = mat_VT(j_point, k_point)/vec_S(j_point)
          end do ! k_point
       end do ! j_point

       ! 2) (U^T)*psi
       call dgemm('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, psi, num_x_nodes, &
            0.0_dp, vec_UT_psi, num_x_nodes)

       ! 3) (V*Sigma^-1) * (U^T*psi)
       call dgemm('N', 'N', num_points, 1, num_x_nodes, 1.0_dp, mat_VT_s, num_points, vec_UT_psi, &
            num_x_nodes, 0.0_dp, weights_work, num_points)

       weights(i_point, :) = weights_work(:)

       ! calculate the maximum error of the fitting
       call calculate_max_error(num_points, tau_points, omega_points, weights_work, num_x_nodes, x_mu, &
            psi, current_point, max_error, transformation_type)
    end do ! i_point

    deallocate (x_mu, psi, mat_A, weights_work, vec_S, mat_U, mat_VT, work, iwork, mat_VT_s, vec_UT_psi)

  end subroutine get_transformation_weights

  !> \brief Calculate the auxiliary matrix for cosine/sin transformation for tau to omega or viceversa
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] tau_points: imaginary time grid points
  !! @param[in] omega_points: imaginary frequency grid points
  !! @param[in] num_x_nodes: Number of node in the interval [e_min,e_max]
  !! @param[in] x_mu : Transition energy (nodes in the interval [e_min,e_max])
  !! @param[inout] psi: corresponding auxiliary function (see transformation type definition)
  !! @param[inout] mat_A: auxiliary matrix (see transformation type definition)  
  !! @param[in] i_point: pointer for the current grid point
  !! @param[inout] current_point:  current grid point ether omega(i_point) or tau_(i_point)
  !! @param[in] transformation type :
  !!        (1) the cosine transform cos(it) -> cos(iw): psi(omega,x), mat_A = cos(omega*tau)*psi(tau,x) 
  !!        (2) the cosine transform cos(iw) -> cos(it): psi(tau,x)  , mat_A = cos(omega*tau)*psi(omega,x)
  !!        (3) the sine transform   sin(it) -> sin(iw): psi(omega,x), mat_A = sin(omega*tau)*psi(tau,x)
  subroutine calculate_psi_and_mat_A(num_points, tau_points, omega_points, num_x_nodes, x_mu, psi, &
       mat_A, i_point, current_point, transformation_type)

    integer, intent(in)                                :: num_points, num_x_nodes, i_point
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                    :: tau_points, omega_points, x_mu
    real(kind=dp), allocatable, dimension(:), &          
         intent(inout)                                 :: psi         
    real(kind=dp), allocatable, dimension(:, :), &
         intent(inout)                                 :: mat_A
    real(kind=dp), intent(inout)                       :: current_point
    integer, intent(in)                                :: transformation_type

    ! Internal variables
    integer                                            :: i_node, j_point
    real(kind=dp)                                      :: tau, omega

    ! Begin work

    ! the cosine transform cos(it) -> cos(iw)
    if (transformation_type == cosine_tw) then
       omega = omega_points(i_point)
       current_point = omega

       ! psi(omega_k,x) = 2x/(x^2+omega_k^2)
       do i_node = 1, num_x_nodes
          psi(i_node) = 2.0_dp*x_mu(i_node)/((x_mu(i_node))**2 + omega**2)
       end do

       ! mat_A = cos(omega_k * tau) psi(tau,x)
       do j_point = 1, num_points
          tau = tau_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = cos(omega*tau)*exp(-x_mu(i_node)*tau)
          end do
       end do

    ! the cosine transform cos(iw) -> cos(it)  
    else if (transformation_type == cosine_wt) then
       tau = tau_points(i_point)
       current_point = tau

       ! psi(tau_k,x) = =exp(-x*|tau_k|)
       do i_node = 1, num_x_nodes
          psi(i_node) = exp(-x_mu(i_node)*tau)
       end do

       ! mat_A = cos(tau_k,omega) psi(omega,x)
       do j_point = 1, num_points
          omega = omega_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = cos(tau*omega)*2.0_dp*x_mu(i_node)/(x_mu(i_node)**2 + omega**2)
          end do
       end do

    ! the sine transform sin(it) -> sin(iw)         
    else if (transformation_type == sine_tw) then
       omega = omega_points(i_point)
       current_point = omega

       ! psi(omega_k,x) = 2*omega_k/(x^2+omega_k^2)
       do i_node = 1, num_x_nodes
          psi(i_node) = 2.0_dp*omega/(x_mu(i_node)**2 + omega**2)
       end do

       ! mat_A = sin(omega_k,tau)*phi(tau,x)
       do j_point = 1, num_points
          tau = tau_points(j_point)
          do i_node = 1, num_x_nodes
             mat_A(i_node, j_point) = sin(omega*tau)*exp(-x_mu(i_node)*tau)
          end do
       end do
    end if

  end subroutine calculate_psi_and_mat_A

  !> /brief calculate the error of the fit function for the cosine/sin transformation for tau to omega or viceversa
  !! @param[in] num_points: Number of mesh points.
  !! @param[in] tau_points: imaginary time grid points
  !! @param[in] omega_points: imaginary frequency grid points
  !! @param[in] weights_work: work vector for the transformation weights
  !! @param[in] num_x_nodes: Number of node in the interval [e_min,e_max]
  !! @param[in] x_mu : Transition energy (nodes in the interval [e_min,e_max])
  !! @param[in] psi: corresponding auxiliary function (see transformation type definition)
  !! @param[in] current_point:  current grid point ether omega(i_point) or tau_(i_point)
  !! @param[out] max_errors: Max error for the three kind of transforms   
  !! @param[in] transformation type : 1 fit function for the cosine transform cos(it) -> cos(iw), psi(omeaga,x)
  !!                                : 2 fit function for the cosine transform cos(iw) -> cos(it), psi(tau,x)
  !!                                : 3 fit function for the sine transform   sin(it) -> sin(iw), psi(omega,x)
  subroutine calculate_max_error(num_points, tau_points, omega_points, weights_work, num_x_nodes, x_mu, &
       psi, current_point, max_error, transformation_type)

    real(kind=dp), intent(out)                       :: max_error
    real(kind=dp), intent(in)                        :: current_point      
    real(kind=dp), allocatable, dimension(:), &
         intent(in)                                  :: tau_points, omega_points, x_mu, psi, &
         weights_work
    integer, intent(in)                              :: num_points, num_x_nodes
    integer, intent(in)                              :: transformation_type      

    ! Internal variables
    integer                                          :: i_node,i_point
    real(kind=dp)                                    :: func_val, func_val_temp, max_error_tmp, &
                                                        tau, omega, x_val

    ! Begin work
    max_error_tmp = 0.0_dp

    ! the cosine transform cos(it) -> cos(iw)
    if (transformation_type == cosine_tw) then
       omega=current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0_dp
          ! calculate value of the fit function f(x) = f(x) + weights(omega)cos(omega*tau)psi(tau.x)
          do i_point = 1, num_points
             tau = tau_points(i_point) 
             func_val = func_val + weights_work(i_point)*cos(omega*tau)*exp(-x_mu(i_node)*tau)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) then
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do

    ! the cosine transform cos(iw) -> cos(it)
    else if (transformation_type == cosine_wt) then
       tau = current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0_dp
          x_val=x_mu(i_node)
          ! calculate value of the fit function f(x) = f(x) + weights(tau)cos(omega*tau)psi(omega.x)
          do i_point = 1, num_points
             omega = omega_points(i_point)
             func_val = func_val +  weights_work(i_point)*cos(tau*omega)*2.0_dp*x_val/(x_val**2 + omega**2)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) then
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do

    ! the sine transform sin(it) -> sin(iw)
    else if (transformation_type == sine_tw) then
       omega = current_point

       do i_node = 1, num_x_nodes
          func_val = 0.0_dp
          ! calculate value of the fit function f(x) = f(x) + weights(omega)sin(omega*tau)psi(tau.x)
          do i_point = 1, num_points
             tau = tau_points(i_point)
             func_val = func_val +  weights_work(i_point)*sin(omega*tau)*exp(-x_mu(i_node)*tau)
          end do

          if (abs(psi(i_node) - func_val) > max_error_tmp) then
             max_error_tmp = abs(psi(i_node) - func_val)
             func_val_temp = func_val
          end if
       end do
    end if

    if (max_error_tmp > max_error) then
       max_error = max_error_tmp
    end if

  end subroutine calculate_max_error

end module minimax_grids
