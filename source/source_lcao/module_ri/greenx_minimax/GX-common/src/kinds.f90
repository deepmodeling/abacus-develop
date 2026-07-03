!> Precision constants
module kinds
   implicit none

   !> Single precision
   integer, parameter, public :: sp = selected_real_kind(6, 30)
   !> Double precision
   integer, parameter, public :: dp = selected_real_kind(14, 200)
   !> Medium length character
   integer, parameter, public :: medium_char = 100
   !> Long length character
   integer, parameter, public :: long_char = 200

end module kinds
