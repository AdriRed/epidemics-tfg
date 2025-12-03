module epidemic
   use net_loader
   use mt19937
   use iso_fortran_env, only: int32, real64
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ik = int32

   type epidemic_simulation
      type(epidemic_net) net
      real(dp) :: infection_rate = 0.0, recovery_rate= 0.0
      real(dp) :: time = 0.0
      integer(ik), allocatable :: infected_node_indexes(:)
      integer(ik), allocatable :: potential_infection_link_indexes(:)

   contains
      procedure, public :: advance_time
      procedure, private :: get_rates
   end type epidemic_simulation

contains

   real(dp) function get_rates(sim) result(retval)
      class(epidemic_simulation), intent(inout) :: sim
      retval = sim%infection_rate + sim%recovery_rate
   end function get_rates

   real(dp) function advance_time(sim) result(retval)
      class(epidemic_simulation), intent(inout) :: sim
      real(dp) :: tau, test, distr, rate
      logical :: found = .false.

      rate = sim%get_rates()
      do while (.not. found)
         tau = grnd()
         test = grnd()
         distr = rate*exp(-rate*tau)
         found = test .ge. distr
      end do
      sim%time = sim%time + tau
      retval = tau
   end function advance_time

end module epidemic
