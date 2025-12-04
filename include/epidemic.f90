module epidemic
   use net_loader
   use mt19937
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter, private :: dp = real64
   integer, parameter, private :: ik = int32
   integer, parameter, private :: bk = int8

   type epidemic_rates
      real(dp) actual_infection_rate
      real(dp) actual_recovery_rate
      real(dp) total_rate
   end type epidemic_rates

   type epidemic_simulation
      type(epidemic_net) net
      type(epidemic_rates) actual_rates

      real(dp) :: infection_rate = 0.0, recovery_rate= 0.0
      real(dp) :: time = 0.0

      integer(ik), allocatable :: infected_nodes(:)
      integer(ik) :: infected_nodes_count = 0

      integer(ik), allocatable :: active_links(:, :)
      integer(ik) :: active_links_count = 0

      integer(bk), allocatable :: node_states(:)

      integer(ik), allocatable :: neighbours_active_links_index(:)
   contains
      procedure, public :: advance_time
      procedure, private :: calculate_actual_rates
      procedure, public :: act
      procedure, private :: infect
      procedure, private :: recover
   end type epidemic_simulation

contains

   type(epidemic_simulation) function initialize_simulation(net) result(retval)
      type(epidemic_net), intent(in) :: net

      retval%net = net
      allocate(retval%active_links(net%links_count, 2), &
         retval%node_states(net%nodes_count), &
         retval%infected_nodes(net%nodes_count), &
         retval%neighbours_active_links_index(2*net%links_count))

   end function initialize_simulation

   subroutine infect(this)
      class(epidemic_simulation), intent(inout) :: this
   end subroutine infect

   subroutine recover(this)
      class(epidemic_simulation), intent(inout) :: this

   end subroutine recover

   subroutine act(this)
      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: numb
      numb = grnd()*this%actual_rates%total_rate

      if (numb < this%actual_rates%actual_infection_rate) then
         call this%infect()
      else
         call this%recover()
      end if
   end subroutine act

   subroutine calculate_actual_rates(this)
      class(epidemic_simulation), intent(inout) :: this
      this%actual_rates%actual_infection_rate = this%active_links_count*this%infection_rate
      this%actual_rates%actual_recovery_rate= this%infected_nodes_count*this%recovery_rate
      this%actual_rates%total_rate = this%actual_rates%actual_infection_rate + this%actual_rates%actual_recovery_rate
   end subroutine calculate_actual_rates

   real(dp) function advance_time(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: tau, test, distr
      logical :: found = .false.
      call this%calculate_actual_rates()
      do while (.not. found)
         tau = grnd()
         test = grnd()
         distr = this%actual_rates%total_rate*exp(-this%actual_rates%total_rate*tau)
         found = test .ge. distr
      end do
      this%time = this%time + tau
      retval = tau
   end function advance_time

end module epidemic
