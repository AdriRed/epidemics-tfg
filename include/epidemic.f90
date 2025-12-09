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
      procedure, private :: remove_active_link
      procedure, private :: add_active_link
      procedure, private :: calculate_actual_rates
      procedure, public :: act
      procedure, private :: infect
      procedure, public :: infect_node
      procedure, private :: recover
      procedure, public :: recover_node
   end type epidemic_simulation

contains

   type(epidemic_simulation) function initialize_simulation(net) result(retval)
      type(epidemic_net), intent(in) :: net

      retval%net = net
      allocate(retval%active_links(net%links_count, 2), &
         retval%node_states(net%nodes_count), &
         retval%infected_nodes(net%nodes_count), &
         retval%neighbours_active_links_index(2*net%links_count))

      retval%active_links(:, :) = 0
      retval%node_states(:) = 0
      retval%infected_nodes(:) = 0
      retval%neighbours_active_links_index(:) = 0

   end function initialize_simulation

   subroutine infect(this)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: chosen_link_idx
      chosen_link_idx = int(grnd()*this%active_links_count, kind=ik)+1
      call this%infect_node(chosen_link_idx)
   end subroutine infect

   subroutine infect_node(this, active_link_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: active_link_idx
      integer(ik) :: neighbor_link_idx
      integer(ik) :: origin_node_idx, node_to_infect_idx
      integer(ik) :: i
      origin_node_idx = this%active_links(active_link_idx, 1)
      node_to_infect_idx = this%active_links(active_link_idx, 2)
      ! change node state
      this%node_states(node_to_infect_idx) = 1
      ! swap with last active link
      call this%remove_active_link(active_link_idx)
      ! add to infected nodes
      this%infected_nodes_count = this%infected_nodes_count+1
      this%infected_nodes(this%infected_nodes_count) = node_to_infect_idx

      ! foreach neighbour
      do i = this%net%starter_ptrs(node_to_infect_idx), this%net%end_ptrs(node_to_infect_idx)
         if (this%node_states(i) == 1) then
            ! remove potential active links between node to infect and already infected node
            neighbor_link_idx = this%neighbours_active_links_index(i)
            if (this%active_links(neighbor_link_idx, 2) == node_to_infect_idx) then
               call this%remove_active_link(neighbor_link_idx)
            end if
         elseif (this%node_states(i) == 0) then
            ! add link between node to infect and node not infected
            neighbor_link_idx = this%add_active_link(node_to_infect_idx, i)
            this%neighbours_active_links_index(i) = neighbor_link_idx
         end if
      end do

   end subroutine infect_node

   integer(ik) function add_active_link(this, origin_node, target_node) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: origin_node, target_node

      this%active_links_count = this%active_links_count+1
      this%active_links(this%active_links_count, 1) = origin_node
      this%active_links(this%active_links_count, 2) = target_node

      retval = this%active_links_count
   end function add_active_link

   subroutine remove_active_link(this, idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: idx
      this%active_links(idx, 1) = this%active_links(this%active_links_count, 1)
      this%active_links(idx, 2) = this%active_links(this%active_links_count, 2)
      this%active_links_count = this%active_links_count - 1
   end subroutine remove_active_link

   subroutine recover(this)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: chosen_node_idx
      chosen_node_idx = int(grnd()*this%infected_nodes_count, kind=ik)+1
      call this%recover_node(chosen_node_idx)
   end subroutine recover

   subroutine recover_node(this, recovered_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: recovered_idx
      integer(ik) :: neighbor_link_idx, i, recovered_node
      recovered_node = this%infected_nodes(recovered_idx)
      ! change node state
      this%node_states(recovered_node) = 0
      ! remove infected node
      this%infected_nodes(recovered_idx) = this%infected_nodes(this%infected_nodes_count)
      this%infected_nodes_count = this%infected_nodes_count - 1

      ! foreach neighbour
      do i = this%net%starter_ptrs(recovered_node), this%net%end_ptrs(recovered_node)
         if (this%node_states(i) == 1) then
            ! add as active link with the self recovered node
            neighbor_link_idx = this%add_active_link(recovered_node, i)
            this%neighbours_active_links_index(i) = neighbor_link_idx
         elseif (this%node_states(i) == 0) then
            ! remove active link between two recovered nodes
            neighbor_link_idx = this%neighbours_active_links_index(i)
            call this%remove_active_link(neighbor_link_idx)
         end if
      end do
   end subroutine recover_node

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
