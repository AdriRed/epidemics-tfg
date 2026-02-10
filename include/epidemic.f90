module epidemic
   use net_loader
   use mt19937_par
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

   type epidemic_step_event
      integer(ik) selected_node
      character(1) action
      real(dp) elapsed_time
   end type epidemic_step_event

   type epidemic_simulation_stats
      real(dp) infected_density
      real(dp) healthy_density
      type(epidemic_rates) rates
   end type epidemic_simulation_stats

   type epidemic_simulation
      type(epidemic_net) net
      type(epidemic_rates) actual_rates
      
      real(dp) :: infection_rate = 0.0, recovery_rate= 0.0
      real(dp) :: time = 0.0
      ! lista de nodos infectados
      integer(ik), allocatable :: infected_nodes(:)
      integer(ik) :: infected_nodes_count = 0
      ! lista de links activos (linx entre un nodo infectado y otro no infectado)
      integer(ik), allocatable :: active_links(:, :)
      integer(ik) :: active_links_count = 0
      ! estados de nodos
      integer(bk), allocatable :: node_states(:)

      ! lista de indices de active_links a los cuales hace referencia el mismo elemento de vecinos
      ! active_links                  = [4-3][1-2]
      ! neighbours                    = [2][3]|[1]|[1][4]|[3][5]...
      ! neighbours_active_links_index = [2][ ]|[ ]|[ ][ ]|[1][ ]...
      ! posición coincidente con neighbours y referencia a índice de active_links
      integer(ik), allocatable :: neighbours_active_links_index(:)

      type(mt19937_state) :: rnd
   contains
      procedure, public :: act
      procedure, public :: get_stats
      procedure, public :: infect_node
      procedure, public :: set_infected_node
      procedure, public :: recover_node
      procedure, public :: clear
      procedure, public :: calculate_actual_rates
      procedure, private :: infect
      procedure, private :: recover
      procedure, private :: advance_time
      procedure, private :: remove_active_link
      procedure, private :: update_active_links_ptrs
      procedure, private :: add_active_link
   end type epidemic_simulation

contains

   type(epidemic_simulation) function initialize_simulation(net, seed) result(retval)
      type(epidemic_net), intent(in) :: net
      integer(ik), intent(in) :: seed

      retval%net = net
      allocate(retval%active_links(2*net%stats%links_count, 2), &
         retval%node_states(net%stats%nodes_count), &
         retval%infected_nodes(net%stats%nodes_count), &
         retval%neighbours_active_links_index(2*net%stats%links_count))

      call retval%rnd%init_genrand(seed)
      retval%active_links(:, :) = 0
      retval%node_states(:) = 0
      retval%infected_nodes(:) = 0
      retval%neighbours_active_links_index(:) = 0
      retval%time = 0
      write(*,*) 'Initialized simulation'
   end function initialize_simulation

   subroutine clear(this)
      class(epidemic_simulation), intent(inout) :: this
      deallocate(this%active_links, this%node_states, this%infected_nodes, this%neighbours_active_links_index)
   end subroutine clear

   type(epidemic_step_event) function act(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: numb
      retval%elapsed_time = this%advance_time()
      numb = this%rnd%grnd()*this%actual_rates%total_rate
      write(*, *) "Choosen random ", numb
      write(*, *) "Rates INF = ", this%actual_rates%actual_infection_rate
      write(*, *) "Rates REC = ", this%actual_rates%actual_recovery_rate
      if (numb < this%actual_rates%actual_infection_rate) then
         write(*, *) "INFECTION"
         retval%selected_node = this%infect()
         retval%action = 'I'
      else
         write(*, *) "RECOVERY"
         retval%selected_node = this%recover()
         retval%action = 'R'
      end if
   end function act

   type(epidemic_simulation_stats) function get_stats(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      retval%rates = this%actual_rates
      retval%infected_density = real(this%infected_nodes_count) / this%net%stats%nodes_count
      retval%healthy_density = 1-retval%infected_density
   end function get_stats

   subroutine calculate_actual_rates(this)
      class(epidemic_simulation), intent(inout) :: this
      this%actual_rates%actual_infection_rate = this%active_links_count*this%infection_rate
      this%actual_rates%actual_recovery_rate= this%infected_nodes_count*this%recovery_rate
      this%actual_rates%total_rate = this%actual_rates%actual_infection_rate + this%actual_rates%actual_recovery_rate
   end subroutine calculate_actual_rates

   real(dp) function advance_time(this) result(retval)
      ! class(epidemic_simulation), intent(inout) :: this
      ! real(dp) :: tau, test, distr
      ! logical :: found
      ! found = .false.
      ! call this%calculate_actual_rates()
      ! do while (.not. found)
      !    tau = this%rnd%grnd()
      !    test = this%rnd%grnd()
      !    distr = this%actual_rates%total_rate*exp(-this%actual_rates%total_rate*tau)
      !    found = test .le. distr ! gillespielle will mark as valid when test <= distr
      ! end do
      ! this%time = this%time + tau
      ! retval = tau

      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: u, R

      call this%calculate_actual_rates()
      R = this%actual_rates%total_rate

      ! Gillespie: tiempo al siguiente evento
      u = max(this%rnd%grnd(), epsilon(real(1.0, kind=dp)))
      retval = -log(u) / R

      this%time = this%time + retval
   end function advance_time

   ! infects a random node
   ! returns: the selected node index
   integer(ik) function infect(this) result(chosen_node)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: chosen_link_idx
      chosen_link_idx = int(this%rnd%grnd()*this%active_links_count, kind=ik)+1
      chosen_node = this%active_links(chosen_link_idx, 2)

      call this%infect_node(chosen_link_idx)
   end function infect

   integer(ik) function recover(this) result(chosen_node)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: chosen_node_idx
      chosen_node_idx = int(this%rnd%grnd()*this%infected_nodes_count, kind=ik)+1
      chosen_node = this%infected_nodes(chosen_node_idx)
      call this%recover_node(chosen_node_idx)
   end function recover

   ! infects a targeted node given by the first parameter
   subroutine infect_node(this, active_link_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: active_link_idx
      integer(ik) :: origin_node_idx, node_to_infect_idx
      ! selection of the active links
      origin_node_idx = this%active_links(active_link_idx, 1) ! the infected node
      node_to_infect_idx = this%active_links(active_link_idx, 2) ! the node to infect
      
      ! if selected active_index is not the last
      if (active_link_idx /= this%active_links_count) then
         call this%update_active_links_ptrs(active_link_idx, this%active_links_count)
      end if
      ! swap with last active link
      call this%remove_active_link(active_link_idx)
      ! add to infected nodes
      call this%set_infected_node(node_to_infect_idx, origin_node_idx)

   end subroutine infect_node

   subroutine set_infected_node(this, node_to_infect_idx, origin_node_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: node_to_infect_idx
      integer(ik), intent(inout), optional :: origin_node_idx
      integer(ik) :: neighbor, i, neighbor_link_idx

      logical :: present_value
      
      present_value = present(origin_node_idx)
      ! change node state
      if (this%node_states(node_to_infect_idx) == 1) then
         write(*, *) "WARNING - Infecting already infected node"
      end if
      this%node_states(node_to_infect_idx) = 1
      
      this%infected_nodes_count = this%infected_nodes_count+1
      this%infected_nodes(this%infected_nodes_count) = node_to_infect_idx
      ! foreach neighbour
      do i = this%net%starter_ptrs(node_to_infect_idx), this%net%end_ptrs(node_to_infect_idx)
         neighbor = this%net%neighbours(i)
         
         if (present_value) then
            if (neighbor == origin_node_idx) cycle
         end if

         if (this%node_states(neighbor) == 1) then
            ! remove potential active links between node to infect and already infected node
            neighbor_link_idx = this%neighbours_active_links_index(i)
            if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
               ! compute last BEFORE removal
               if (neighbor_link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
               end if
               call this%remove_active_link(neighbor_link_idx)
               ! clear the mapping for this neighbour position
               this%neighbours_active_links_index(i) = 0
            end if

            ! counterpart (other direction)
            neighbor_link_idx = this%neighbours_active_links_index(this%net%neighbour_counterpart_ptrs(i))
            if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
               if (neighbor_link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
               end if
               call this%remove_active_link(neighbor_link_idx)
               this%neighbours_active_links_index(this%net%neighbour_counterpart_ptrs(i)) = 0
            end if
         elseif (this%node_states(neighbor) == 0) then
            ! add link between node to infect and node not infected
            neighbor_link_idx = this%add_active_link(node_to_infect_idx, neighbor)
            this%neighbours_active_links_index(i) = neighbor_link_idx
         end if
      end do

   end subroutine set_infected_node

   subroutine recover_node(this, recovered_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: recovered_idx
      integer(ik) :: neighbor_link_idx, i, recovered_node, neighbor
      recovered_node = this%infected_nodes(recovered_idx)
      if (this%node_states(recovered_node) == 0) then
         write(*, *) "WARNING - Recovering already recovered node"
      end if
      ! change node state
      this%node_states(recovered_node) = 0
      ! remove infected node
      this%infected_nodes(recovered_idx) = this%infected_nodes(this%infected_nodes_count)
      this%infected_nodes_count = this%infected_nodes_count - 1

      ! foreach neighbour
      do i = this%net%starter_ptrs(recovered_node), this%net%end_ptrs(recovered_node)
         neighbor = this%net%neighbours(i)
         if (this%node_states(neighbor) == 1) then
            ! add as active link with the self recovered node
            neighbor_link_idx = this%add_active_link(neighbor, recovered_node)
            this%neighbours_active_links_index(this%net%neighbour_counterpart_ptrs(i)) = neighbor_link_idx
         elseif (this%node_states(neighbor) == 0) then
            ! remove active link between two recovered/susceptible nodes
            neighbor_link_idx = this%neighbours_active_links_index(i)
            if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
               if (neighbor_link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
               end if
               call this%remove_active_link(neighbor_link_idx)
               this%neighbours_active_links_index(i) = 0
            end if

            ! counterpart
            neighbor_link_idx = this%neighbours_active_links_index(this%net%neighbour_counterpart_ptrs(i))
            if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
               if (neighbor_link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
               end if
               call this%remove_active_link(neighbor_link_idx)
               this%neighbours_active_links_index(this%net%neighbour_counterpart_ptrs(i)) = 0
            end if
         end if
      end do
   end subroutine recover_node

   
   ! hace swap de dos 
   subroutine update_active_links_ptrs(this, new_idx, old_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: new_idx, old_idx
      integer(ik) :: origin_node, i

      origin_node = this%active_links(old_idx, 1)

      ! para cada uno de los vecinos de un nodo
      do i = this%net%starter_ptrs(origin_node), this%net%end_ptrs(origin_node)
         ! si el vecino apunta al index antiguo
         if (this%neighbours_active_links_index(i) == old_idx) then
            ! haz que apunte al nuevo index
            this%neighbours_active_links_index(i) = new_idx   
         end if
      end do

   end subroutine update_active_links_ptrs

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
      this%active_links(idx, :) = this%active_links(this%active_links_count, :)
      this%active_links_count = this%active_links_count - 1
   end subroutine remove_active_link




end module epidemic
