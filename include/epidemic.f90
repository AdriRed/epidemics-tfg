module epidemic
   use net_loader
   use mt19937_par
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter, private :: dp = real64
   integer, parameter, private :: ik = int32
   integer, parameter, private :: bk = int8

   integer(bk), parameter, public :: SIS_MODEL = 1
   integer(bk), parameter, public :: SIR_MODEL = 2

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
      real(dp) recovered_density
      type(epidemic_rates) rates
   end type epidemic_simulation_stats

   type epidemic_simulation
      type(epidemic_net) net
      type(epidemic_rates) actual_rates
      integer(bk) :: model_type ! 1 = SIS, 2 = SIR
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

      integer(ik) :: recovered_nodes_count
      type(mt19937_state) :: rnd
   contains
      procedure, public :: act
      procedure, public :: get_stats
      procedure, public :: infect_node
      procedure, public :: set_infected_node
      procedure, public :: recover_node
      procedure, public :: clear
      procedure, public :: calculate_actual_rates
      procedure, public :: verify_consistency
      procedure, private :: infect
      procedure, private :: recover
      procedure, private :: advance_time
      procedure, private :: remove_active_link
      procedure, private :: update_active_links_ptrs
      procedure, private :: add_active_link
   end type epidemic_simulation

contains

   type(epidemic_simulation) function initialize_simulation(net, seed, mode) result(retval)
      type(epidemic_net), intent(in) :: net
      integer(ik), intent(in) :: seed
      integer(bk), intent(in) :: mode
      retval%model_type = mode
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
      retval%recovered_nodes_count = 0
      write(*,*) 'Initialized simulation'
   end function initialize_simulation

   subroutine clear(this)
      class(epidemic_simulation), intent(inout) :: this
      deallocate(this%active_links, this%node_states, this%infected_nodes, this%neighbours_active_links_index)
   end subroutine clear

   type(epidemic_step_event) function act(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: numb
      call this%calculate_actual_rates()

      if (this%actual_rates%total_rate == 0.) then
         retval%selected_node = -1
         retval%action = 'E'
         return
      end if

      retval%elapsed_time = this%advance_time()
      numb = this%rnd%grnd()*this%actual_rates%total_rate
      if (numb < this%actual_rates%actual_infection_rate) then
         retval%selected_node = this%infect()
         retval%action = 'I'
      else
         retval%selected_node = this%recover()
         retval%action = 'R'
      end if
   end function act

   type(epidemic_simulation_stats) function get_stats(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      retval%rates = this%actual_rates
      retval%infected_density = real(this%infected_nodes_count, kind=dp) / this%net%stats%nodes_count
      if (this%model_type == SIS_MODEL) then
         retval%healthy_density = 1.-retval%infected_density
         retval%recovered_density = 0.
      else if (this%model_type == SIR_MODEL) then
         retval%recovered_density = real(this%recovered_nodes_count, kind=dp) / this%net%stats%nodes_count
         retval%healthy_density = 1.-retval%infected_density-retval%recovered_density
      end if
   end function get_stats

   subroutine calculate_actual_rates(this)
      class(epidemic_simulation), intent(inout) :: this
      this%actual_rates%actual_infection_rate = this%active_links_count*this%infection_rate
      this%actual_rates%actual_recovery_rate= this%infected_nodes_count*this%recovery_rate
      this%actual_rates%total_rate = this%actual_rates%actual_infection_rate + this%actual_rates%actual_recovery_rate
   end subroutine calculate_actual_rates

   real(dp) function advance_time(this) result(retval)
      class(epidemic_simulation), intent(inout) :: this
      real(dp) :: u

      ! Gillespie: tiempo al siguiente evento
      u = max(this%rnd%grnd(), epsilon(real(1.0, kind=dp)))
      retval = -log(u) / this%actual_rates%total_rate

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
      integer(ik) :: i, j, link_idx

      origin_node_idx = this%active_links(active_link_idx, 1)
      node_to_infect_idx = this%active_links(active_link_idx, 2)

      ! easy verification
      if (this%node_states(origin_node_idx) /= 1 .or. &
         this%node_states(node_to_infect_idx) /= 0) then
         write(*,*) "ERROR in infect_node: Invalid state combination"
         write(*,*) "  Origin state:", this%node_states(origin_node_idx)
         write(*,*) "  Target state:", this%node_states(node_to_infect_idx)
         stop
      end if

      ! delete the active selected link

      ! remove neighbour's active link pointer
      do i = this%net%starter_ptrs(origin_node_idx), this%net%end_ptrs(origin_node_idx)
         if (this%neighbours_active_links_index(i) == active_link_idx) then
            this%neighbours_active_links_index(i) = 0
            exit
         end if
      end do

      ! if not last active link in list, override this position with last active link
      if (active_link_idx /= this%active_links_count) then
         call this%update_active_links_ptrs(active_link_idx, this%active_links_count)
      end if

      call this%remove_active_link(active_link_idx)

      call this%set_infected_node(node_to_infect_idx, origin_node_idx)

   end subroutine infect_node

   subroutine set_infected_node(this, node_to_infect_idx, origin_node_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: node_to_infect_idx
      integer(ik), intent(inout), optional :: origin_node_idx
      integer(ik) :: neighbor, i, j, neighbor_link_idx
      logical :: present_value

      present_value = present(origin_node_idx)

      ! easy verification
      if (this%node_states(node_to_infect_idx) /= 0) then
         write(*, *) "ERROR in set_infected_node: Node already infected"
         write(*, *) "  Node:", node_to_infect_idx
         stop
      end if

      ! change status
      this%node_states(node_to_infect_idx) = 1

      ! add to infected nodes list
      this%infected_nodes_count = this%infected_nodes_count + 1
      this%infected_nodes(this%infected_nodes_count) = node_to_infect_idx

      ! foreach neigbour of the infected node
      do i = this%net%starter_ptrs(node_to_infect_idx), this%net%end_ptrs(node_to_infect_idx)
         neighbor = this%net%neighbours(i)

         ! if node was infected from other node skip original node
         if (present_value) then
            if(neighbor == origin_node_idx) cycle
         end if

         ! if neighbour of infected node is also infected
         if (this%node_states(neighbor) == 1) then

            ! search for any active link
            do j = this%net%starter_ptrs(neighbor), this%net%end_ptrs(neighbor)
               ! if current neighbour link is against the infected node
               if (this%net%neighbours(j) == node_to_infect_idx) then
                  neighbor_link_idx = this%neighbours_active_links_index(j)
                  ! if current neighbour link is valid
                  if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
                     ! if current neighbour link is not last update all pointers
                     if (neighbor_link_idx /= this%active_links_count) then
                        call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
                     end if
                     ! remove active link
                     call this%remove_active_link(neighbor_link_idx)
                     ! reset index
                     this%neighbours_active_links_index(j) = 0
                  end if
                  exit
               end if
            end do

            ! fast check for inverse links (should not exist, warning included)
            neighbor_link_idx = this%neighbours_active_links_index(i)
            if (neighbor_link_idx > 0 .and. neighbor_link_idx <= this%active_links_count) then
               write(*,*) "WARNING: Active link between two infected nodes found"
               if (neighbor_link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(neighbor_link_idx, this%active_links_count)
               end if
               call this%remove_active_link(neighbor_link_idx)
               this%neighbours_active_links_index(i) = 0
            end if

         else if (this%node_states(neighbor) == 0) then ! if neighbor is susceptible
            ! add to infection list
            neighbor_link_idx = this%add_active_link(node_to_infect_idx, neighbor)
            this%neighbours_active_links_index(i) = neighbor_link_idx
         end if
      end do

   end subroutine set_infected_node

   subroutine recover_node(this, recovered_idx)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik), intent(in) :: recovered_idx
      integer(ik) :: recovered_node, neighbor, i, j, link_idx

      recovered_node = this%infected_nodes(recovered_idx)

      ! easy check
      if (this%node_states(recovered_node) /= 1) then
         write(*, *) "ERROR in recover_node: Node not infected"
         write(*, *) "  Node:", recovered_node
         write(*, *) "  State:", this%node_states(recovered_node)
         stop
      end if

      ! change status
      if (this%model_type == SIS_MODEL) then
         this%node_states(recovered_node) = 0
      else if (this%model_type == SIR_MODEL) then
         this%node_states(recovered_node) = 2
         this%recovered_nodes_count = this%recovered_nodes_count + 1
      end if

      ! remove from infected list (swap with last)
      this%infected_nodes(recovered_idx) = this%infected_nodes(this%infected_nodes_count)
      this%infected_nodes_count = this%infected_nodes_count - 1

      ! remove active links from the recovered node
      do i = this%net%starter_ptrs(recovered_node), this%net%end_ptrs(recovered_node)
         link_idx = this%neighbours_active_links_index(i)
         ! if link is valid
         if (link_idx > 0 .and. link_idx <= this%active_links_count) then
            ! check consistency
            if (this%active_links(link_idx, 1) == recovered_node .and. &
               this%active_links(link_idx, 2) == this%net%neighbours(i)) then
               ! remove active link
               if (link_idx /= this%active_links_count) then
                  call this%update_active_links_ptrs(link_idx, this%active_links_count)
               end if
               call this%remove_active_link(link_idx)
            end if
            ! reset index
            this%neighbours_active_links_index(i) = 0
         end if
      end do

      ! add active link from infected nodes to recovered node if SIS
      if (this%model_type == SIS_MODEL) then
         do i = this%net%starter_ptrs(recovered_node), this%net%end_ptrs(recovered_node)
            neighbor = this%net%neighbours(i)
            ! if neighbour is infected
            if (this%node_states(neighbor) == 1) then
               ! search for link between neighbour and active link
               do j = this%net%starter_ptrs(neighbor), this%net%end_ptrs(neighbor)
                  if (this%net%neighbours(j) == recovered_node) then
                     ! add active link and exit loop
                     link_idx = this%add_active_link(neighbor, recovered_node)
                     this%neighbours_active_links_index(j) = link_idx
                     exit
                  end if
               end do
            end if
         end do
      end if

   end subroutine recover_node


   ! update nodes to new link
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


   subroutine verify_consistency(this)
      class(epidemic_simulation), intent(inout) :: this
      integer(ik) :: i, j, origin, target, ptr_idx, neighbour_idx
      integer(ik) :: count_infected_in_array, count_infected_in_states
      logical :: found

      ! ============================================
      ! 1. Verificar coherencia de infected_nodes
      ! ============================================

      ! Contar nodos en estado infectado
      count_infected_in_states = count(this%node_states == 1)

      if (count_infected_in_states /= this%infected_nodes_count) then
         write(*, *) 'ERROR: Infected count mismatch!'
         write(*, *) '  Nodes with state=1:', count_infected_in_states
         write(*, *) '  infected_nodes_count:', this%infected_nodes_count
         stop
      end if

      ! Verificar que todos los nodos en infected_nodes están realmente infectados
      do i = 1, this%infected_nodes_count
         if (this%infected_nodes(i) < 1 .or. this%infected_nodes(i) > this%net%stats%nodes_count) then
            write(*, *) 'ERROR: Invalid node index in infected_nodes!'
            write(*, *) '  Position:', i
            write(*, *) '  Node index:', this%infected_nodes(i)
            stop
         end if

         if (this%node_states(this%infected_nodes(i)) /= 1) then
            write(*, *) 'ERROR: Node in infected_nodes is not infected!'
            write(*, *) '  Position:', i
            write(*, *) '  Node index:', this%infected_nodes(i)
            write(*, *) '  State:', this%node_states(this%infected_nodes(i))
            stop
         end if
      end do

      ! Verificar que no hay duplicados en infected_nodes
      do i = 1, this%infected_nodes_count
         do j = i+1, this%infected_nodes_count
            if (this%infected_nodes(i) == this%infected_nodes(j)) then
               write(*, *) 'ERROR: Duplicate in infected_nodes!'
               write(*, *) '  Positions:', i, j
               write(*, *) '  Node:', this%infected_nodes(i)
               stop
            end if
         end do
      end do

      ! ============================================
      ! 2. Verificar coherencia de active_links
      ! ============================================

      do i = 1, this%active_links_count
         origin = this%active_links(i, 1)
         target = this%active_links(i, 2)

         ! Verificar índices válidos
         if (origin < 1 .or. origin > this%net%stats%nodes_count) then
            write(*, *) 'ERROR: Invalid origin node in active_link!'
            write(*, *) '  Link index:', i
            write(*, *) '  Origin:', origin
            stop
         end if

         if (target < 1 .or. target > this%net%stats%nodes_count) then
            write(*, *) 'ERROR: Invalid target node in active_link!'
            write(*, *) '  Link index:', i
            write(*, *) '  Target:', target
            stop
         end if

         ! Verificar que el origen está infectado
         if (this%node_states(origin) /= 1) then
            write(*, *) 'ERROR: Active link origin not infected!'
            write(*, *) '  Link index:', i
            write(*, *) '  Origin node:', origin
            write(*, *) '  Origin state:', this%node_states(origin)
            stop
         end if

         ! Verificar que el destino es susceptible
         if (this%node_states(target) /= 0) then
            write(*, *) 'ERROR: Active link target not susceptible!'
            write(*, *) '  Link index:', i
            write(*, *) '  Target node:', target
            write(*, *) '  Target state:', this%node_states(target)
            stop
         end if

         ! Verificar que origin y target son vecinos en la red
         found = .false.
         do j = this%net%starter_ptrs(origin), this%net%end_ptrs(origin)
            if (this%net%neighbours(j) == target) then
               found = .true.
               exit
            end if
         end do

         if (.not. found) then
            write(*, *) 'ERROR: Active link between non-neighbours!'
            write(*, *) '  Link index:', i
            write(*, *) '  Origin:', origin
            write(*, *) '  Target:', target
            stop
         end if
      end do

      ! Verificar que no hay duplicados en active_links
      do i = 1, this%active_links_count
         do j = i+1, this%active_links_count
            if (this%active_links(i, 1) == this%active_links(j, 1) .and. &
               this%active_links(i, 2) == this%active_links(j, 2)) then
               write(*, *) 'ERROR: Duplicate active link!'
               write(*, *) '  Positions:', i, j
               write(*, *) '  Origin:', this%active_links(i, 1)
               write(*, *) '  Target:', this%active_links(i, 2)
               stop
            end if
         end do
      end do

      ! ============================================
      ! 3. Verificar neighbours_active_links_index
      ! ============================================

      do i = 1, this%net%stats%nodes_count
         do j = this%net%starter_ptrs(i), this%net%end_ptrs(i)
            ptr_idx = this%neighbours_active_links_index(j)
            neighbour_idx = this%net%neighbours(j)

            if (ptr_idx > 0) then
               ! Si hay un puntero, debe ser válido
               if (ptr_idx > this%active_links_count) then
                  write(*, *) 'ERROR: neighbours_active_links_index out of bounds!'
                  write(*, *) '  Node:', i
                  write(*, *) '  Neighbour position:', j
                  write(*, *) '  Neighbour:', neighbour_idx
                  write(*, *) '  Pointer:', ptr_idx
                  write(*, *) '  active_links_count:', this%active_links_count
                  stop
               end if

               ! El enlace debe existir entre node i y su vecino
               if (this%active_links(ptr_idx, 1) /= i .or. &
                  this%active_links(ptr_idx, 2) /= neighbour_idx) then
                  write(*, *) 'ERROR: neighbours_active_links_index points to wrong link!'
                  write(*, *) '  Node:', i
                  write(*, *) '  Neighbour:', neighbour_idx
                  write(*, *) '  Pointer:', ptr_idx
                  write(*, *) '  Link origin:', this%active_links(ptr_idx, 1)
                  write(*, *) '  Link target:', this%active_links(ptr_idx, 2)
                  stop
               end if
            end if
         end do
      end do

      ! ============================================
      ! 4. Verificar que TODOS los enlaces activos tienen punteros
      ! ============================================

      do i = 1, this%active_links_count
         origin = this%active_links(i, 1)
         target = this%active_links(i, 2)
         found = .false.

         ! Buscar el puntero en neighbours_active_links_index
         do j = this%net%starter_ptrs(origin), this%net%end_ptrs(origin)
            if (this%net%neighbours(j) == target) then
               if (this%neighbours_active_links_index(j) == i) then
                  found = .true.
                  exit
               end if
            end if
         end do

         if (.not. found) then
            write(*, *) 'ERROR: Active link has no pointer in neighbours_active_links_index!'
            write(*, *) '  Link index:', i
            write(*, *) '  Origin:', origin
            write(*, *) '  Target:', target
            ! Mostrar qué punteros tiene el origen
            write(*, *) '  Pointers from origin:'
            do j = this%net%starter_ptrs(origin), this%net%end_ptrs(origin)
               write(*, *) '    Neighbour:', this%net%neighbours(j), &
                  ' Ptr:', this%neighbours_active_links_index(j)
            end do
            stop
         end if
      end do

      ! ============================================
      ! 5. Verificar coherencia de los rates
      ! ============================================

      call this%calculate_actual_rates()

      if (this%actual_rates%actual_infection_rate /= &
         this%active_links_count * this%infection_rate) then
         write(*, *) 'ERROR: actual_infection_rate inconsistent!'
         write(*, *) '  actual_infection_rate:', this%actual_rates%actual_infection_rate
         write(*, *) '  active_links_count * infection_rate:', &
            this%active_links_count * this%infection_rate
         stop
      end if

      if (this%actual_rates%actual_recovery_rate /= &
         this%infected_nodes_count * this%recovery_rate) then
         write(*, *) 'ERROR: actual_recovery_rate inconsistent!'
         write(*, *) '  actual_recovery_rate:', this%actual_rates%actual_recovery_rate
         write(*, *) '  infected_nodes_count * recovery_rate:', &
            this%infected_nodes_count * this%recovery_rate
         stop
      end if

      ! ============================================
      ! 6. Verificar que no hay enlaces "fantasma"
      ! ============================================

      ! Para cada nodo infectado, verificar que tiene los enlaces activos correctos
      do i = 1, this%infected_nodes_count
         origin = this%infected_nodes(i)

         ! Contar cuántos enlaces activos debería tener
         do j = this%net%starter_ptrs(origin), this%net%end_ptrs(origin)
            neighbour_idx = this%net%neighbours(j)

            if (this%node_states(neighbour_idx) == 0) then
               ! Debería haber un enlace activo
               ptr_idx = this%neighbours_active_links_index(j)

               if (ptr_idx == 0 .or. ptr_idx > this%active_links_count) then
                  write(*, *) 'ERROR: Missing active link for infected-susceptible pair!'
                  write(*, *) '  Infected node:', origin
                  write(*, *) '  Susceptible neighbour:', neighbour_idx
                  write(*, *) '  Pointer value:', ptr_idx
                  stop
               end if
            elseif (this%node_states(neighbour_idx) == 1) then
               ! NO debería haber enlace activo
               ptr_idx = this%neighbours_active_links_index(j)

               if (ptr_idx > 0 .and. ptr_idx <= this%active_links_count) then
                  write(*, *) 'ERROR: Active link between two infected nodes!'
                  write(*, *) '  Infected node 1:', origin
                  write(*, *) '  Infected node 2:', neighbour_idx
                  write(*, *) '  Pointer:', ptr_idx
                  stop
               end if
            elseif (this%node_states(neighbour_idx) == 2) then
               ! En SIR: NUNCA debe haber enlace activo hacia recuperados
               if (this%model_type == SIR_MODEL) then
                  ptr_idx = this%neighbours_active_links_index(j)
                  if (ptr_idx > 0 .and. ptr_idx <= this%active_links_count) then
                     write(*, *) 'ERROR: Active link to recovered node in SIR!'
                     write(*, *) '  Infected node:', origin
                     write(*, *) '  Recovered node:', neighbour_idx
                     write(*, *) '  Pointer:', ptr_idx
                     stop
                  end if
               end if
            end if
         end do
      end do

      ! ============================================
      ! 7. Estadísticas de debugging
      ! ============================================

      write(*, *) '--- Consistency Check Passed ---'
      write(*, *) '  Time:', this%time
      write(*, *) '  Infected nodes:', this%infected_nodes_count
      write(*, *) '  Active links:', this%active_links_count
      write(*, *) '  Infection rate:', this%actual_rates%actual_infection_rate
      write(*, *) '  Recovery rate:', this%actual_rates%actual_recovery_rate

   end subroutine verify_consistency

end module epidemic
