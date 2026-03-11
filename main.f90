program main
   use net_loader
   use epidemic
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ik = int32
   integer, parameter :: bk = int8

   type(epidemic_net) :: net
   integer :: i_lambdas
   character(:), allocatable :: name
   
   ! open(unit=11, file='./nets/out.moreno_beach_beach', action='read')
   open(unit=11, file='./nets/musae_git_edges.csv', action='read')

   net = initialize_net(11, weighted=.false.)
   call net%hashmap%clear()
   call net%print_stats()
   close(unit=11)


   ! call sis_prevalence(net)
   name = 'musae_git_edges'
   call execute_simulation(net, 2._dp, 1._dp, 42072, 500._dp, SIR_MODEL, stats_unit=21, events_unit=22, net_name=name)

   ! !$omp parallel do private(i_lambdas) schedule(dynamic)
   ! do i_lambdas = 100, 1, -1
   !    call execute_simulation(net, i_lambdas/real(100, dp), 1._dp, 1._dp, 42069, real(2, dp), stats_unit=(i_lambdas+50))
   ! end do
   ! !$omp end parallel do



contains

   subroutine sir_model_evolution(init_net)
      implicit none
      type(epidemic_net), intent(inout) :: init_net
      integer(ik) :: i_sim, j_sim
      !$omp parallel do private(i_sim, j_sim) schedule(dynamic) collapse(2)
      do i_sim = 1, 100
         do j_sim = 1, 60
            call execute_simulation(init_net, real(i_sim*2, dp)/100, real(1., dp) &
               , 42069+j_sim, real(100, dp), SIR_MODEL, stats_unit=(i_sim*100+j_sim))
         end do
      end do
      !$omp end parallel do
   end subroutine sir_model_evolution

   subroutine sis_prevalence(init_net)
      implicit none
      type(epidemic_net), intent(inout) :: init_net
      type(epidemic_simulation) :: simulation
      type(epidemic_step_event) :: event
      type(epidemic_simulation_stats) :: stats
      real(dp) :: time_limit, eval_time, relax_time
      integer(ik) :: i_sim, node_id
      simulation = initialize_simulation(init_net, 42070, SIS_MODEL)
      ! set all nodes infected

      do i_sim = 1, init_net%stats%nodes_count
         call simulation%set_infected_node(i_sim)
      end do


      eval_time = 0.
      time_limit = 600.
      relax_time = 50.
      ! delta
      simulation%recovery_rate = 1
      simulation%infection_rate = 1
      open(unit=22, file='density_by_rate_weighted.dat', action='write')
      ! open(unit=23, file='density_by_rate_weighted_events.dat', action='write')
      do i_sim = 100, 1, -1
         ! lambda
         simulation%infection_rate = real(i_sim, kind=dp)/100.
         write(*, *) "Starting for infection rate = ", simulation%infection_rate
         do while (eval_time < time_limit)

            event = simulation%act()
            eval_time = eval_time + event%elapsed_time
            stats = simulation%get_stats()
            ! call simulation%verify_consistency()
            if (relax_time < eval_time) then
               write(22, *) simulation%time, simulation%infection_rate, &
                  stats%infected_density
            end if
            write(*, "(A, F5.3, A, F15.5, A, F5.3)") "Infection rate = ", simulation%infection_rate, ", Time = ", simulation%time, &
               ", Density = ", stats%infected_density
            call init_net%rev_hashmap%get(event%selected_node, node_id)
            ! write(23, "(F15.5, F15.5, A1, I10)") simulation%time, simulation%infection_rate, event%action, node_id
            if (stats%infected_density == 0.) then
               write(*, *) "Reached healthy state"
               exit
            end if
         end do
         eval_time = eval_time - time_limit
         if (stats%infected_density == 0.) then
            exit
         end if
      end do
      close(22)
      call simulation%clear()
      ! close(23)
   end subroutine sis_prevalence

   subroutine execute_simulation(initialized_net, infection_rate, recovery_rate, &
      seed, limit_time, model_type, stats_unit, events_unit, net_name)
      implicit none
      class(epidemic_net), intent(inout) :: initialized_net
      real(dp), intent(in) :: infection_rate, recovery_rate, limit_time
      integer(ik), intent(in) :: seed
      integer(ik), intent(in), optional :: stats_unit, events_unit
      integer(bk), intent(in) :: model_type
      character(:), intent(in), allocatable, optional :: net_name
      type(epidemic_simulation) :: sim
      type(epidemic_step_event) :: sim_event
      type(epidemic_simulation_stats) :: sim_stats
      logical :: should_write_stats, should_write_events

      character(70) :: name1, name2, name3
      integer(ik) :: i_sim, node_id, max_degree, max_degree_node_index
      character(len=:), allocatable :: filename
      filename = ''

      write(name1, '(A,F10.5,A,F10.5,A,I5)') 'I=', infection_rate, '-R=', recovery_rate, '-S=', seed

      if (model_type == SIR_MODEL) then
         write(name2, '(A, A)') 'SIR-', trim(name1)
      else
         write(name2, '(A, A)') 'SIS-', trim(name1)
      end if

      if (net%weighted) then
         write(name3, '(A, A)') 'w', trim(name2)
      else
         write(name3, '(A)') trim(name2)
      end if

      if (present(net_name)) then
         filename= trim(net_name) // '-' // trim(name3)
      else
         filename= trim(name3)
      end if
      !$omp critical(name_write)
      write(*, '(A, A)') 'Filename will be ', filename
      !$omp end critical(name_write)

      sim = initialize_simulation(initialized_net, seed, model_type)

      max_degree = 0
      do i_sim = 1, initialized_net%stats%nodes_count
         if (initialized_net%degree(i_sim) > max_degree) then
            max_degree = initialized_net%degree(i_sim)
            max_degree_node_index = i_sim
         end if
      end do

      should_write_stats = present(stats_unit)
      should_write_events = present(events_unit)

      sim%infection_rate = infection_rate
      sim%recovery_rate = recovery_rate
      i_sim = 0
      call sim%set_infected_node(max_degree_node_index)
      call initialized_net%rev_hashmap%get(1, node_id)

      !$omp critical(file_write)
      if (should_write_stats) then
         open(unit=stats_unit, file='./output/stats-'// filename //'.dat', action='write')
         write(stats_unit, *) '# Epidemic simulation - stats'
         if (initialized_net%weighted) then
            write(stats_unit, *) '# Weighted: YES'

         else
            write(stats_unit, *) '# Weighted: NO'
         end if
         if (model_type == SIR_MODEL) then
            write(stats_unit, *) '# Model: SIR'
         else
            write(stats_unit, *) '# Model: SIS'
         end if

         if (present(net_name)) then
            write(stats_unit, *) '# Net name: ', net_name
         end if
         write(stats_unit, *) '# Using mt19993 random generator with seed: ', seed
         write(stats_unit, *) '# Infection rate: ', infection_rate
         write(stats_unit, *) '# Recovery rate: ', recovery_rate
         write(stats_unit, *) '# Start node: ', node_id
         write(stats_unit, *) '# Start node degree: ', max_degree
         write(stats_unit, *) '# -----------------------------------------------------'
         write(stats_unit, *) '# time, infected_density, recovered_density, actual_infection_rate, actual_recovery_rate'

      end if
      if (should_write_events) then
         open(unit=events_unit, file='./output/events-'// filename //'.dat', action='write')

         write(events_unit, *) '# Epidemic simulation - events'
         if (initialized_net%weighted) then
            write(events_unit, *) '# Weighted: YES'

         else
            write(events_unit, *) '# Weighted: NO'
         end if
         if (model_type == SIR_MODEL) then
            write(events_unit, *) '# Model: SIR'
         else
            write(events_unit, *) '# Model: SIS'
         end if

         if (present(net_name)) then
            write(events_unit, *) '# Net name: ', net_name
         end if
         write(events_unit, *) '# Using mt19993 random generator with seed: ', seed
         write(events_unit, *) '# Infection rate: ', infection_rate
         write(events_unit, *) '# Recovery rate: ', recovery_rate
         write(events_unit, *) '# Start node: ', node_id
         write(events_unit, *) '# Start node degree: ', max_degree
         write(events_unit, *) '# -----------------------------------------------------'


         write(events_unit, *) '# time, node_id, event'
      end if
      !$omp end critical(file_write)
      write(*, '(A,F10.5,A,I5,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=start'
      do
         sim_event = sim%act()

         if (should_write_stats) then
            sim_stats = sim%get_stats()

            !$omp critical(file_write)
            write(stats_unit, "(E20.10, E20.10, E20.10, E20.10, E20.10)") sim%time, &
               sim_stats%infected_density, sim_stats%recovered_density, &
               sim_stats%rates%actual_infection_rate, sim_stats%rates%actual_recovery_rate
            !$omp end critical(file_write)
         end if

         if (should_write_events) then
            call initialized_net%rev_hashmap%get(sim_event%selected_node, node_id)
            !$omp critical(file_write)
            write(events_unit, "(E20.10, I10, A2)") sim%time, node_id, sim_event%action
            !$omp end critical(file_write)
         end if

         if (sim_event%action == 'E') then
            write(*, '(A,F10.5,A,I5,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=dead'
            exit
         end if
         ! write(*, '(A,F10.5,A,I5,A,F10.5)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=',sim%time

         if (sim%time > limit_time) then
            write(*, '(A,F10.5,A,I5,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=max'
            exit
         end if

      end do
      !$omp critical(file_write)
      if (should_write_stats) close(unit=stats_unit)
      if (should_write_events) close(unit=events_unit)
      !$omp end critical(file_write)
      call sim%clear()
   end subroutine execute_simulation
end program main
