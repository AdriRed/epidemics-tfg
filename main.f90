program main
   use net_loader
   use epidemic
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ik = int32
   integer, parameter :: bk = int8

   type(epidemic_net) :: net
   type(epidemic_simulation) :: simulation
   type(epidemic_step_event) :: event
   type(epidemic_simulation_stats) :: stats
   integer(ik) :: i, j
   real(dp) :: time_limit, eval_time, relax_time
   open(unit=11, file='./ignore-files/musae_git_edges.csv', action='read')

   net = initialize_net(11)
   call net%hashmap%clear()
   call net%print_stats()
   close(unit=11)

   ! simulation = initialize_simulation(net, 42069, SIS_MODEL)
   ! set all nodes infected
   ! do i = 1, net%stats%nodes_count
   !    call simulation%set_infected_node(i)
   ! end do
   ! eval_time = 0.
   ! time_limit = 2.
   ! relax_time = 0.5
   ! ! delta
   ! simulation%recovery_rate = 1
   ! simulation%infection_rate = 1
   ! open(unit=12, file='density_by_rate.dat', action='write')
   ! do i = 100, 1, -1
   !    ! lambda
   !    simulation%infection_rate = real(i, kind=dp)/100.
   !    write(*, *) "Starting for infection rate = ", simulation%infection_rate
   !    do while (eval_time < time_limit)

   !       event = simulation%act()
   !       eval_time = eval_time + event%elapsed_time
   !       stats = simulation%get_stats()
   !       ! call simulation%verify_consistency()
   !       if (relax_time < eval_time) then
   !          write(12, *) simulation%time, simulation%infection_rate, &
   !          stats%infected_density
   !       end if
         
   !       write(*, "(A, F5.3, A, F10.5, A, F5.3)") "Infection rate = ", simulation%infection_rate, ", Time = ", simulation%time, &
   !          ", Density = ", stats%infected_density
   !    end do
   !    eval_time = eval_time - time_limit

   ! end do
   ! close(12)



   ! !$omp parallel do private(i_lambdas) schedule(dynamic)
   ! do i_lambdas = 100, 1, -1
   !    call execute_simulation(net, i_lambdas/real(100, dp), real(1, dp),real(1, dp), i_lambdas+50, 42069, real(2, dp))
   ! end do
   ! !$omp end parallel do
   ! close(12)

   call execute_simulation(net, real(0.3, dp), real(1., dp), 21, 2, real(1000, dp), SIR_MODEL)

contains

   subroutine execute_simulation(initialized_net, infection_rate, recovery_rate, &
      unit, seed, limit_time, model_type)
      implicit none
      class(epidemic_net), intent(in) :: initialized_net
      real(dp), intent(in) :: infection_rate, recovery_rate, limit_time
      integer(ik), intent(in) :: unit, seed
      integer(bk), intent(in) :: model_type
      type(epidemic_simulation) :: sim
      type(epidemic_step_event) :: sim_event
      type(epidemic_simulation_stats) :: sim_stats

      character(70) :: name
      integer(ik) :: i_sim
      character(len=:), allocatable :: filename

      if (model_type == SIR_MODEL) then
         write(name, '(A,F10.5,A,F10.5,A,I5)') 'SIR-I=', infection_rate, '-R=', recovery_rate, '-S=', seed
      else if (model_type == SIS_MODEL) then
         write(name, '(A, A,F10.5,A,F10.5,A,I5)') 'SIS-I=', infection_rate, '-R=', recovery_rate, '-S=', seed
      end if
      
         filename = trim(adjustl(name))
      !$omp critical(name_write)
      write(*, *) name
      !$omp end critical(name_write)

      sim = initialize_simulation(initialized_net, seed, model_type)
      ! set all nodes infected
      ! do i_sim = 1, initialized_net%stats%nodes_count
      !    call sim%set_infected_node(i_sim)
      ! end do
      call sim%set_infected_node(1)

      sim%infection_rate = infection_rate
      sim%recovery_rate = recovery_rate

      !$omp critical(file_write)
      open(unit=unit, file='./output/stats-'// filename //'.dat', action='write')
      !$omp end critical(file_write)
      do
         sim_event = sim%act()

         if (sim_event%action == 'E') then
            write(*, '(A,F10.5,A,I5,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=dead'
            exit
         end if
         write(*, '(A,F10.5,A,I5,A,F10.5)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=',sim%time
         sim_stats = sim%get_stats()
         !$omp critical(file_write)
         write(unit, "(E20.10, E20.10, E20.10, E20.10, E20.10)") sim%time, &
            sim_stats%rates%actual_infection_rate, sim_stats%rates%actual_recovery_rate, &
            sim_stats%infected_density, sim_stats%recovered_density
         !$omp end critical(file_write)

         if (sim%time > limit_time) then
            write(*, '(A,F10.5,A,I5,A)') 'I/R=', infection_rate/recovery_rate, '-S=', seed, '-t=max'
            exit
         end if
      end do


      close(unit=unit)
      call sim%clear()
   end subroutine execute_simulation
end program main
