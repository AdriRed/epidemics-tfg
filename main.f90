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
   type(epidemic_stats) :: stats
   integer :: i

   call init_genrand(42069)

   ! open(unit=11, file='./ignore-files/musae_git_edges.csv', action='read')
   ! open(unit=11, file='./files/test-file-1.txt', action='read')
   open(unit=11, file='./ignore-files/large_twitch_edges.csv', action='read')
   net = initialize_net(11)
   call net%hashmap%clear()
   close(unit=11)

   call execute_simulation(net, 0.5_dp, 2.71_dp, int(1E7, kind=ik), 10000)


contains

   subroutine execute_simulation(initialized_net, infection_rate, recovery_rate, limit_steps, output_file_steps)
      class(epidemic_net), intent(inout) :: initialized_net
      real(dp), intent(in) :: infection_rate, recovery_rate
      integer(ik), intent(in) :: limit_steps, output_file_steps
      character(70) :: name
      integer(ik) :: percentage_steps
      character(len=:), allocatable :: filename

      percentage_steps = int(real(limit_steps, kind=dp)/1.E2, kind=ik)

      write(name, '(A,F12.9,A,F12.9)') 'I=', infection_rate, '-R=', recovery_rate
      filename = trim(adjustl(name))
      write(*, *) name
      simulation = initialize_simulation(initialized_net)
      simulation%infection_rate = infection_rate
      simulation%recovery_rate = recovery_rate

      open(unit=12, file='./output/events-'// filename //'.dat', action='write')
      open(unit=13, file='./output/stats-'// filename //'.dat', action='write')
      write(12, "(E20.10, A2, I10)") 0._dp, 'I', 1
      call simulation%set_infected_node(1)
      do i = 1, limit_steps
         event = simulation%act()
         if (mod(i,percentage_steps) == 0) write(*, "(I3.3, A)") i/percentage_steps, '%'
         if (mod(i, output_file_steps) == 0) then
            stats = simulation%get_stats()
            write(12, "(E20.10, A2, I10)") simulation%time, event%action, event%selected_node
            write(13, "(E20.10, E20.10, E20.10, E20.10)") simulation%time, &
               stats%rates%actual_infection_rate, stats%rates%actual_recovery_rate, &
               stats%infected_density
         end if
         if (simulation%infected_nodes_count == simulation%net%nodes_count) then
            write(*, *) '100% infection reached'
            exit
         elseif (simulation%infected_nodes_count == 0) then
            write(*, *) '100% health reached'
            exit
         end if
      end do

      if (mod(i, 100) /= 0) then
         write(12, "(E20.10, A2, I10)") simulation%time, event%action, event%selected_node
         write(13, "(E20.10, E20.10, E20.10, E20.10)") simulation%time, &
            stats%rates%actual_infection_rate, stats%rates%actual_recovery_rate, stats%infected_density
      end if
      close(unit=12)
      close(unit=13)
      call simulation%clear()
   end subroutine execute_simulation

end program main
