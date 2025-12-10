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

   ! open(unit=11, file='./ignore-files/large_twitch_edges.csv', action='read')
   open(unit=11, file='./files/test-file-1.txt', action='read')
   net = initialize_net(11)
   call net%hashmap%clear()
   close(unit=11)

   simulation = initialize_simulation(net)
   simulation%infection_rate = 0.1
   simulation%recovery_rate = 0.07

   write(12, *) 'I', 1
   call simulation%set_infected_node(1)
   open(unit=12, file='./builds/events.dat', action='write')
   open(unit=13, file='./builds/stats.dat', action='write')
   do i = 1, 100000
      event = simulation%act()
      stats = simulation%get_stats()
      write(12, *) simulation%time, event%action, event%selected_node
      write(13, *) simulation%time, stats%rates%actual_infection_rate, stats%healthy_density
   end do

   close(unit=12)
   close(unit=13)
end program main
