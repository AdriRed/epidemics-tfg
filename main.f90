program main
   use net_loader
   use epidemic
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter :: dp = real64
   integer, parameter :: ik = int32
   integer, parameter :: bk = int8

   type(epidemic_net) :: net

   character(:), allocatable :: option
   character(len=256) :: buf
   real(dp), allocatable :: rates(:,:)
   integer(ik) :: i, j, rate_size

   rate_size = 10000

   allocate(rates(rate_size, 3))
   do i = 1, 100
      do j = 1, 100
         rates((i-1)*100+j, 1) = 0.01*i
         rates((i-1)*100+j, 2) = 1
         rates((i-1)*100+j, 3) = j
      end do
   end do

   open(unit=11, file='./ignore-files/musae_git_edges.csv', action='read')
   ! open(unit=11, file='./ignore-files/ia-infect-dublin.mtx', action='read')
   ! open(unit=11, file='./files/test-file-1.txt', action='read')
   ! open(unit=11, file='./ignore-files/large_twitch_edges.csv', action='read')
   ! open(unit=11, file='./ignore-files/soc-epinions.mtx', action='read')
   net = initialize_net(11)
   call net%hashmap%clear()
   call net%print_stats()
   close(unit=11)

   ! option = ' '

   ! do
   !    write(*, *) 'Write stats in the following format: "{infection rate} {recovery rate}" or "E" to end'
   !    read(*,'(A)') buf             ! lee la línea completa en option
   !    option = trim(adjustl(buf))
   !    if (option == 'E' .or. option == 'e') exit
   !    write(*, "(A, A)") 'Using option: ', option
   !    ! internal read desde la cadena option para obtener los números
   !    read(option, *) inf_rate, rec_rate
   !    call execute_simulation(net, inf_rate, rec_rate, int(1E6, kind=ik), 100)
   ! end do

   !$omp parallel do private(i) schedule(dynamic)
   do i = 1, rate_size
      call execute_simulation(net, rates(i, 1), rates(i, 2), int(1E6, kind=ik), 100, 10+i, &
         int(rates(i, 3), kind=ik), real(1E-10, kind=dp))
   end do
   !$omp end parallel do
contains

   subroutine execute_simulation(initialized_net, infection_rate, recovery_rate, &
      limit_steps, output_file_steps, unit, seed, epsilon)
      implicit none
      class(epidemic_net), intent(in) :: initialized_net
      real(dp), intent(in) :: infection_rate, recovery_rate, epsilon
      integer(ik), intent(in) :: limit_steps, output_file_steps, unit, seed
      real(dp) :: last_density

      character(70) :: name
      integer(ik) :: percentage_steps, i_step
      character(len=:), allocatable :: filename
      type(epidemic_simulation) :: simulation
      type(epidemic_step_event) :: event
      type(epidemic_simulation_stats) :: stats

      percentage_steps = int(real(limit_steps, kind=dp)/1.E2, kind=ik)

      write(name, '(A,F10.5,A,F10.5,A,I5)') 'I=', infection_rate, '-R=', recovery_rate, '-S=', seed
      filename = trim(adjustl(name))
      ! $omp critical(name_write)
      write(*, *) name
      ! $omp end critical(name_write)


      simulation = initialize_simulation(initialized_net, seed)
      simulation%infection_rate = infection_rate
      simulation%recovery_rate = recovery_rate

      ! $omp critical(file_write)
      open(unit=unit, file='./output/stats-'// filename //'.dat', action='write')
      ! $omp end critical(file_write)
      last_density = -1.
      call simulation%set_infected_node(1)
      do i_step = 0, limit_steps
         event = simulation%act()

         if (mod(i_step,percentage_steps) == 0) then
            ! $omp critical(output)
            write(*, "(A, F4.2, A, F4.2, A, I5, A, I3.3, A)") &
               'I=',infection_rate, ', R=',recovery_rate, ', S=', seed, ' - ', i_step/percentage_steps, '%'
            ! $omp end critical(output)

         end if
         ! write(*, *) 'Last density difference = ', last_density-stats%infected_density
         if (output_file_steps == 1 .or. mod(i_step, output_file_steps) == 0) then
            stats = simulation%get_stats()
            !$omp critical(file_write)
            write(unit, "(E20.10, E20.10, E20.10, E20.10)") simulation%time, &
               stats%rates%actual_infection_rate, stats%rates%actual_recovery_rate, &
               stats%infected_density
            flush(unit)
            !$omp end critical(file_write)

            if (abs(last_density-stats%infected_density) < epsilon) then
               write(*, *) 'Stationary state reached, density: ',  stats%infected_density
               exit
            end if

            last_density = stats%infected_density

         end if
         if (simulation%infected_nodes_count == simulation%net%stats%nodes_count) then
            ! $omp critical(output)
            write(*, *) '100% infection reached'
            ! $omp end critical(output)
            exit
         elseif (simulation%infected_nodes_count == 0) then
            ! $omp critical(output)
            write(*, *) '100% health reached'
            ! $omp end critical(output)
            exit
         end if
      end do
      
      stats = simulation%get_stats()

      ! if (mod(i, 100) /= 0) then
      ! $omp critical(file_write)
      write(unit, "(E20.10, E20.10, E20.10, E20.10)") simulation%time, &
         stats%rates%actual_infection_rate, stats%rates%actual_recovery_rate, stats%infected_density
      ! $omp end critical(file_write)

      ! end if
      close(unit=unit)
      call simulation%clear()
   end subroutine execute_simulation
end program main
