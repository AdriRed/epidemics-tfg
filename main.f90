program main
   use net_loader
   implicit none

   type(epidemic_net) net

   open(unit=11, file='../files/test-file-1.txt', action='read')

   net = initialize_net(11)

   write(*, *) 'NEIGHBOUR LENGTHS', size(net%neighbours)

   close(unit=11)
end program main
