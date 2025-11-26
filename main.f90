program main
   use net_loader
   implicit none

   type(epidemic_net) net

   open(unit=11, file='../files/test-file-1.txt', action='read')

   net = initialize_net(11)

   write(*, *) 'NEIGHBOUR LENGTHS', size(net%neighbours)

   write(*, *) 1, "->", net%get_neighbours(1)
   write(*, *) 2, "->", net%get_neighbours(2)
   write(*, *) 3, "->", net%get_neighbours(3)
   write(*, *) 4, "->", net%get_neighbours(4)
   write(*, *) 5, "->", net%get_neighbours(5)
   write(*, *) 6, "->", net%get_neighbours(6)
   write(*, *) 7, "->", net%get_neighbours(7)
   write(*, *) 8, "->", net%get_neighbours(8)
   write(*, *) 9, "->", net%get_neighbours(9)
   write(*, *) 10, "->", net%get_neighbours(10)
   write(*, *) 11, "->", net%get_neighbours(11)

   close(unit=11)
end program main
