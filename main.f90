program main
   use net_loader
   use epidemic
   implicit none

   type(epidemic_net) net

   open(unit=11, file='./files/test-file-2.txt', action='read')

   net = initialize_net(11)
   call net%hashmap%clear()
   write(*, *) "1 -> ", net%get_neighbours_by_index(1)
   write(*, *) "2 -> ", net%get_neighbours_by_index(2)
   write(*, *) "3 -> ", net%get_neighbours_by_index(3)
   write(*, *) "4 -> ", net%get_neighbours_by_index(4)
   write(*, *) "5 -> ", net%get_neighbours_by_index(5)
   write(*, *) "6 -> ", net%get_neighbours_by_index(6)
   write(*, *) "7 -> ", net%get_neighbours_by_index(7)
   write(*, *) "8 -> ", net%get_neighbours_by_index(8)
   write(*, *) "9 -> ", net%get_neighbours_by_index(9)
   write(*, *) "10 -> ", net%get_neighbours_by_index(10)
   write(*, *) "11 -> ", net%get_neighbours_by_index(11)


   close(unit=11)
end program main
