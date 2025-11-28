program main
   use net_loader
   implicit none

   type(epidemic_net) net

   open(unit=11, file='../ignore-files/large_twitch_edges.csv', action='read')

   net = initialize_net(11)
   write(*, *) 'Memory used (bytes): ', net_memory_bytes(net)

   close(unit=11)
end program main
