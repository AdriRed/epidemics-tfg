module net_loader
   use fhash
   implicit none

   type epidemic_net
      integer, allocatable :: neighbours(:), starter_ptrs(:), end_ptrs(:), degree(:)
      integer :: nodes_count, links_count
      type(int_int_hashmap), private :: hashmap
   end type epidemic_net
contains
   type(epidemic_net) function initialize_net(unit) result(net)
      integer, intent(in) :: unit
      call init_hashmap(unit, net)
      call init_degrees_pointers(unit, net)
   end function initialize_net

   subroutine init_degrees_pointers(unit, net)
      integer, intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer :: iostat, node_a, node_b, index_node_a, index_node_b, i
      rewind(unit)

      loop: do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit loop
         else
            call net%hashmap%get(node_a, index_node_a)
            call net%hashmap%get(node_b, index_node_b)
            net%degree(index_node_a) = net%degree(index_node_a)+1
            net%degree(index_node_b) = net%degree(index_node_b)+1
         end if
      end do loop

      net%starter_ptrs(1) = 1
      net%end_ptrs(1) = 0
      do i = 1, net%nodes_count-1
         net%starter_ptrs(i+1) = net%starter_ptrs(i) + net%degree(i)
         net%end_ptrs(i+1) = net%starter_ptrs(i+1)-1
      end do

   end subroutine init_degrees_pointers

   subroutine init_hashmap(unit, net)
      integer, intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer :: i, iostat, dummy, node_a, node_b, lines
      logical :: exists
      rewind(unit)
      i = 1
      lines = 0
      call net%hashmap%reserve(100000000)
      loop: do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit loop
         else
            lines = lines+1
            call net%hashmap%get(node_a, dummy, exists)
            if (.not. exists) then
               call net%hashmap%set(node_a, i)
               i = i+1
            end if
            call net%hashmap%get(node_b, dummy, exists)
            if (.not. exists) then
               call net%hashmap%set(node_b, i)
               i = i+1
            end if

         end if
      end do loop

      net%nodes_count = net%hashmap%key_count()
      net%links_count = lines


      allocate(net%neighbours(net%links_count), &
         net%starter_ptrs(net%nodes_count), &
         net%end_ptrs(net%nodes_count), &
         net%degree(net%nodes_count))

   end subroutine init_hashmap
end module net_loader
