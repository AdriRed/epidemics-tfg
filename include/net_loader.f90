module net_loader
   use fhash
   use iso_fortran_env, only: int32
   implicit none
   integer, parameter, private :: ik = int32
   type epidemic_net
      integer(ik), allocatable :: neighbour_counterpart_ptrs(:), neighbours(:), starter_ptrs(:), end_ptrs(:), degree(:)
      integer(ik) :: nodes_count = 0, links_count= 0
      type(int_int_hashmap) :: hashmap
   contains
      procedure, public :: get_neighbours_by_index
   end type epidemic_net
contains

   function get_neighbours_by_index(this, node_index) result(retval)
      integer(ik), allocatable:: retval(:)
      integer(ik), intent(in) :: node_index
      class(epidemic_net), intent(inout) :: this
      allocate(retval(this%degree(node_index)))
      retval = this%neighbours(this%starter_ptrs(node_index):this%end_ptrs(node_index))
   end function get_neighbours_by_index



   type(epidemic_net) function initialize_net(unit) result(net)
      integer(ik), intent(in) :: unit
      integer(ik) :: initial_links
      call init_hashmap(unit, net)
      write(*, *) 'Initialized hash map'
      write(*, *) '---- nodes count -> ', net%nodes_count
      write(*, *) '---- links count -> ', net%links_count
      initial_links = net%links_count
      allocate(net%neighbour_counterpart_ptrs(2*net%links_count), &
         net%neighbours(2*net%links_count), &
         net%starter_ptrs(net%nodes_count), &
         net%end_ptrs(net%nodes_count), &
         net%degree(net%nodes_count))

      call init_degrees_pointers(unit, net)
      write(*, *) 'Initialized degrees and pointers'

      call init_neighbours(unit, net)
      write(*, *) 'Initialized neighbour array'

      call clean_repeated_negibours(net)
      write(*, *) 'Cleaned neighbours. Reduced neighbours by ', net%links_count-initial_links
   end function initialize_net


   subroutine init_degrees_pointers(unit, net)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: iostat, node_a, node_b, index_node_a, index_node_b, i
      rewind(unit)

      do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit
         else
            if (node_a == node_b) cycle ! skip autolinks

            call net%hashmap%get(node_a, index_node_a)
            call net%hashmap%get(node_b, index_node_b)
            net%degree(index_node_a) = net%degree(index_node_a)+1
            net%degree(index_node_b) = net%degree(index_node_b)+1
         end if
      end do

      net%starter_ptrs(1) = 1
      net%end_ptrs(1) = 0
      do i = 1, net%nodes_count-1
         net%starter_ptrs(i+1) = net%starter_ptrs(i) + net%degree(i)
         net%end_ptrs(i+1) = net%starter_ptrs(i+1)-1
      end do

   end subroutine init_degrees_pointers

   integer(ik) function count_lines(unit) result(retval)
      integer(ik), intent(in) :: unit
      integer(ik) :: iostat
      retval = 0
      rewind(unit)
      iostat = 0

      do
         read(unit, *, iostat=iostat)
         if (iostat < 0) then
            exit
         else
            retval = retval + 1
         end if
      end do

      return
   end function count_lines

   subroutine init_hashmap(unit, net)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: i, iostat, dummy, node_a, node_b, lines
      logical :: exists
      i = 1
      lines = count_lines(unit)
      rewind(unit)
      iostat = 0
      call net%hashmap%reserve(lines) ! reserve N approx E nodes
      do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit
         else
            if (node_a == node_b) cycle ! skip autolinks

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
      end do

      net%nodes_count = net%hashmap%key_count()
      net%links_count = lines
   end subroutine init_hashmap

   subroutine init_neighbours(unit, net)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: iostat, node_a, node_b, index_node_a, index_node_b
      iostat = 0

      rewind(unit)
      do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat <0) exit

         if (node_a == node_b) cycle ! skip autolinks

         call net%hashmap%get(node_a, index_node_a)
         call net%hashmap%get(node_b, index_node_b)

         net%end_ptrs(index_node_a) = net%end_ptrs(index_node_a) + 1
         net%end_ptrs(index_node_b) = net%end_ptrs(index_node_b) + 1
         net%neighbours(net%end_ptrs(index_node_a)) = index_node_b
         net%neighbours(net%end_ptrs(index_node_b)) = index_node_a
         net%neighbour_counterpart_ptrs(net%end_ptrs(index_node_a)) = net%end_ptrs(index_node_b)
         net%neighbour_counterpart_ptrs(net%end_ptrs(index_node_b)) = net%end_ptrs(index_node_a)

      end do
   end subroutine init_neighbours

   subroutine clean_repeated_negibours(net)
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: i, j, k
      integer(ik) :: startp, endp, curr_neigh
      logical :: found


      do i = 1, net%nodes_count
         startp = net%starter_ptrs(i)
         endp   = net%end_ptrs(i)

         ! empty neighbours
         if (endp < startp) cycle

         j = startp
         do while (j <= endp)
            curr_neigh = net%neighbours(j)
            found = .false.

            ! search repeated neighbours
            do k = startp, j-1
               if (net%neighbours(k) == curr_neigh) then ! found repeated neighbour
                  found = .true.
                  exit
               end if
            end do

            if (found) then ! if repeated neighbour found
               ! replace with the latest neighbour
               net%neighbours(j) = net%neighbours(endp)
               ! reduce end pointer and degree
               endp = endp - 1
               net%degree(i) = net%degree(i) - 1
               net%links_count = net%links_count-1
               ! continue do while
               cycle
            else ! if not found, next neighbour
               j = j + 1
            end if
         end do

         net%end_ptrs(i) = endp
      end do
   end subroutine clean_repeated_negibours

end module net_loader
