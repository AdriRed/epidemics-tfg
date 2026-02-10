module net_loader
   use fhash
   use iso_fortran_env, only: int32, real64
   implicit none
   integer, parameter, private :: ik = int32
   integer, parameter, private :: dp = real64
   

   type epidemic_net_stats
      integer(ik) :: nodes_count, links_count
      real(dp) :: average_degree, average_degree_sqrd
   end type epidemic_net_stats
   type epidemic_net
      integer(ik), allocatable :: neighbour_counterpart_ptrs(:), neighbours(:), starter_ptrs(:), end_ptrs(:), degree(:)
      type(epidemic_net_stats) :: stats
      type(int_int_hashmap) :: hashmap
   contains
      procedure, public :: get_neighbours_by_index
      procedure, public :: print_stats
   end type epidemic_net
contains

   function get_neighbours_by_index(this, node_index) result(retval)
      integer(ik), allocatable:: retval(:)
      integer(ik), intent(in) :: node_index
      class(epidemic_net), intent(inout) :: this
      allocate(retval(this%degree(node_index)))
      retval = this%neighbours(this%starter_ptrs(node_index):this%end_ptrs(node_index))
   end function get_neighbours_by_index

   subroutine print_stats(this)
      class(epidemic_net), intent(inout) :: this
      write(*, *) '--- stats ---'
      write(*, "(A, I8)") 'N     = ', this%stats%nodes_count
      write(*, "(A, I8)") 'E     = ', this%stats%links_count
      write(*, "(A, F17.8)") '<k>   = ', this%stats%average_degree
      write(*, "(A, F17.8)") '<k^2> = ', this%stats%average_degree_sqrd
      write(*, "(A, F17.8)") 'Var k = ', sqrt(this%stats%average_degree_sqrd - &
         this%stats%average_degree*this%stats%average_degree)

   end subroutine print_stats
   


   type(epidemic_net) function initialize_net(unit) result(net)
      integer(ik), intent(in) :: unit
      integer(ik) :: initial_links
      call init_hashmap(unit, net)
      write(*, *) 'Initialized hash map'
      initial_links = net%stats%links_count
      allocate(net%neighbour_counterpart_ptrs(2*net%stats%links_count), &
         net%neighbours(2*net%stats%links_count), &
         net%starter_ptrs(net%stats%nodes_count), &
         net%end_ptrs(net%stats%nodes_count), &
         net%degree(net%stats%nodes_count))

      net%neighbour_counterpart_ptrs(:) = 0
      net%neighbours(:) = 0
      net%starter_ptrs(:) = 0
      net%end_ptrs(:) = 0
      net%degree(:) = 0
      

      call init_degrees_pointers(unit, net)
      write(*, *) 'Initialized degrees and pointers'

      call init_neighbours(unit, net)
      write(*, *) 'Initialized neighbour array'

      call clean_repeated_neighbours(net)
      write(*, *) 'Cleaned neighbours. Reduced neighbours by ', net%stats%links_count-initial_links

      call calculate_stats(net)


   end function initialize_net

   subroutine calculate_stats(net)
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: i
      net%stats%average_degree = 0.
      net%stats%average_degree_sqrd = 0.

      do i = 1, net%stats%nodes_count
         net%stats%average_degree = net%stats%average_degree + net%degree(i)
         net%stats%average_degree_sqrd = net%stats%average_degree_sqrd + net%degree(i)*net%degree(i)
      end do

      net%stats%average_degree = net%stats%average_degree / net%stats%nodes_count
      net%stats%average_degree_sqrd = net%stats%average_degree_sqrd / net%stats%nodes_count

   end subroutine calculate_stats


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
      do i = 1, net%stats%nodes_count-1
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

      net%stats%nodes_count = net%hashmap%key_count()
      net%stats%links_count = lines
   end subroutine init_hashmap

   subroutine init_neighbours(unit, net)
      integer(ik), intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: iostat, node_a, node_b, index_node_a, index_node_b
      iostat = 0

      rewind(unit)
      do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) exit

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

   subroutine clean_repeated_neighbours(net)
      type(epidemic_net), intent(inout) :: net
      integer(ik) :: i, j, k
      integer(ik) :: startp, endp, curr_neigh
      logical :: found


      do i = 1, net%stats%nodes_count
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
               net%stats%links_count = net%stats%links_count-1
               ! continue do while
               cycle
            else ! if not found, next neighbour
               j = j + 1
            end if
         end do

         net%end_ptrs(i) = endp
      end do
   end subroutine clean_repeated_neighbours

end module net_loader
