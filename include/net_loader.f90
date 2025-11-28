module net_loader
   use fhash
   use iso_fortran_env, only: int64
   implicit none

   type epidemic_net
      integer, allocatable :: neighbours(:), starter_ptrs(:), end_ptrs(:), degree(:)
      integer :: nodes_count, links_count
      type(int_int_hashmap), private :: hashmap
   contains
      procedure, public :: get_neighbours


   end type epidemic_net
contains

   function get_neighbours(this, node_id) result(retval)
      integer, allocatable:: retval(:)
      integer, intent(in) :: node_id
      class(epidemic_net), intent(inout) :: this
      integer :: index
      call this%hashmap%get(node_id, index)
      allocate(retval(this%end_ptrs(index)-this%starter_ptrs(index)))
      retval = this%neighbours(this%starter_ptrs(index):this%end_ptrs(index))
   end function get_neighbours

   type(epidemic_net) function initialize_net(unit) result(net)
      integer, intent(in) :: unit
      call init_hashmap(unit, net)
      write(*, *) 'Initialized hash map'
      write(*, *) '---- nodes count -> ', net%nodes_count
      write(*, *) '---- links count -> ', net%links_count

      allocate(net%neighbours(2*net%links_count), &
         net%starter_ptrs(net%nodes_count), &
         net%end_ptrs(net%nodes_count), &
         net%degree(net%nodes_count))

      call init_degrees_pointers(unit, net)
      write(*, *) 'Initialized degrees and pointers'

      call init_neighbours(unit, net)
      write(*, *) 'Initialized neighbour array'

   end function initialize_net

   integer(int64) function net_memory_bytes(net) result(total)
      type(epidemic_net), intent(in) :: net
      integer(int64) :: bytes_element

      total = 0

      ! neighbours
      bytes_element = storage_size(net%neighbours) / 8
      total = total + size(net%neighbours, kind=int64) * bytes_element

      ! starter_ptrs
      bytes_element = storage_size(net%starter_ptrs) / 8
      total = total + size(net%starter_ptrs, kind=int64) * bytes_element

      ! end_ptrs
      bytes_element = storage_size(net%end_ptrs) / 8
      total = total + size(net%end_ptrs, kind=int64) * bytes_element

      ! degree
      bytes_element = storage_size(net%degree) / 8
      total = total + size(net%degree, kind=int64) * bytes_element

      ! Hashmap → depende de fhash (tienes que mirar su implementación)
      ! Si tiene allocatables internos, habrá que sumarlos igual que arriba

   end function net_memory_bytes


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

   integer function count_lines(unit) result(retval)
      integer, intent(in) :: unit
      integer :: iostat
      retval = 0
      rewind(unit)

      loop: do
         read(unit, *, iostat=iostat)
         if (iostat < 0) then
            exit loop
         else
            retval = retval + 1
         end if
      end do loop

      return
   end function count_lines

   subroutine init_hashmap(unit, net)
      integer, intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer :: i, iostat, dummy, node_a, node_b, lines
      logical :: exists
      i = 1
      lines = count_lines(unit)
      rewind(unit)
      call net%hashmap%reserve(lines) ! reserve N approx E nodes
      loop: do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit loop
         else
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
   end subroutine init_hashmap

   subroutine init_neighbours(unit, net)
      integer, intent(in) :: unit
      type(epidemic_net), intent(inout) :: net
      integer :: iostat, node_a, node_b, index_node_a, index_node_b

      rewind(unit)
      loop: do
         read(unit, *, iostat=iostat) node_a, node_b
         if (iostat < 0) then
            exit loop
         else
            call net%hashmap%get(node_a, index_node_a)
            call net%hashmap%get(node_b, index_node_b)

            net%end_ptrs(index_node_a) = net%end_ptrs(index_node_a) + 1
            net%end_ptrs(index_node_b) = net%end_ptrs(index_node_b) + 1
            net%neighbours(net%end_ptrs(index_node_a)) = index_node_b
            net%neighbours(net%end_ptrs(index_node_b)) = index_node_a

         end if
      end do loop
   end subroutine init_neighbours
end module net_loader
