module reversed_skiplist
   use mt19937_par
   use iso_fortran_env, only: int8, int32, real64
   implicit none

   integer, parameter, private :: dp = real64
   integer, parameter, private :: ik = int32
   integer, parameter, private :: bk = int8


   type skiplist_entry_ptr
      type(skiplist_entry), pointer :: ptr => null()
   end type skiplist_entry_ptr

   type skiplist_entry
      real(dp) :: weight, total_weight
      integer(ik), allocatable :: indexes(:)
      integer(ik) :: indexes_count
      type(skiplist_entry_ptr), allocatable :: next(:)
   end type skiplist_entry

   type skiplist
      type(mt19937_state) :: rng
      integer(ik) :: max_level
      type(skiplist_entry_ptr), pointer :: head
      real(dp) :: promote_probability
      integer(ik) :: max_indexes
   contains
      procedure :: clear => skiplist_clear
      procedure :: add => skiplist_add
      procedure :: debug_print => skiplist_print
      procedure :: remove_entry => skiplist_remove_entry
      procedure :: update_index => skiplist_update_index
   end type skiplist
contains

   subroutine skiplist_print(this)
      class(skiplist), intent(inout) :: this
      type(skiplist_entry), pointer :: previous
      integer(ik) :: i
      print*, '-------------SKIPLIST-------------'
      do i = 1, this%max_level
         previous => this%head%ptr
         write(*, "(A11, I1, A5)") "---- Level ", i, " ----"
         do while(associated(previous))
            call print_node(previous)
            previous => previous%next(i)%ptr
         end do
         print *, 'NULL'
      end do
   end subroutine

   subroutine print_node(previous)
      type(skiplist_entry), pointer :: previous
      integer(ik) :: j
      write(*, "(F8.1, A2)", advance="no") previous%weight, ' ('
      if (previous%indexes_count > 0) then
         write(*, "(I3)", advance="no") previous%indexes(1)
         do j = 2, previous%indexes_count
            write(*, "(A1, I3)", advance="no") ',', previous%indexes(j)
         end do
      end if
      write(*, "(A5)", advance="no") ') ->'
   end subroutine print_node

   function init_skiplist(maxlevels, seed, probability, max_indexes) result(retval)
      integer(ik), intent(in) :: maxlevels, max_indexes, seed
      real(dp), intent(in) :: probability
      type(skiplist) :: retval
      call retval%rng%init_genrand(seed)
      retval%max_level = maxlevels
      retval%promote_probability = probability
      retval%max_indexes = max_indexes
      allocate(retval%head)

   end function init_skiplist

   subroutine skiplist_update_index(this, weight, position, new_value)
      class(skiplist), intent(inout) :: this
      real(dp), intent(in) :: weight
      integer(ik), intent(in) :: position
      integer(ik), intent(in) :: new_value
      type(skiplist_entry), pointer :: current
      integer(ik) :: level

      ! Buscar el nodo con el peso dado (recorrido simple por nivel 1)
      current => this%head%ptr
      level = 1
      do
         if (.not. associated(current%next(level)%ptr)) then ! end of the line

            level = level+1
            if (level > this%max_level) exit ! final level reached
            cycle ! if not final level reached keep searching
         else
            if (weight < current%next(level)%ptr%weight) then ! if next lesser
               current => current%next(level)%ptr ! get next in line
            else if (weight == current%next(level)%ptr%weight ) then ! insertion in node only
               current => current%next(level)%ptr
               exit
            else ! if next isnt lesser, go to level below
               level = level+1
               if (level > this%max_level) then
                  exit
               end if
            end if
         end if
      end do

      current%indexes(position) = new_value

   end subroutine

   subroutine skiplist_add(this, weight, index, index_position)
      class(skiplist), intent(inout) :: this
      real(dp), intent(in) :: weight
      integer(ik), intent(in) :: index
      integer(ik), intent(inout), optional :: index_position
      integer(ik) :: level, i
      real(dp) :: r
      type(skiplist_entry), pointer :: data, previous, old_head
      type(skiplist_entry_ptr), allocatable :: previous_ptrs(:)

      logical :: insert_head_link
      allocate(data)
      data%weight = weight

      allocate(data%indexes(this%max_indexes), data%next(this%max_level))
      data%indexes_count = 1
      data%indexes(1) = index
      data%total_weight = weight
      if (present(index_position)) index_position = 1
      if (.not. associated(this%head%ptr)) then ! case where list is empty
         this%head%ptr => data
      else if (this%head%ptr%weight <= weight) then ! case where is going to be first
         if (this%head%ptr%weight == weight) then !! if same weight
            deallocate(data%indexes, data%next)
            deallocate(data)
            this%head%ptr%indexes_count = this%head%ptr%indexes_count +1
            this%head%ptr%indexes(this%head%ptr%indexes_count) = index
            this%head%ptr%total_weight = this%head%ptr%total_weight + weight
            if (present(index_position)) then

               index_position = this%head%ptr%indexes_count
            end if
         else !insert new head
            insert_head_link = .true.
            old_head => this%head%ptr
            this%head%ptr => data

            do i = this%max_level, 1, -1 ! linking new head the same as old head
               if (insert_head_link) then
                  this%head%ptr%next(i)%ptr => old_head
                  r = this%rng%grnd()
                  insert_head_link = r < this%promote_probability
               else
                  this%head%ptr%next(i)%ptr => old_head%next(i)%ptr
                  old_head%next(i)%ptr => null()
               end if
            end do
         end if
      else ! normal case
         level = 1
         previous => this%head%ptr
         allocate(previous_ptrs(this%max_level))

         do
            if (.not. associated(previous%next(level)%ptr)) then ! end of the line

               previous_ptrs(level)%ptr => previous
               level = level+1
               if (level > this%max_level) exit ! final level reached
               cycle ! if not final level reached keep searching
            else
               if (weight < previous%next(level)%ptr%weight) then ! if next lesser
                  previous => previous%next(level)%ptr ! get next in line
               else if (weight == previous%next(level)%ptr%weight ) then ! insertion in node only
                  previous => previous%next(level)%ptr
                  exit
               else ! if next isnt lesser, go to level below
                  previous_ptrs(level)%ptr => previous
                  level = level+1
                  if (level > this%max_level) then
                     exit
                  end if
               end if
            end if
         end do

         if (previous%weight == weight) then
            deallocate(data%indexes, data%next)
            deallocate(data)
            previous%indexes_count = previous%indexes_count +1
            previous%indexes(previous%indexes_count) = index
            previous%total_weight = previous%total_weight + weight
            if (present(index_position)) then
               index_position = previous%indexes_count
            end if
         else !insert new node

            ! promoting
            level = this%max_level
            do
               data%next(level)%ptr => previous_ptrs(level)%ptr%next(level)%ptr ! new node pointing now to old next node
               previous_ptrs(level)%ptr%next(level)%ptr => data
               r = this%rng%grnd()
               if (r < this%promote_probability) then
                  level = level - 1
                  if (level == 0) exit
               else
                  exit
               end if
            end do
         end if
         deallocate(previous_ptrs)


      end if

   end subroutine skiplist_add

   subroutine skiplist_clear(this)
      class(skiplist), intent(inout) :: this
      type(skiplist_entry), pointer :: current
      type(skiplist_entry), pointer :: next_ptr

      if (.not. associated(this%head)) return

      current => this%head%ptr

      do while (associated(current))
         next_ptr => current%next(this%max_level)%ptr
         deallocate(current%indexes)
         deallocate(current%next)
         deallocate(current)
         current => next_ptr
      end do

      nullify(this%head%ptr)
      deallocate(this%head)

      this%max_level = 0
   end subroutine skiplist_clear

   ! subroutine skiplist_remove(this, weight, index_value)

   ! end subroutine skiplist_remove

   subroutine skiplist_remove_entry(this, weight, index_position, update_value, node_deletion)
      class(skiplist), intent(inout) :: this
      real(dp), intent(in) :: weight
      integer(ik), intent(in) :: index_position
      integer(ik), intent(out), optional :: update_value
      logical, intent(out), optional :: node_deletion
      type(skiplist_entry_ptr), allocatable :: previous_ptrs(:)
      type(skiplist_entry), pointer :: before_ptr, current, delete_node
      integer(ik) :: curr_lvl, i
      if (.not. associated(this%head%ptr)) return
      curr_lvl = 1
      before_ptr => this%head%ptr
      allocate(previous_ptrs(this%max_level))
      delete_node => before_ptr
      if (delete_node%weight == weight) then ! case where head is going to be deleted

         if (delete_node%indexes_count /= 1_ik) then ! not last index to remove
            current => this%head%ptr
            ! No es el último índice -> solo reestructura interna
            if (present(update_value)) update_value = current%indexes(current%indexes_count)
            current%indexes(index_position) = current%indexes(current%indexes_count)
            current%indexes(current%indexes_count) = 0
            current%indexes_count = current%indexes_count - 1
            current%total_weight = current%total_weight - current%weight
            node_deletion = .false.

         else
            this%head%ptr => this%head%ptr%next(this%max_level)%ptr ! pointer points to next node
            node_deletion = .true.
            do i = 1, this%max_level
               if (.not. associated(this%head%ptr, delete_node%next(i)%ptr)) then
                  this%head%ptr%next(i)%ptr => delete_node%next(i)%ptr
               end if
            end do

            deallocate(delete_node%next, delete_node%indexes)
            deallocate(delete_node)
         end if

         deallocate(previous_ptrs)
         return
      end if

      do
         if (.not. associated(before_ptr%next(curr_lvl)%ptr)) then ! end of the line
            previous_ptrs(curr_lvl)%ptr => before_ptr
            curr_lvl = curr_lvl+1
            if (curr_lvl > this%max_level) then
               exit ! final level reached
            end if
            cycle ! if not final level reached keep searching
         else
            if (weight < before_ptr%next(curr_lvl)%ptr%weight) then ! if next lesser
               before_ptr => before_ptr%next(curr_lvl)%ptr ! get next in line
               delete_node => before_ptr

            else if (weight == before_ptr%next(curr_lvl)%ptr%weight ) then ! if next is weight to remove
               delete_node => before_ptr%next(curr_lvl)%ptr
               if (delete_node%indexes_count /= 1_ik) then ! not last index to remove
                  current => before_ptr%next(curr_lvl)%ptr
                  ! No es el último índice -> solo reestructura interna
                  if (present(update_value)) update_value = current%indexes(current%indexes_count)
                  current%indexes(index_position) = current%indexes(current%indexes_count)
                  current%indexes(current%indexes_count) = 0
                  current%indexes_count = current%indexes_count - 1
                  current%total_weight = current%total_weight - current%weight
                  deallocate(previous_ptrs)
                  node_deletion = .false.

                  return
               end if

               previous_ptrs(curr_lvl)%ptr => before_ptr

               curr_lvl = curr_lvl+1
               if (curr_lvl > this%max_level) then
                  exit
               end if
            else ! if next isnt lesser, go to level below
               previous_ptrs(curr_lvl)%ptr => before_ptr
               curr_lvl = curr_lvl+1
               if (curr_lvl > this%max_level) then
                  exit
               end if
            end if
         end if
      end do

      if (delete_node%weight /= weight) then
         deallocate(previous_ptrs)
         return
      end if
      node_deletion = .true.

      if (delete_node%indexes_count == 1_ik) then ! delete node
         do i = 1, this%max_level
            if (associated(previous_ptrs(i)%ptr)) then
               current => previous_ptrs(i)%ptr
               if (associated(current%next(i)%ptr)) then
                  current%next(i)%ptr => current%next(i)%ptr%next(i)%ptr !skip node
               else
                  current%next(i)%ptr => null()
               end if
            end if
         end do
         if (associated(delete_node)) then
            if (associated(delete_node, this%head%ptr)) then
               ! if (associated(this%head%ptr%next(this%max_level)%ptr)) then
               this%head%ptr => this%head%ptr%next(this%max_level)%ptr
               ! else
               !    this%head%ptr => null()
               ! end if
            end if
            deallocate(delete_node%indexes, delete_node%next)
            deallocate(delete_node)

         end if
         if (present(update_value)) update_value = 0
      end if

      deallocate(previous_ptrs)

   end subroutine skiplist_remove_entry

end module reversed_skiplist
