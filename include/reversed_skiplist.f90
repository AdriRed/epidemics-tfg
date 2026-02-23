
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
      real(dp) :: weight
      integer(ik) :: index
      type(skiplist_entry_ptr), allocatable :: next(:)
   end type skiplist_entry

   type skiplist
      type(mt19937_state), pointer :: rng
      integer(ik) :: max_level
      type(skiplist_entry_ptr), pointer :: head
      real(dp) :: promote_probability

   contains
      procedure :: clear => skiplist_clear
      procedure :: add => insert_skiplist
      procedure :: search_greater_or_equal => find_previous_node_weight_greater_or_equal
      procedure :: debug_print => print_skiplist
   end type skiplist
contains

   subroutine print_skiplist(this)
      class(skiplist), intent(inout) :: this
      type(skiplist_entry), pointer :: previous
      integer(ik) :: i
      do i = 1, this%max_level
         previous => this%head%ptr
         write(*, "(A11, I1, A5)") "---- Level ", i, " ----"
         do while(associated(previous))
            write(*, "(F5.1, A4)", advance="no") previous%weight, ' -> '
            previous => previous%next(i)%ptr
         end do
         print *, 'NULL'
      end do
   end subroutine

   function init_skiplist(maxlevels, rng, probability) result(retval)
      integer(ik), intent(in) :: maxlevels
      real(dp), intent(in) :: probability
      type(mt19937_state), intent(in), target :: rng
      type(skiplist) :: retval
      retval%rng => rng
      retval%max_level = maxlevels
      retval%promote_probability = probability
      allocate(retval%head)

   end function init_skiplist

   subroutine insert_skiplist(this, weight, index)
      class(skiplist), intent(inout) :: this
      real(dp), intent(in) :: weight
      integer(ik), intent(in) :: index
      integer(ik) :: level, i
      real(dp) :: r
      type(skiplist_entry), pointer :: data, previous, old_head
      logical :: insert_head_link
      allocate(data)
      data%index = index
      data%weight = weight
      allocate(data%next(this%max_level))
      insert_head_link = .true.
      if (.not. associated(this%head%ptr)) then
         this%head%ptr => data
      else if (this%head%ptr%weight <= weight) then ! case where is going to be first
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

      else
         previous => find_previous_node_weight_greater_or_equal(this, weight)
         ! promoting
         level = this%max_level
         do
            data%next(level)%ptr => previous%next(level)%ptr ! new node pointing now to old next node
            previous%next(level)%ptr => data ! old node pointing to new node
            r = this%rng%grnd()
            if (r < this%promote_probability) then
               level = level - 1
               if (level == 0) exit
            else
               exit
            end if
         end do
      end if

   end subroutine insert_skiplist

   function find_previous_node_weight_greater_or_equal(this, weight) result(before_ptr)
      class(skiplist), intent(inout) :: this
      real(dp), intent(in) :: weight
      type(skiplist_entry), pointer :: before_ptr
      integer(ik) :: curr_lvl


      curr_lvl = 1
      before_ptr => this%head%ptr
      do
         if (.not. associated(before_ptr%next(curr_lvl)%ptr)) then ! end of the line
            curr_lvl = curr_lvl+1
            if (curr_lvl > this%max_level) exit ! final level reached
            cycle ! if not final level reached keep searching
         else
            if (weight <= before_ptr%next(curr_lvl)%ptr%weight) then ! if next lesser
               before_ptr => before_ptr%next(curr_lvl)%ptr ! get next in line
            else ! if next isnt lesser, go to level below
               curr_lvl = curr_lvl+1
               if (curr_lvl > this%max_level) then
                  exit
               end if
            end if
         end if
      end do

   end function find_previous_node_weight_greater_or_equal

   subroutine skiplist_clear(this)
      class(skiplist), intent(inout) :: this
      type(skiplist_entry), pointer :: current
      type(skiplist_entry), pointer :: next_ptr

      if (.not. associated(this%head)) return

      current => this%head%ptr

      do while (associated(current))
         next_ptr => current%next(this%max_level)%ptr
         deallocate(current%next)
         deallocate(current)
         current => next_ptr
      end do

      nullify(this%head%ptr)
      deallocate(this%head)

      this%max_level = 0
   end subroutine skiplist_clear



end module reversed_skiplist
