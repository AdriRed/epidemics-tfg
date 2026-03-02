program test_skiplist
   use iso_fortran_env, only: int8, int32, real64
   use reversed_skiplist
   implicit none


   integer, parameter :: dp = real64
   integer, parameter :: ik = int32
   integer, parameter :: bk = int8

   type(skiplist) :: list
   integer, parameter :: max_levels = 4
   integer, parameter :: seed = 12345
   real(dp), parameter :: probability = 0.5_dp
   integer, parameter :: max_indexes = 3
   integer(ik) :: i, position, new_value, update_value
   real(dp) :: weight
   logical :: is_empty, full_deletion

! Inicializar la lista con parámetros extremos
   print *, "Testing initialization with extreme parameters..."
   list = init_skiplist(max_levels, seed, probability, max_indexes)

   ! Prueba de adición de nodos con pesos extremos
   print *, "Adding nodes with extreme weights..."
   weight = huge(1.0_dp)  ! Peso muy grande
   call list%add(weight, 1, position)
   print *, "Added very large weight:", weight, "at position:", position

   weight = tiny(1.0_dp)  ! Peso muy pequeño
   call list%add(weight, 2, position)
   print *, "Added very small weight:", weight, "at position:", position

   ! Imprimir la lista después de añadir nodos con pesos extremos
   print *, "List after adding nodes with extreme weights:"
   call list%debug_print()

   ! Verificar que los pesos están en orden decreciente después de añadir nodos con pesos extremos
   print *, "Verifying order after adding nodes with extreme weights..."
   if (is_sorted_descending(list)) then
      print *, "List is correctly sorted in descending order."
   else
      print *, "List is NOT correctly sorted in descending order."
   end if

   ! Prueba de adición de nodos con el número máximo de índices
   print *, "Adding nodes with maximum number of indexes..."
   do i = 1, max_levels
      weight = 10.0_dp - i * 1.0_dp
      call list%add(weight, i, position)
      print *, "Added weight:", weight, "with", max_indexes, "indexes at position:", position
   end do

   ! Imprimir la lista después de añadir nodos con el número máximo de índices
   print *, "List after adding nodes with maximum indexes:"
   call list%debug_print()

   ! Verificar que los índices son correctos
   print *, "Verifying indexes..."
   if (are_indexes_correct(list)) then
      print *, "Indexes are correct."
   else
      print *, "Indexes are NOT correct."
   end if

   ! Prueba de adición de nodos con pesos duplicados
   print *, "Adding nodes with duplicate weights..."
   weight = 5.0_dp
   call list%add(weight, 10, position)
   print *, "Added duplicate weight:", weight, "at position:", position

   ! Imprimir la lista después de añadir nodos con pesos duplicados
   print *, "List after adding nodes with duplicate weights:"
   call list%debug_print()

   ! Prueba de eliminación de nodos que no existen
   print *, "Attempting to remove non-existent nodes..."
   weight = 100.0_dp  ! Peso que no existe en la lista
   call list%remove_entry(weight, position, full_deletion, update_value=update_value)
   print *, "Attempted to remove non-existent weight:", weight, "at position:", position, "update value:", update_value

   ! Imprimir la lista después de intentar eliminar un nodo que no existe
   print *, "List after attempting to remove non-existent node:"
   call list%debug_print()

   ! Prueba de eliminación de todos los nodos
   print *, "Removing all nodes..."
   ! Primero, eliminar los nodos con pesos extremos
   weight = huge(1.0_dp)
   call list%remove_entry(weight, position, full_deletion, update_value=update_value)
   print *, "Removed weight:", weight, "at position:", position, "update value:", update_value

   weight = tiny(1.0_dp)
   call list%remove_entry(weight, position, full_deletion, update_value=update_value)
   print *, "Removed weight:", weight, "at position:", position, "update value:", update_value
    call list%debug_print()

   ! Luego, eliminar los nodos con pesos intermedios
   do i = 1, max_levels
      weight = 10.0_dp - i * 1.0_dp
      call list%remove_entry(weight, position, full_deletion, update_value=update_value)
      print *, "Removed weight:", weight, "at position:", position, "update value:", update_value
      call list%debug_print()
   end do

   ! Finalmente, eliminar el nodo con peso duplicado
   weight = 5.0_dp
   call list%remove_entry(weight, position, full_deletion, update_value=update_value)
   print *, "Removed weight:", weight, "at position:", position, "update value:", update_value

   ! Imprimir la lista después de eliminar todos los nodos
   print *, "List after removing all nodes:"
   call list%debug_print()

   ! Verificar que la lista está vacía
   print *, "Checking if list is empty..."
   is_empty = .not. associated(list%head%ptr)
   if (is_empty) then
      print *, "List is empty."
   else
      print *, "List is NOT empty."
   end if
contains

   ! Función para verificar que los pesos están en orden decreciente
   logical function is_sorted_descending(list) result(sorted)
      type(skiplist), intent(in) :: list
      type(skiplist_entry), pointer :: current
      real(dp) :: prev_weight
      integer(ik) :: level

      sorted = .true.
      level = 1  ! Verificar en el nivel más bajo (slow lane)
      current => list%head%ptr
      if (.not. associated(current)) then
         sorted = .true.
         return
      end if
      prev_weight = current%weight
      do while (associated(current%next(level)%ptr))
         current => current%next(level)%ptr
         if (current%weight > prev_weight) then
            sorted = .false.
            exit
         end if
         prev_weight = current%weight
      end do
   end function is_sorted_descending

   ! Función para verificar que los índices son correctos
   logical function are_indexes_correct(list) result(correct)
      type(skiplist), intent(in) :: list
      type(skiplist_entry), pointer :: current
      integer(ik) :: i, j

      correct = .true.
      current => list%head%ptr
      do while (associated(current))
         ! Verificar que los índices están dentro del rango permitido
         do i = 1, current%indexes_count
            if (current%indexes(i) < 1 .or. current%indexes(i) > max_indexes) then
               correct = .false.
               return
            end if
         end do
         ! Verificar que no hay índices duplicados en un nodo
         do i = 1, current%indexes_count
            do j = i+1, current%indexes_count
               if (current%indexes(i) == current%indexes(j)) then
                  correct = .false.
                  return
               end if
            end do
         end do
         current => current%next(1)%ptr
      end do
   end function are_indexes_correct

end program test_skiplist
