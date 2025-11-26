module net_loader
    use fhash
    implicit none

    type epidemic_net
        integer, allocatable :: neighbours(:), starter_ptrs(:), end_ptrs(:), degree(:)
        type(int_int_hashmap) :: hashmap
    end type epidemic_net 
contains
    
    

    type(epidemic_net) function init_arrays(unit) result(retval)
        integer, intent(in) :: unit
        integer :: lines_count

        lines_count = count_lines(unit)
        allocate(retval%neighbours(2*lines_count), &
            retval%starter_ptrs(lines_count), &
            retval%end_ptrs(lines_count), &
            retval%degree(lines_count))



    end function init_arrays

    ! subroutine init_nodes_and_degrees(unit, net)
    !     integer, intent(in) :: unit
    !     type(epidemic_net), intent(inout) :: net
        
    !     rewind(unit)



        
    ! end subroutine init_nodes_and_degrees

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



end module net_loader