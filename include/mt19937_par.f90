module mt19937_par
   implicit none
   private

   integer, parameter :: dp = selected_real_kind(12, 60)

   ! MT19937 parameters
   integer, parameter :: n = 624, m = 397
   integer, parameter :: n1 = n + 1
   integer, parameter :: mata  = -1727483681
   integer, parameter :: umask = -2147483647 - 1
   integer, parameter :: lmask =  2147483647
   integer, parameter :: tmaskb = -1658038656
   integer, parameter :: tmaskc = -272236544

   type, public :: mt19937_state
      private
      integer :: mt(0:n-1)
      integer :: mti = n1
   contains
      procedure :: sgrnd
      procedure :: init_genrand
      procedure :: grnd
   end type mt19937_state

contains

!=========================================================
   subroutine sgrnd(this, seed)
      class(mt19937_state), intent(inout) :: this
      integer, intent(in)                 :: seed
      integer :: i

      this%mt(0) = iand(seed, -1)
      do i = 1, n-1
         this%mt(i) = iand(69069 * this%mt(i-1), -1)
      end do
      this%mti = n
   end subroutine sgrnd

!=========================================================
   subroutine init_genrand(this, seed)
      class(mt19937_state), intent(inout) :: this
      integer, intent(in)                 :: seed
      integer :: i, latest

      this%mt(0) = seed
      latest = seed
      do i = 1, n-1
         latest = ieor(latest, ishft(latest, -30))
         latest = latest * 1812433253 + i
         this%mt(i) = latest
      end do
      this%mti = n
   end subroutine init_genrand

!=========================================================
   function grnd(this) result(r)
      class(mt19937_state), intent(inout) :: this
      real(dp) :: r

      integer, parameter :: mag01(0:1) = (/ 0, mata /)
      integer :: y, kk

      if (this%mti >= n) then
         if (this%mti == n1) then
            call this%sgrnd(4357)
         end if

         do kk = 0, n-m-1
            y = ior(iand(this%mt(kk), umask), &
                    iand(this%mt(kk+1), lmask))
            this%mt(kk) = ieor(ieor(this%mt(kk+m), ishft(y,-1)), &
                               mag01(iand(y,1)))
         end do

         do kk = n-m, n-2
            y = ior(iand(this%mt(kk), umask), &
                    iand(this%mt(kk+1), lmask))
            this%mt(kk) = ieor(ieor(this%mt(kk+(m-n)), ishft(y,-1)), &
                               mag01(iand(y,1)))
         end do

         y = ior(iand(this%mt(n-1), umask), &
                 iand(this%mt(0), lmask))
         this%mt(n-1) = ieor(ieor(this%mt(m-1), ishft(y,-1)), &
                             mag01(iand(y,1)))

         this%mti = 0
      end if

      y = this%mt(this%mti)
      this%mti = this%mti + 1

      ! Tempering
      y = ieor(y, ishft(y, -11))
      y = ieor(y, iand(ishft(y, 7),  tmaskb))
      y = ieor(y, iand(ishft(y, 15), tmaskc))
      y = ieor(y, ishft(y, -18))

      if (y < 0) then
         r = (dble(y) + 2.0d0**32) / (2.0d0**32 - 1.0d0)
      else
         r = dble(y) / (2.0d0**32 - 1.0d0)
      end if
   end function grnd

end module mt19937_par
