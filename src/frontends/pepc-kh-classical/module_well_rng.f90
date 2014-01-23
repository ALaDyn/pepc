!> WELL 1024a pseudo RNG
!> http://dx.doi.org/10.1145/1132973.1132974
module module_rng
  implicit none

  private

  integer, parameter :: W = 32, R = 32, M1 = 3, M2 = 24, M3 = 10

  type well_rng_t
    integer(kind = 4), dimension(0:(R - 1)) :: state
    integer :: i

    contains
    procedure :: init => well_rng_t_init
    procedure :: next => well_rng_t_next
    procedure :: next_real => well_rng_t_next_real
  end type well_rng_t

  type(well_rng_t) :: rng

  public :: rng_init, rng_next, rng_next_real

  contains

  subroutine rng_init(i)
    implicit none

    integer(kind=4), intent(in) :: i

    rng = well_rng_t(state = 0, i = 0)
    call rng%init(i)
  end subroutine rng_init

  function rng_next()
    implicit none

    integer(kind=4) :: rng_next

    rng_next = rng%next()
  end function rng_next

  function rng_next_real()
    implicit none

    real :: rng_next_real

    rng_next_real = rng%next_real()
  end function rng_next_real

  function well_rng_t_next(this)
    class(well_rng_t), intent(inout) :: this
    integer(kind = 4) :: z0, z1, z2, well_rng_t_next

    z0 = this%state(mod(this%i + R - 1, R))
    z1 = ieor(this%state(this%i), MAT0(8, this%state(mod(this%i + M1, R))))
    z2 = ieor(MAT0(-19, this%state(mod(this%i + M2, R))), MAT0(-14, this%state(mod(this%i + M3, R))))
    this%state(this%i) = ieor(z1, z2)
    this%state(mod(this%i + R - 1, R)) = ieor(ieor(MAT0(-11, z0), MAT0(-7, z1)), MAT0(-13, z2))
    this%i = mod(this%i + R - 1, R)

    well_rng_t_next = this%state(this%i)

    contains
    pure function MAT0(t, v)
      integer, intent(in) :: t
      integer(kind = 4), intent(in) :: v
      integer(kind = 4) :: MAT0

      MAT0 = ieor(v, ishft(v, -t))
    end function MAT0

  end function well_rng_t_next

  function well_rng_t_next_real(this)
    class(well_rng_t), intent(inout) :: this
    real :: well_rng_t_next_real

    integer(kind = 4) :: i
    real, parameter :: IMI = 0.5 / huge(i)

    i = this%next()
    well_rng_t_next_real = i * IMI + 0.5
  end function well_rng_t_next_real

  subroutine well_rng_t_init(this, y0)
    class(well_rng_t), intent(inout) :: this
    integer(kind = 4), intent(in) :: y0
    integer(kind = 4) :: y
    integer :: i

    y = y0

    do i = 0, (R - 1)
      y = ieor(y, ishft(y, 13))
      y = ieor(y, ishft(y, -17))
      y = ieor(y, ishft(y, 5))

      this%state(i) = y
    end do
  end subroutine well_rng_t_init

end module module_rng
