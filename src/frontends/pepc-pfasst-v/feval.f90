module feval
contains

  subroutine init_feval(filename)
    implicit none
    character(len=*), intent(in) :: filename
  end subroutine init_feval

  ! Evaluate the explicit part.
  subroutine eval_f1(y, t, Nvar, level, f1)
    !use RHSFunctions

    implicit none

    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y(Nvar), t
    real(kind=8), intent(out) :: f1(Nvar)

    integer :: counter

    ! reshape y to vortex_particles
    counter = 1
    do i = 1,np

       vortex_particles(i)%x(1) = y(counter+0)
       vortex_particles(i)%x(2) = y(counter+1)
       vortex_particles(i)%x(3) = y(counter+2)
       vortex_particles(i)%data%alpha(1) = y(counter+3)
       vortex_particles(i)%data%alpha(2) = y(counter+4)
       vortex_particles(i)%data%alpha(3) = y(counter+5)

       counter = counter + 6

    end do

    call pepc_grow_and_traverse(np, n, vortex_particles, itime, .false., .false.)

    do i=1,np
       vortex_particles(i)%results%u( 1:3) = vortex_particles(i)%results%u( 1:3) * force_const
       vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
    end do

    ! reshape vortex_particles to f1
    counter = 1
    do i = 1,np

       f1(counter+0) = vortex_particles(i)%results%u(1)
       f1(counter+1) = vortex_particles(i)%results%u(2)
       f1(counter+2) = vortex_particles(i)%results%u(3)
       f1(counter+3) = vortex_particles(i)%results%af(1)
       f1(counter+4) = vortex_particles(i)%results%af(2)
       f1(counter+5) = vortex_particles(i)%results%af(3)

       counter = counter + 6

    end do

  end subroutine eval_f1

  ! Evaluate the implicit part.
  subroutine eval_f2(y, t, Nvar, level, f2)
    implicit none

    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y(Nvar), t
    real(kind=8), intent(out) :: f2(Nvar)

    f2 = 0.0

  end subroutine eval_f2

  ! Solve the implicit part.
  subroutine comp_f2(y, t, Dt, rhs, Nvar, level, f2)
    implicit none
    integer,      intent(in   ) :: Nvar, level
    real(kind=8), intent(inout) :: y(Nvar)
    real(kind=8), intent(in   ) :: rhs(Nvar), Dt, t
    real(kind=8), intent(out  ) :: f2(Nvar)

    f2 = 0.0
    y  = rhs

  end subroutine comp_f2

  subroutine exact(y0, t, Nvar, level, yexact)
    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y0(Nvar), t
    real(kind=8), intent(out) :: yexact(Nvar)

    yexact = 0.0

  end subroutine exact

end module feval
