! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

module pfasst_calc_module




contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit part.
  subroutine eval_f1(y, t, Nvar, level, f1)

    use physvars
    use module_interaction_specific
    use module_pepc
    use pfasst_helper_module
    use manipulate_particles
    implicit none

    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y(Nvar), t
    real(kind=8), intent(out) :: f1(Nvar)

    integer :: i

    ! reshape y to vortex_particles
    call pfasst_to_pepc(vortex_particles(1:np), np, y(1:Nvar))

    call pepc_particleresults_clear(vortex_particles, np)

    call direct_sum(np, vortex_particles, vortex_particles%results, my_rank_space, n_cpu_space)

    !call pepc_grow_and_traverse(np, n, vortex_particles, 1, .false., .false.)

    do i=1,np
       vortex_particles(i)%results%u( 1:3) = vortex_particles(i)%results%u( 1:3) * force_const
       vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
    end do

    ! reshape vortex_particles to f1
    call pepc_to_pfasst_res(vortex_particles(1:np), np, f1(1:Nvar))

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit part.
  subroutine eval_f2(y, t, Nvar, level, f2)
    implicit none

    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y(Nvar), t
    real(kind=8), intent(out) :: f2(Nvar)

    f2 = 0.0

  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(y_start, t, Nvar, level, yexact)
    integer,      intent(in ) :: Nvar, level
    real(kind=8), intent(in ) :: y_start(Nvar), t
    real(kind=8), intent(out) :: yexact(Nvar)

    yexact = 0.0

  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sdcsweep(y_start, t0, dt, ySDC, fSDC, Nvar, Nnodes, level, tau)

    use pfasst_helper_module
    implicit none

    !  Arguments
    integer,      intent(in   ) :: Nvar, Nnodes, level
    real(kind=8), intent(in   ) :: dt, t0
    real(kind=8), intent(in   ) :: y_start(Nvar)
    real(kind=8), intent(inout) :: ySDC(Nvar,Nnodes)
    real(kind=8), intent(inout) :: fSDC(Nvar,Nnodes,2)
    real(kind=8), intent(in   ), optional :: tau(Nvar,Nnodes-1)

    !  Local stuff
    integer      :: m
    real(kind=8) :: t
    real(kind=8) :: dtsdc(1:Nnodes-1)
    real(kind=8) :: Smat(1:Nnodes-1,1:Nnodes)
    real(kind=8) :: Setil(1:Nnodes-1,1:Nnodes)
    real(kind=8) :: Sitil(1:Nnodes-1,1:Nnodes)
    real(kind=8) :: qnodes(Nnodes)
    real(kind=8) :: Irhs(Nvar,1:Nnodes-1)
    real(kind=8) :: f1(Nvar), f2(Nvar), rhs(Nvar), ynew(Nvar)

    !  Set integration matrices and nodes
    select case(level)
    case (1)
       Smat   = SmatF
       Setil  = StilLF
       Sitil  = StilnoLF
       qnodes = qnodesF
    case (2)
       Smat   = SmatG
       Setil  = StilLG
       Sitil  = StilnoLG
       qnodes = qnodesG
    end select

    !  Do the integration
    Irhs = dt*( matmul(fSDC(:,:,1),transpose(Smat-Setil)) + &
                matmul(fSDC(:,:,2),transpose(Smat-Sitil)) )

    !  Add FAS correction if provided
    if (present(tau)) then
       Irhs = Irhs + tau
    end if

    !  Do the timstepping
    ySDC(:,1) = y_start

    call eval_f1(y_start, t0, Nvar, level, f1)
    fSDC(:,1,1) = f1

    call eval_f2(y_start, t0, Nvar, level, f2)
    fSDC(:,1,2) = f2

    t = t0
    dtsdc = dt*( qnodes(2:Nnodes)-qnodes(1:Nnodes-1) )
    do m = 1, Nnodes-1
       t = t + dtsdc(m)

       rhs = ySDC(:,m) + dtsdc(m)*f1 + Irhs(:,m)

       ynew = ySDC(:,m+1)
       call comp_f2(ynew, t, dtsdc(m), rhs, Nvar, level, f2)
       call eval_f1(ynew, t, Nvar, level, f1)

       ySDC(:,m+1)   = ynew
       fSDC(:,m+1,1) = f1
       fSDC(:,m+1,2) = f2
    end do

  end subroutine sdcsweep

end module pfasst_calc_module
