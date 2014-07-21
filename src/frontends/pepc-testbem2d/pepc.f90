module bem
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  private

  type(t_particle), allocatable, target :: bem_el(:)
  real*8, allocatable :: bem_tmp(:)
  integer :: bem_iteration

  public :: bem_init
  public :: bem_solve
  public :: bem_uninit
  public :: bem_matmult

  contains

  subroutine bem_init()
    implicit none

    #include <finclude/petsc.h90>

    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  end subroutine bem_init


  subroutine bem_uninit()
    implicit none

    #include <finclude/petsc.h90>

    PetscErrorCode :: ierr

    call PetscFinalize(ierr)
  end subroutine bem_uninit


  subroutine bem_solve(el, part)
    use module_pepc
    use module_pepc_types
    use module_interaction_specific_types
    implicit none

    #include <finclude/petsc.h90>

    type(t_particle), intent(inout) :: el(:)
    type(t_particle), optional, allocatable, target, intent(inout) :: part(:)

    Mat :: A
    Vec :: x, b
    PetscInt :: m
    PetscErrorCode :: ierr
    KSP :: ksp
    PetscScalar, pointer :: sol(:)

    m = size(el)
    allocate(bem_el(m))
    allocate(bem_tmp(m))
    bem_el = el
    bem_iteration = 0

    call MatCreateShell(PETSC_COMM_SELF, m, m, m, m, PETSC_NULL_OBJECT, A, ierr)
    call MatShellSetOperation(A, MATOP_MULT, bem_matmult, ierr)

    call VecCreateSeq(PETSC_COMM_SELF, m, x, ierr)
    call VecDuplicate(x, b, ierr)

    ! initialize rhs
    call bem_rhs(b, part)

    ! initialize solution vector to previous/user-supplied solution
    call VecGetArrayF90(x, sol, ierr)
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET) sol = bem_el(:)%data%q
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)   sol = bem_el(:)%data%phi
    call VecRestoreArrayF90(x, sol, ierr)

    ! mask out knowns
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET)
      bem_tmp = bem_el(:)%data%phi
      bem_el(:)%data%phi = 0.0_8
    end where
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)
      bem_tmp = bem_el(:)%data%q
      bem_el(:)%data%q = 0.0_8
    end where

    call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
    call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
    call KSPSetFromOptions(ksp, ierr)
    call KSPSolve(ksp, b, x, ierr)

    ! write solution to boundary element particles and restore knowns
    call VecGetArrayF90(x, sol, ierr)
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET)
      bem_el(:)%data%phi = bem_tmp
      bem_el(:)%data%q   = sol
    end where
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)
      bem_el(:)%data%phi = sol
      bem_el(:)%data%q   = bem_tmp
    end where
    call VecRestoreArrayF90(x, sol, ierr)

    call KSPDestroy(ksp, ierr)
    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)

    el = bem_el
    deallocate(bem_tmp)
    deallocate(bem_el)
  end subroutine bem_solve


  subroutine bem_rhs(b, part)
    use module_pepc
    use module_interaction_specific
    use module_interaction_specific_types
    implicit none

    #include <finclude/petsc.h90>

    Vec :: b
    type(t_particle), optional, allocatable, target, intent(inout) :: part(:)

    real*8, parameter :: pi = 3.1415926535897932385_8

    PetscScalar, pointer :: res(:)
    PetscErrorCode :: ierr
    PetscInt :: i

    call pepc_particleresults_clear(bem_el)

    ! compute particular part RHS
    if (present(part)) then
      call pepc_grow_and_traverse_for_others(part, bem_el, -1, no_dealloc = .false., no_restore = .false.)
    end if

    ! mask out the unknowns for RHS computation
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET)
      bem_tmp = bem_el(:)%data%q
      bem_el(:)%data%q = 0.0_8
    end where
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)
      bem_tmp = bem_el(:)%data%phi
      bem_el(:)%data%phi = 0.0_8
    end where

    call pepc_grow_and_traverse(bem_el, -1, no_dealloc = .false., no_restore = .false.)

    ! restore unknowns to particles
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET) bem_el(:)%data%q   = bem_tmp
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)   bem_el(:)%data%phi = bem_tmp

    call VecGetArrayF90(b, res, ierr)
    res = bem_el(:)%results%pot
    call VecRestoreArrayF90(b, res, ierr)
  end subroutine bem_rhs


  subroutine bem_matmult(A, x, y, ierr)
    use module_pepc
    use module_interaction_specific
    use module_interaction_specific_types
    use module_walk, only: interactions_local
    implicit none

    #include <finclude/petsc.h90>

    Mat :: A
    Vec :: x
    Vec :: y
    PetscErrorCode :: ierr

    real*8, parameter :: pi = 3.1415926535897932385_8

    PetscInt :: i
    PetscScalar, pointer :: lhs(:), res(:)

    write (*, '("bem_matmult iteration: ", I0)', advance = 'NO') bem_iteration

    !write (*, *) "LHS vector:"
    !call VecView(x, PETSC_VIEWER_STDOUT_SELF, ierr)

    ! assign values given by LHS vec x to bem particles
    call VecGetArrayF90(x, lhs, ierr)
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_DIRICHLET) bem_el(:)%data%q   = lhs
    where (bem_el(:)%data%source_kind == CALC_FORCE_SOURCE_KIND_NEUMANN)   bem_el(:)%data%phi = lhs
    call VecRestoreArrayF90(x, lhs, ierr)

    call pepc_particleresults_clear(bem_el)
    call pepc_grow_and_traverse(bem_el, bem_iteration, no_dealloc = .false., no_restore = .false.)

    write (*, *) ", interactions: ", interactions_local

    ! store results in vector y
    call VecGetArrayF90(y, res, ierr)
    res = -bem_el(:)%results%pot
    call VecRestoreArrayF90(y, res, ierr)

    !write (*, *) "Result vector:"
    !call VecView(y, PETSC_VIEWER_STDOUT_SELF, ierr)

    bem_iteration = bem_iteration + 1
  end subroutine bem_matmult
end module bem



module writers
  use module_interaction_specific_types
  use module_vtk
  implicit none

  contains

  subroutine write_results(d, r, vtkf)
    implicit none

    type(t_particle_data), intent(in) :: d(:)
    type(t_particle_results), intent(in) :: r(:)
    type(vtkfile_unstructured_grid), intent(inout) :: vtkf

    call vtkf%write_data_array("phi", r(:)%pot)
    call vtkf%write_data_array("ex", r(:)%e(1))
    call vtkf%write_data_array("ey", r(:)%e(2))
  end subroutine write_results

  subroutine write_bc(d, r, vtkf)
    implicit none

    type(t_particle_data), intent(in) :: d(:)
    type(t_particle_results), intent(in) :: r(:)
    type(vtkfile_unstructured_grid), intent(inout) :: vtkf

    call vtkf%write_data_array("phi", d(:)%phi)
    call vtkf%write_data_array("q", d(:)%q)
  end subroutine write_bc
end module writers



program pepc
  use, intrinsic :: iso_fortran_env
  use module_pepc
  use module_pepc_types
  use module_interaction_specific
  use module_interaction_specific_types
  use treevars
  use module_walk
  use module_vtk
  use module_vtk_helpers

  use writers
  use bem
  implicit none

  include 'mpif.h'

  integer(kind_particle), parameter :: NX = 128, NY = NX, NSIDE = 128, NPART = 128

  type(t_particle), target, allocatable :: pgrid(:), pside(:), part(:)

  integer(kind_particle) :: ix, iy, is, ip
  integer(kind_default) :: mpi_comm
  integer(kind_pe) :: mpi_rank, mpi_size
  integer :: fd

  real*8, parameter :: pi = 3.1415926535897932385_8
  real*8 :: l, dx, dy

  integer, dimension(2,3) :: globaldims, mydims
  real*8, allocatable :: xcoords(:), ycoords(:), zcoords(:)
  real*8, allocatable :: scalarvalues(:, :, :), vectorvalues(:, :, :, :)

  call pepc_initialize("pepc-testbem2d", mpi_rank, mpi_size, .true., comm = mpi_comm)
  num_threads = 8
  theta2 = 0.09_8
  call pepc_prepare(2_kind_dim)
  call pepc_write_parameters(output_unit)

  call bem_init()

  allocate(part(NPART))
  open (file = 'part.dat', status = 'replace', newunit = fd)
  do ip = 1, NPART
    associate (p => part(ip))
      call random_number(p%x(1:2))
      p%x(3) = 0.0_8
      p%work = 1.0_8
      p%data%v = [ 0.0_8, 0.0_8, 0.0_8 ]
      p%data%m = 0.0_8

      p%data%phi = 0.0_8
      p%data%q = 1.0_8 / NPART
      if (0 == mod(ip, 2)) then; p%data%q = p%data%q * (-1.0_8); end if

      p%data%source_kind = CALC_FORCE_SOURCE_KIND_PARTICULAR
      write (fd, '(3(g0,:,","))') p%x(1), p%x(2), p%data%q
    end associate
  end do
  close (fd)

  ! Set up sides
  allocate(pside(4 * NSIDE))
  open (file = 'boundary.dat', status = 'replace', newunit = fd)
  l = 1.0_8 / NSIDE
  do is = 1, 4
    do ip = 1, NSIDE
      associate (p => pside(ip + (is - 1) * NSIDE))
        p%work = 1.0_8
        p%data%v = [ 0.0_8, 0.0_8, 0.0_8 ]
        p%data%m = 0.0_8

        select case (is)
        case (1)
          p%data%ra = [ (ip - 1) * l, 0.0_8 ]
          p%data%rb = [ ip * l, 0.0_8 ]
          p%data%phi = 0.0_8
          p%data%q = -1.0_8
          p%data%source_kind = CALC_FORCE_SOURCE_KIND_DIRICHLET
        case (2)
          p%data%ra = [ 1.0_8, (ip - 1) * l ]
          p%data%rb = [ 1.0_8, ip * l ]
          p%data%phi = 1.0_8
          p%data%q = 0.0_8
          p%data%source_kind = CALC_FORCE_SOURCE_KIND_NEUMANN
        case (3)
          p%data%ra = [ (NSIDE - ip + 1) * l, 1.0_8 ]
          p%data%rb = [ (NSIDE - ip) * l, 1.0_8 ]
          p%data%phi = 1.0_8
          p%data%q = 1.0_8
          p%data%source_kind = CALC_FORCE_SOURCE_KIND_DIRICHLET
        case (4)
          p%data%ra = [ 0.0_8, (NSIDE - ip + 1) * l ]
          p%data%rb = [ 0.0_8, (NSIDE - ip) * l ]
          p%data%phi = 1.0_8
          p%data%q = 0.0_8
          p%data%source_kind = CALC_FORCE_SOURCE_KIND_NEUMANN
        end select

        p%x(1:2) = (p%data%ra + p%data%rb) / 2.0_8
        p%x(3) = 0.0_8
        write (fd, '(7(g0,:,","))') p%data%ra, p%data%rb, p%data%phi, p%data%q, p%data%source_kind
      end associate
    end do
  end do
  close (fd)

  call bem_solve(pside, part)

  open (file = 'boundary_solution.dat', status = 'replace', newunit = fd)
  do ip = 1, 4 * NSIDE
    associate (p => pside(ip))
      write (fd, '(7(g0,:,","))') p%data%ra, p%data%rb, p%data%phi, p%data%q, 0
    end associate
  end do
  close (fd)

  ! Set up grid
  allocate(pgrid(NX * NY))
  dx = 1.0_8 / NX
  dy = 1.0_8 / NY
  do ix = 1, NX
    do iy = 1, NY
      associate (p => pgrid(iy + NY * (ix - 1)))
        p%work = 1.0_8
        p%data%phi = 0.0_8
        p%data%q = 0.0_8
        p%data%v = [ 0.0_8, 0.0_8, 0.0_8 ]
        p%data%m = 0.0_8

        p%data%source_kind = CALC_FORCE_SOURCE_KIND_PARTICULAR

        p%x(1) = (ix - 0.5_8) * dx
        p%x(2) = (iy - 0.5_8) * dy
        p%x(3) = 0.0_8
      end associate
    end do
  end do

  call pepc_particleresults_clear(pgrid)
  call pepc_grow_tree(pside)
  call pepc_traverse_tree(pgrid)
  call vtk_write_leaves(0, 0.0_8, VTK_STEP_FIRST, global_tree)
  call pepc_timber_tree()

  call pepc_grow_tree(part)
  call pepc_traverse_tree(pgrid)
  call pepc_timber_tree()

  pgrid(:)%results%pot = pgrid(:)%results%pot / (2 * pi)
  pgrid(:)%results%e(1) = pgrid(:)%results%e(1) / (2 * pi)
  pgrid(:)%results%e(2) = pgrid(:)%results%e(2) / (2 * pi)

  open (file = 'result.dat', status = 'replace', newunit = fd)
  do ip = 1, NX * NY
    associate (p => pgrid(ip))
      write (fd, '(3(g0,:,","))') p%x(1:2), p%results%pot
    end associate
  end do
  close (fd)

  call vtk_write_particles("grid", mpi_comm, 0, 0.0_8, VTK_STEP_FIRST, pgrid, write_results)
  call vtk_write_particles("boundary", mpi_comm, 0, 0.0_8, VTK_STEP_FIRST, pside, write_bc)
  call vtk_write_particles("boundary_results", mpi_comm, 0, 0.0_8, VTK_STEP_FIRST, pside, write_results)

  allocate (xcoords(NX), ycoords(NY), zcoords(1))
  allocate (scalarvalues(NX, NY, 1), vectorvalues(NX, NY, 1, 3))
  globaldims(1, :) = [ 1, 1, 1 ]
  globaldims(2, :) = [ NX, NY, 1_8 ]
  mydims = globaldims

  xcoords = [ ((ix - 0.5_8) * dx, ix = 1, NX) ]
  ycoords = [ ((iy - 0.5_8) * dy, iy = 1, NY) ]
  zcoords = [ 0.0 ]

  do ix = 1, NX
    do iy = 1, NY
      associate (p => pgrid(iy + NY * (ix - 1)))
        scalarvalues(ix, iy, 1) = p%results%pot
        vectorvalues(ix, iy, 1, 1:2) = p%results%e(1:2)
        vectorvalues(ix, iy, 1, 3) = 0.0_8
      end associate
    end do
  end do

  call vtk_write_field_on_grid("grid", 0, 0.0_8, VTK_STEP_FIRST, globaldims, mydims, xcoords, ycoords, zcoords, &
    scalarvalues, "phi", vectorvalues, "e", MPI_COMM_WORLD)
  deallocate (xcoords, ycoords, zcoords)
  deallocate (scalarvalues, vectorvalues)

  deallocate(pside, pgrid, part)

  call bem_uninit()
  call pepc_finalize(mpi_comm)
end program
