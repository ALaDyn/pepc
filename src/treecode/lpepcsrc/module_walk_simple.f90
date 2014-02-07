! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!>
!> A (not so) simple force calculation / tree traversal
!>
!> The algorithm uses OpenMP tasks to hide MPI communication latency.
!> The force calculation / tree traversal is performed for several particles grouped in a "tile" at the same time.
!>
module module_walk
  !use module_domains, only: t_decomposition
  use module_tree, only: t_tree
  use module_pepc_types, only: t_particle, kind_particle
  implicit none
  private

  type t_walk_tile
    type(t_particle), pointer :: p(:)
  end type

  type(t_tree), pointer :: walk_tree

  !type(t_decomposition) :: walk_decomposition

  integer(kind_particle), parameter :: TILE_SIZE = 64
  type(t_walk_tile), allocatable :: walk_tiles(:)
  integer(kind_particle) :: num_walk_tiles

  integer :: num_walk_threads

  real*8, public :: interactions_local, mac_evaluations_local

  public :: tree_walk_run
  public :: tree_walk_init
  public :: tree_walk_uninit
  public :: tree_walk_prepare
  public :: tree_walk_finalize
  public :: tree_walk_read_parameters
  public :: tree_walk_write_parameters
  public :: tree_walk_statistics

  contains

  subroutine tree_walk_statistics(u)
    use treevars, only: me
    implicit none

    integer, intent(in) :: u

    if (0 == me) then; write (u, *) "module_walk_simple: no statistics for now."; end if
  end subroutine tree_walk_statistics


  subroutine tree_walk_read_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine tree_walk_read_parameters


  subroutine tree_walk_write_parameters(fh)
    implicit none
    integer, intent(in) :: fh
  end subroutine tree_walk_write_parameters


  subroutine tree_walk_prepare()
    implicit none
  end subroutine tree_walk_prepare


  subroutine tree_walk_finalize()
    implicit none
  end subroutine tree_walk_finalize


  subroutine tree_walk_init(t, p_, num_threads)
    !use module_box, only: t_box, box_create
    !use module_comm_env, only: comm_env_self
    !use module_domains, only: domain_decompose
    use module_tree, only: TREE_KEY_ROOT
    !use module_spacefilling, only: level_from_key, compute_particle_keys, KEY_INVALID
    use module_spacefilling, only: level_from_key, KEY_INVALID
    use module_debug
    implicit none

    type(t_tree), target, intent(in) :: t
    type(t_particle), target, intent(in) :: p_(:)
    integer, intent(in) :: num_threads

    !type(t_box) :: b
    !type(t_particle) :: dummy(2)
    integer(kind_particle) :: ip

    call pepc_status('WALK INIT')

    walk_tree => t
    num_walk_threads = num_threads

    !if (p_(1)%key == KEY_INVALID) then
    !  call box_create(b, p_)
    !  call compute_particle_keys(b, p_)
    !  call domain_decompose(walk_decomposition, b, size(p_, kind = kind_particle), p_, dummy, comm_env_self())
    !end if

    ! worst case: every particle in its own tile
    allocate(walk_tiles(size(p_)))

    num_walk_tiles = 0
    if (p_(1)%key == KEY_INVALID) then
      num_walk_tiles = size(p_, kind = kind_particle)
      do ip = 1, num_walk_tiles
        walk_tiles(ip)%p => p_(ip:ip)
      end do
    else
      call tree_walk_init_aux(TREE_KEY_ROOT, level_from_key(TREE_KEY_ROOT), p_)
    end if
    DEBUG_ASSERT(num_walk_tiles <= size(p_))
    print *, "walk_simple using ", num_walk_tiles, " tiles."

    contains

    recursive subroutine tree_walk_init_aux(k, l, p)
      use module_pepc_types, only: kind_key, kind_level, kind_node
      use module_spacefilling, only: child_key_from_parent_key, is_ancestor_of_particle
      use treevars, only: idim, nlev
      implicit none

      integer(kind_key), intent(in) :: k
      integer(kind_level), intent(in) :: l
      type(t_particle), target :: p(:)

      integer(kind_particle) :: np, pstart, pend
      integer(kind_key) :: childkey
      integer(kind_level) :: childlevel
      integer :: ichild

      if (l == nlev) then ! no more levels left, cannot split
        DEBUG_ERROR(*, "walk_simple: insufficient key resolution in tile formation.")
      end if

      np = size(p, kind = kind_particle)

      if (np == 0) then
        ! do nothing
      else if (np <= TILE_SIZE) then
        num_walk_tiles = num_walk_tiles + 1
        walk_tiles(num_walk_tiles)%p => p(:)
      else
        childlevel = l + 1
        pstart = 1
        pend = pstart - 1
        do ichild = 0, 2**idim - 1
          childkey = child_key_from_parent_key(k, ichild)

          pend = pstart - 1
          do
            if (pend == np) exit
            if (.not. is_ancestor_of_particle(p(pend + 1)%key, childkey, childlevel)) exit
            pend = pend + 1
          end do

          call tree_walk_init_aux(childkey, childlevel, p(pstart: pend))
          pstart = pend + 1
        end do
      end if
    end subroutine tree_walk_init_aux
  end subroutine tree_walk_init


  subroutine tree_walk_uninit(t, p)
    use module_domains, only: decomposition_allocated, domain_restore
    use module_spacefilling, only: KEY_INVALID
    use module_debug
    implicit none

    type(t_tree), intent(in) :: t
    type(t_particle), intent(inout) :: p(:)

    call pepc_status('WALK UNINIT')

    deallocate(walk_tiles)

    !if (decomposition_allocated(walk_decomposition)) then
    !  call domain_restore(walk_decomposition, p)
    !  p(:)%key = KEY_INVALID
    !end if
  end subroutine tree_walk_uninit


  subroutine tree_walk_run(vbox)
    use pthreads_stuff, only: THREAD_TYPE_WORKER, THREAD_TYPE_MAIN, place_thread
    use module_debug
    use omp_lib
    implicit none

    real*8, intent(in) :: vbox(3) !< lattice vector

    integer(kind_particle) :: i

    call pepc_status('WALK SIMPLE')

    interactions_local = 0.0_8
    mac_evaluations_local = 0.0_8

    !$omp parallel default(shared) num_threads(num_walk_threads)
    call place_thread(THREAD_TYPE_WORKER, omp_get_thread_num() + 1)
    !$omp do
    do i = 1, num_walk_tiles
      !$omp task untied default(shared) firstprivate(i)
      call tree_walk_single(walk_tiles(i), vbox)
      !$omp end task
    end do
    !$omp end do nowait
    !$omp end parallel

    call place_thread(THREAD_TYPE_MAIN, 0)
  end subroutine tree_walk_run


  subroutine tree_walk_single(tl, vbox)
    use module_pepc_types, only: kind_node, kind_level
    use module_interaction_specific, only: pack_particle_list, unpack_particle_list
    use module_interaction_specific_types, only: t_particle_pack
    use module_debug
    use treevars, only: nlev
    implicit none

    type(t_walk_tile), intent(inout) :: tl
    real*8, intent(in) :: vbox(3)

    type(t_particle_pack) :: p
    real*8, allocatable :: x(:,:)
    real*8, allocatable :: r(:,:)
    real*8, allocatable :: r2(:)
    integer(kind_node) :: ni
    integer(kind_particle) :: ip, np
    integer(kind_level) :: i
    real*8 :: b2(0:nlev), num_int, num_mac

    num_int = 0.0_8
    num_mac = 0.0_8

    b2(0) = maxval(walk_tree%bounding_box%boxsize)**2
    do i = 1, nlev
      b2(i) = b2(i - 1) / 4.0_8
    end do

    np = size(tl%p)

    allocate(x(np, 3), r(np, 3), r2(np))

    do ip = 1, np
      x(ip,:) = tl%p(ip)%x - vbox
    end do

    call pack_particle_list(tl%p, p)

    DEBUG_DO(ni = 0_kind_node)
    call tree_walk_single_aux(walk_tree%node_root)
    DEBUG_ASSERT(ni == walk_tree%npart)

    call unpack_particle_list(p, tl%p)

    deallocate(x, r, r2)

    do ip = 1, np
      tl%p(ip)%work = tl%p(ip)%work + num_int
    end do

    !$omp critical
    interactions_local = interactions_local + num_int
    mac_evaluations_local = mac_evaluations_local + num_mac
    !$omp end critical

  contains

    recursive subroutine tree_walk_single_aux(n)
      use module_pepc_types, only: t_tree_node
      use module_atomic_ops, only: atomic_read_barrier
      use module_interaction_specific
      use module_tree_communicator, only: tree_node_fetch_children
      use module_tree_node, only: tree_node_children_available, tree_node_is_leaf, &
        tree_node_get_first_child, tree_node_get_next_sibling, NODE_INVALID
      use pthreads_stuff, only: pthreads_sched_yield
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
      use module_mirror_boxes, only: spatial_interaction_cutoff
      #endif
      implicit none

      integer(kind_node), intent(in) :: n

      type(t_tree_node), pointer :: node
      integer(kind_node) :: ns
      integer :: id

      node => walk_tree%nodes(n)

      r2 = 0.0_8
      do id = 1, 3
        do ip = 1, np
          r(ip, id) = x(ip, id) - node%interaction_data%coc(id)
          r2(ip) = r2(ip) + r(ip, id) * r(ip, id)
        end do
      end do

      if (tree_node_is_leaf(node)) then
        DEBUG_DO(ni = ni + 1)
        num_int = num_int + 1.0_8
        call calc_force_per_interaction_with_leaf(r, r2, p, node%interaction_data)
      else ! not a leaf, evaluate MAC
        num_mac = num_mac + 1.0_8
        if (multimac(tl%p, node%interaction_data, r2, b2(node%level))) then ! MAC OK: interact
          DEBUG_DO(ni = ni + node%leaves)
          num_int = num_int + 1.0_8
          call calc_force_per_interaction_with_twig(r, r2, p, node%interaction_data)
        else ! MAC fails: resolve
          if (.not. tree_node_children_available(node)) then
            call tree_node_fetch_children(walk_tree, node, n, tl%p(1), tl%p(1)%x - vbox)
            do ! loop and yield until children have been fetched
              ERROR_ON_FAIL(pthreads_sched_yield())
              if (tree_node_children_available(node)) exit
              !$omp taskyield
            end do
          end if

          call atomic_read_barrier()

          ns = tree_node_get_first_child(node)
          DEBUG_ASSERT_MSG(ns /= NODE_INVALID, *, "walk_simple: unexpectedly, this twig had no children")
          do
            call tree_walk_single_aux(ns)
            ns = tree_node_get_next_sibling(walk_tree%nodes(ns))
            if (ns == NODE_INVALID) exit
          end do
        end if
      end if
    end subroutine tree_walk_single_aux


    function multimac(p, n, r2, b2)
      use module_interaction_specific_types, only: t_tree_node_interaction_data
      use module_interaction_specific, only: mac
      implicit none

      logical :: multimac
      type(t_particle), intent(in) :: p(:)
      type(t_tree_node_interaction_data), intent(in) :: n
      real*8, intent(in) :: r2(:)
      real*8, intent(in) :: b2

      multimac = .true.
      do ip = 1, size(p)
        multimac = multimac .and. mac(IF_MAC_NEEDS_PARTICLE(p(ip)) n, r2(ip), b2)
        if (.not. multimac) return
      end do
    end function
  end subroutine tree_walk_single
end module module_walk
