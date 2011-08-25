!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all lpepcsrc-global allocation and deallocation routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_allocation
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Allocates all tree- and multipole-specific arrays,
        !> Initializes several estimations on maximum used memory
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine allocate_tree(theta)
		  use treevars
		  use timings
          use module_debug
          use module_htable
          use module_branching
		  implicit none

		  real, intent(in) :: theta

		  integer :: nintest

		  call timer_start(t_allocate)

		  if (allocated(htable)) call deallocate_tree(nppm)

		  nppm=npp
		  ! Estimate of interaction list length - Hernquist expression
		  if (theta > 0.01 ) then
		     nintest = int(35.*log(1.*npartm)/(theta**2))
		  else
		     nintest = npartm
		  endif

		  nintmax=min(nintest,npartm)

		  !  Space for # table and tree arrays
		  !  TODO: need good estimate for max # branches
		  size_tree = max(30*nintmax+4*npp,10000)

		  if (np_mult>0) then
		      nbaddr = int(max(log(1.*size_tree)/log(2.),15.))
		      maxaddress = size_tree
		  else
		     maxaddress = int(abs(np_mult)*10000)
		     nbaddr = int(max(log(1.*maxaddress)/log(2.) ,15.))
		  endif

		  if (num_pe==1) then
		    maxleaf = npart
		    maxtwig = maxaddress/2
		  else if (num_pe.lt.1024) then
		    maxleaf = maxaddress/3
		    maxtwig = 2*maxleaf
		  else
		!  required # twigs increases with P because of branches
		    maxleaf = int(maxaddress/(1.+log(1.*num_pe)/3.))
		    maxtwig = maxaddress-maxleaf
		  endif

		  hashconst = 2**nbaddr-1

		  allocate ( htable(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
		       treekey(maxaddress), branch_key(branch_max_global), branch_owner(branch_max_global), &
		       pebranch(branch_max_global), twig_key(maxtwig) )

		  free_addr = 0

		  ! Empty hashtable
		  call htable_clear()

		  ! Allocate memory for tree node properties

		  allocate ( node_level(-maxtwig:maxleaf) )

		  allocate ( tree_nodes(-maxtwig:maxleaf) )

		  call timer_stop(t_allocate)

		end subroutine allocate_tree



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Deallocates all tree- and multipole-specific arrays
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine deallocate_tree(nppm_ori)
		  use treevars
          use module_debug
          use module_htable
		  implicit none
          integer, intent(in) :: nppm_ori

		  nppm = nppm_ori

		  if (me==0 .and. tree_debug) write(*,*) 'Deallocating multipole fields'

		  deallocate ( htable, free_addr, point_free, &
		       treekey, branch_key, branch_owner, &
		       pebranch, twig_key )

		  deallocate ( node_level )

		! multipole moments
		  deallocate ( tree_nodes )

		end subroutine deallocate_tree



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Allocates all particle and results-arrays
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine allocate_particles(nppm_ori)
          use treevars
          implicit none
		  integer, intent(out) :: nppm_ori

          if (allocated(x)) call deallocate_particles()

          npartm = npart

		  if (num_pe.eq.1) then
		    nppm=int(1.5*npart + 1000)  ! allow for additional ghost particles for field plots
		  else
		    nppm = 2*max(npartm/num_pe,1000) ! allow 50% fluctuation
		  endif

		  nppm_ori = nppm

		  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)

		  ! array allocation

		  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), &
		       q(nppm), work(nppm), &
		       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

        allocate (nbranches(num_pe+2) )

        end subroutine allocate_particles



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Deallocates all particle and results-arrays
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine deallocate_particles()
		  use treevars
		  implicit none

		  ! particle array deallocation
		  deallocate ( x, y, z, ux, uy, uz, &
		       q, work, pepid, pelabel, pekey )

		  deallocate ( nbranches )

        end subroutine deallocate_particles




end module module_allocation
