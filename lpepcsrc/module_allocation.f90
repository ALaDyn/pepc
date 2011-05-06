!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all lpepcsrc-global allocation and deallocation routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_allocation
      use module_debug
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
		  implicit none

		  real, intent(in) :: theta

		  integer :: k, nintest

		  call timer_start(t_allocate)

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

		  if (num_pe > 1) then
		     nbranch_max = int(.75*maxaddress)         ! Global max # branches
		     nbranch_local_max = 2*nbranch_max/num_pe  ! Local max # branches
		  else
		     nbranch_max = 5*nintmax
		     nbranch_local_max = 5*nintmax
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

		  if (me==0 .and. tree_debug) then
		    write(*,'(//a/)') 'Allocating new multipole fields'
		    write(*,*) '# procs',num_pe
		    write(*,*) 'npart=',npart
		    write(*,*) 'N/P=',npart/num_pe
		    write(*,*) 'npp= ',npp
		    write(*,*) 'nppm= ',nppm
		    write(*,*) 'nintest/max=',nintest,nintmax
		    write(*,*) 'size_tree= ',size_tree
		    write(*,*) 'max address = ',maxaddress
		    write(*,*) 'address bits = ',nbaddr
		    write(*,*) '# const = ',hashconst
		    write(*,*) 'max leaf = ',maxleaf
		    write(*,*) 'max twig = ',maxtwig
		    write(*,*) 'max branches = ',nbranch_max
		    write(*,*) 'np_mult= ',np_mult
		    write(*,'(a/)') '... done'
		  endif

		  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
		       treekey(maxaddress), branch_key(nbranch_max), branch_owner(nbranch_max), &
		       pebranch(nbranch_max), leaf_key(maxaddress), twig_key(maxaddress) )

		  all_addr = (/ (k,k=0,maxaddress) /)      ! List of all possible # table addresses

		  free_addr = 0

		  ! Empty hashtable
		  htable = hash(0,0_8,-1,0,0,0_8,0)

		  ! Allocate memory for tree node properties

		  allocate ( first_child(-maxtwig:maxleaf), node_level(-maxtwig:maxleaf) )

		  allocate ( charge(-maxtwig:maxleaf), &                    ! charge
		       abs_charge(-maxtwig:maxleaf), &                ! absolute charge
		       xcoc(-maxtwig:maxleaf), ycoc(-maxtwig:maxleaf), zcoc(-maxtwig:maxleaf), &    ! centre of charge
		       xshift(-maxtwig:maxleaf), yshift(-maxtwig:maxleaf), zshift(-maxtwig:maxleaf), &    ! shift vector
		       size_node(-maxtwig:maxleaf), & ! rms node size
		       xdip(-maxtwig:maxleaf), ydip(-maxtwig:maxleaf), zdip(-maxtwig:maxleaf), &          ! dipole moment
		       xxquad(-maxtwig:maxleaf), yyquad(-maxtwig:maxleaf), zzquad(-maxtwig:maxleaf), &       ! quadrupole moment
		       xyquad(-maxtwig:maxleaf), yzquad(-maxtwig:maxleaf), zxquad(-maxtwig:maxleaf))

		  call timer_stop(t_allocate)

		end subroutine allocate_tree



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Deallocates all tree- and multipole-specific arrays
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine deallocate_tree(nppm_ori)
		  use treevars
		  implicit none
          integer, intent(in) :: nppm_ori

		  nppm = nppm_ori

		  if (me==0 .and. tree_debug) write(*,*) 'Deallocating multipole fields'

		  deallocate ( htable, all_addr, free_addr, point_free, &
		       treekey, branch_key, branch_owner, &
		       pebranch, leaf_key, twig_key )

		  deallocate ( first_child, node_level )

		! multipole moments
		  deallocate ( charge, &                    ! charge
		       abs_charge, &                ! absolute charge
		       xcoc, ycoc, zcoc, size_node, &    ! centre of charge
		       xshift, yshift, zshift, &    ! shift vector
		       xdip, ydip, zdip, &          ! dipole moment
		       xxquad, yyquad, zzquad, &       ! quadrupole moment
		       xyquad, yzquad, zxquad)

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

          npartm = npart

		  if (num_pe.eq.1) then
		    nppm=int(1.5*npart + 1000)  ! allow for additional ghost particles for field plots
		  else
		    nppm = 2*max(npartm/num_pe,1000) ! allow 50% fluctuation
		  endif

		  nppm_ori = nppm

		  nlev = 20                     ! max refinement level
		  iplace = 2_8**(3*nlev)           ! place holder bit
		  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)

		  ! array allocation

		  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), &
		       q(nppm), m(nppm), work(nppm), &
		       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

		  allocate (nbranches(num_pe+2), igap(num_pe+3))

        end subroutine allocate_particles



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Deallocates all particle and results-arrays
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine deallocate_particles()
		  use treevars
		  implicit none

		  if (me==0) then
		     write(*,'(a)') 'LPEPC | De-allocating particle arrays ...'
		  endif

		  ! particle array deallocation
		  deallocate ( x, y, z, ux, uy, uz, &
		       q, m, work, &
		       pepid, pelabel, pekey )

		  deallocate ( nbranches, igap )

        end subroutine deallocate_particles




end module module_allocation
