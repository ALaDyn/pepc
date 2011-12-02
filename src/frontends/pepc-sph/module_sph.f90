!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_sph
  implicit none
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer, parameter :: num_neighbour_particles = 50      !< TODO: set this variable from pepc.f90
  
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


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief compute the kernel for sph
  !> 
  !> 

  !> \author Andreas Breslau
  !> \date 2011.12.01

  !   param[in,out]  Name      Description
  !> \param[in]      distance  
  !> \param[in]      h         
  !> \param[out]     kernel    

  subroutine sph_kernel(distance, h, kernel)

    implicit none
    include 'mpif.h'

    ! \bug ab: using idim for dimensions, be carefull with idim .ne. 3 because of pepc
    real*8, intent(in) :: distance                          !< scalar distance between two particles
    real*8, intent(in) :: h                                 !< sph smoothing length h
    real*8, intent(out) :: kernel                           !< scalar sph kernel, W

    real*8, parameter :: pi = 3.141592653589793
    real*8 :: r                                             !< absolute distance
    real*8 :: q                                             !< r/h
    real*8 :: factor

    integer :: idim                                         !< dimension, change to module variable

    integer :: ierr



    r = distance  ! should always be positive
    q = r/h

    ! TODO: idim as parameter or module variable

    idim = 3

    if( idim .eq. 3 ) then
       factor = 1._8/pi/h**3                                     ! for 3d , see monaghan 1992, S. 554, 555
    else if( idim .eq. 2 ) then
       factor = 10._8/7._8/pi/h**2
    else if( idim .eq. 1 ) then
       factor = 2._8/3._8/h
    else 
       write (*,*) "idim not one of 1, 2, 3. Terminating."

       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if


    if (q >= 0._8 .and. q <= 1._8) then
       kernel       = factor * (1._8 - 3._8/2._8 *(r/h)**2 + 3._8/4._8 * (r/h)**3)
    else if ( q <= 2._8 ) then
       kernel       = factor * (1._8/4._8 * (2._8 - (r/h) )**3)
    else ! q > 2, or negative
       write (*,*) "SPH kernel: q not in [0,2]. Should never happen:", q

       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if

  end subroutine sph_kernel


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief compute sph density
  !> 
  !> 

  !> \author Andreas Breslau
  !> \date 2011.12.01

  !   param[in,out]  Name      Description
  !> \param[in]      
  !> \param[in]      
  !> \param[in]      
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sph_density(np_local, particles, &
       itime, num_neighbour_boxes, neighbour_boxes)

    use treetypes, only: &
         t_particle

    use treevars, only: &
         tree_nodes
    
!    use physvars, only: &
!         n_cpu, &
!         my_rank
    
!    use module_htable, only: & 
!         key2addr
    
!    use module_spacefilling, only: &
!         coord_to_key_lastlevel, &
!         key_to_coord

    implicit none
    include 'mpif.h'

    
    integer, intent(in) :: np_local    !< # particles on this CPU
    type(t_particle), intent(inout) :: particles(:)
    integer, intent(in) :: itime  ! timestep
    integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
    integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
    
    integer :: ierr

    integer :: local_particle_index
    integer :: actual_neighbour
    real*8 :: h
    real*8 :: kernel
    integer*8 :: actual_node

!     IF( me == 0) then 
       
!        write(87, *) me, idim
      
!        DO actual_particle = 1, num_part
!           WRITE(87, *) part_indices( actual_particle ), r_nn( part_indices( actual_particle ) ), pelabel( part_indices( actual_particle ) )
!           DO actual_neighbour = 1, n_nn
!              WRITE(87, *) charge( nodelist( actual_neighbour, actual_particle ) ), dist2_list( actual_neighbour, actual_particle )
!           END DO
!        END DO
       
!     end if


    ! compute smoothing length (h) for all local particles
    ! TODO: move this out of sph_density
    ! TODO: parallelize with OpenMP
    do local_particle_index = 1, np_local
       particles(local_particle_index)%results%h = sqrt(particles(local_particle_index)%results%maxdist2)/2._8
    end do


    do local_particle_index = 1, np_local

       h = particles(local_particle_index)%results%h

       call sph_kernel( 0._8, h, kernel) ! particle self
       particles(local_particle_index)%results%rho = particles(local_particle_index)%data%q *kernel

       do actual_neighbour = 1, num_neighbour_particles
          actual_node = particles(local_particle_index)%results%neighbour_nodes(actual_neighbour)

          call sph_kernel(  sqrt( particles(local_particle_index)%results%dist2(actual_neighbour) ), h, kernel)

          particles(local_particle_index)%results%rho = particles(local_particle_index)%results%rho + &
               tree_nodes(actual_node)%q

       end do

    end do

  end subroutine sph_density


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \brief compute nabla_kernel for sph
  !> 
  !> 

  !> \author Andreas Breslau
  !> \date 2011.12.01

  !   param[in,out]  Name      Description
  !> \param[in]      distance   
  !> \param[in]      h         
  !> \param[out]     nabla_kernel

  subroutine sph_grad_kernel(distance, h, grad_kernel)

    implicit none
    include 'mpif.h'

    ! \bug ab: using idim for dimensions, be carefull with idim .ne. 3 because of pepc
    real*8, intent(in) :: distance                          !< scalar distance between two particles
    real*8, intent(in) :: h                                 !< sph smoothing length h
    real*8, intent(out) :: grad_kernel                      !< scalar part of gradient of sph kernel, Nabla W,
    ! has to be multiplied by the distance vector by the calling function

    real*8, parameter :: pi = 3.141592653589793
    real*8 :: r                                             !< absolute distance
    real*8 :: q                                             !< r/h
    real*8 :: factor
    
    integer :: ierr
    integer :: idim

    ! TODO: make idim module variable
    idim = 3

    r = distance             ! should always be positive
    q = r/h
    
    ! TODO: make kernel factor a module variable and initialise it once!
    if( idim .eq. 3 ) then
       factor = 1._8/pi/h**3                                     ! for 3d , see monaghan 1992, S. 554, 555
    else if( idim .eq. 2 ) then
       factor = 10._8/7._8/pi/h**2
    else if( idim .eq. 1 ) then
       factor = 2._8/3._8/h
    else 
       write (*,*) "idim not one of 1, 2, 3. Terminating."

       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if


    if (q >= 0._8 .and. q <= 1._8) then
       grad_kernel = factor * (9._8 * r - 12._8 * h )/(4._8 * h**3 )
    else if ( q > 1._8 .and. q <= 2._8 ) then
       grad_kernel = factor * (-3._8) * (2._8 -  r/h )**2 / (4._8 * h * r )
    else ! q > 2, or negative
       write (*,*) "SPH grad kernel: q not in [0,2]. Should never happen:", q

       ! \bug ab: contribution of last particle is always 0 besause h is defined as half of the distance to this particle 

       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if

  end subroutine sph_grad_kernel




  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

!  subroutine sph_sum_forces(num_part, part_indices, nppm_ori, ex_tmp, ey_tmp, ez_tmp, thermal_constant, art_vis_alpha, art_vis_beta, temperature_change_tmp, kappa, sph_factor)

    use treetypes, only: &
         t_particle

    use treevars, only: &
         tree_nodes


    implicit none
    include 'mpif.h'

    
    integer, intent(in) :: np_local    !< # particles on this CPU
    type(t_particle), intent(inout) :: particles(:)
    integer, intent(in) :: itime  ! timestep
    integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
    integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
    
    integer :: ierr

    integer :: local_particle_index
    integer :: actual_neighbour
    real*8 :: h
    real*8 :: grad_kernel
    integer*8 :: actual_node

    INTEGER :: actual_particle
    REAL*8 :: const
    REAL*8 :: T
    REAL*8 :: xdist, ydist, zdist
    REAL*8 :: force_contribution                  ! \bug ab: only for testing
    REAL*8 :: artificial_viscosity
    REAL*8 :: vr, rr, mu, eta
    REAL*8 :: sound_speed
    REAL*8 :: scalar_force
    REAL*8 :: thermal_energy_sum
    REAL*8 :: thermal_energy_factor
    REAL*8 :: dvx, dvy, dvz
    REAL*8 :: energy_factor

    CHARACTER(100) :: forcefile
    LOGICAL :: use_artificial_viscosity

    real*8 :: thermal_constant             ! TODO: make this a parameter
    real*8 :: kappa                        ! TODO: make this a parameter
    integer :: idim                        ! TODO: make this a parameter
    real*8 :: art_vis_alpha                ! TODO: make this a parameter
    real*8 :: art_vis_beta                 ! TODO: make this a parameter



    thermal_constant = 1.
    kappa = 1.
    idim = 1

    
    const = real(thermal_constant,8)  ! m**2 s**-2 K
    ! p = rho *const*T , quick n dirty


!write(forcefile,'(a,i6.6,a,i6.6,a)') "sph_", itime-1, "_", me, ".list"
!open(50, FILE=forcefile)



    ! thermal energy change:
    ! du_i/dt = 1/2 SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) (v_i -v_j ) grad_kernel ) (Monaghan 1992 S. 549)
    ! dT/dt = (gamma -1) /k_B du/dt
    ! gamma = kappa (adiabatic exponent)
    ! \todo ab: correction factor f_i (siehe Springel 2010, SEREN 2011) ??

    ! acceleration:
    ! dv_i/dt = -SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) grad_kernel ) (Monaghan 1992 S. 546)

    ! for artificial viscosity (art_vis) see Monaghan 1992, S. 550


    ! TODO: make this a parameter
    use_artificial_viscosity = .true.

    
    
    do local_particle_index = 1, np_local
       
       h = particles(local_particle_index)%results%h
       
!write(50+me, *) "new part:", part_indices( actual_particle ), pelabel( part_indices( actual_particle ) ), tree_x( part_indices( actual_particle) ), h
       
       thermal_energy_sum = 0._8


       do actual_neighbour = 1, num_neighbour_particles
          
          actual_node = particles(local_particle_index)%results%neighbour_nodes(actual_neighbour)
          
          call sph_grad_kernel(  sqrt( particles(local_particle_index)%results%dist2(actual_neighbour) ), h, grad_kernel)
          ! grad_kernel is negative, so the minus in the force turns to a plus

          
!write(50+me, *) '  ', sqrt( dist2_list(actual_neighbour, actual_particle))
          
          
          ! PEPCs sum_force returns the electric field components ex, ey, ez. In velocities the accelerations are then computed as
          ! a = charge / mass * e(x|y|z). For gravity charge/mass = 1, so the fields and accelerations are equal. 
          ! Because of this the accelerations caused by sph_force can be added to the fields.
          
          ! always at least one dimensional
          
          xdist = tree_nodes(actual_node)%coc(1) &
               ! TODO:       + periodic_shift_vectors( shift_list( actual_neighbour, actual_particle ) , 1 ) &
               - particles(local_particle_index)%x(1)

          dvx = tree_nodes(actual_node)%v(1) - particles(local_particle_index)%data%v(1)
          
          ! vr: art_vis, scalar product of velocity difference and distance
          vr =  xdist * dvx
          
          ! rr: art_vis, distance squared
          rr =  particles(local_particle_index)%results%dist2(actual_neighbour)
          
          thermal_energy_factor = vr
          
          if( idim > 1 ) then
             ydist = tree_nodes(actual_node)%coc(2) &
                  ! TODO:       + periodic_shift_vectors( shift_list( actual_neighbour, actual_particle ) , 1 ) &
                  - particles(local_particle_index)%x(2)
             
             dvy = tree_nodes(actual_node)%v(2) - particles(local_particle_index)%data%v(2)
             
             vr = vr + ydist * dvy
             
             thermal_energy_factor = thermal_energy_factor + dvy *ydist
             
             
             if (idim > 2 ) then
                zdist = tree_nodes(actual_node)%coc(3) &
                     ! TODO:       + periodic_shift_vectors( shift_list( actual_neighbour, actual_particle ) , 1 ) &
                     - particles(local_particle_index)%x(3)
                
                dvz = tree_nodes(actual_node)%v(3) - particles(local_particle_index)%data%v(3)
                
                vr = vr + zdist * dvz
                
                thermal_energy_factor = thermal_energy_factor + dvz * zdist
                
             end if
          end if
          
          
          eta = real(0.1,8) * h                                                ! art_vis
          sound_speed = ( sqrt( const * particles(local_particle_index)%data%temperature )  &
               + sqrt( const * tree_nodes(actual_node)%temperature ) )/ 2. ! mean sound_speed 
          
          if ( use_artificial_viscosity .and. (vr < 0._8) ) then
             
             mu = ( h * vr ) / ( rr + eta * eta )                        ! art_vis
             
             artificial_viscosity = ( - art_vis_alpha * sound_speed * mu + art_vis_beta * mu * mu ) / &
                  ( ( tree_nodes(actual_node)%temperature +  particles(local_particle_index)%results%rho )/2. )
             
          else 
             
             artificial_viscosity = 0._8
          end if
          
          
          ! compute scalar part of the force: mass * ( p1/rho1^2 + p2/rho2^2 + art_vis ) * grad_kernel
          scalar_force = tree_nodes(actual_node)%q * &
               ( &
               const * tree_nodes(actual_node)%temperature / tree_nodes(actual_node)%rho + &
               const * particles(local_particle_index)%data%temperature / particles(local_particle_index)%results%rho + &
               artificial_viscosity &
               ) * grad_kernel
          ! \todo ab: use both kernels here (i, j)
          
          particles(local_particle_index)%results%sph_force(1) = &
               particles(local_particle_index)%results%sph_force(1) + scalar_force * xdist
          
          if( idim > 1 ) then
             
             particles(local_particle_index)%results%sph_force(2) = &
                  particles(local_particle_index)%results%sph_force(2) + scalar_force * ydist
             
             if (idim > 2 ) then
                
                particles(local_particle_index)%results%sph_force(3) = &
                     particles(local_particle_index)%results%sph_force(3) + scalar_force * ydist
                
             end if
          end if
          
!write(50+me, *) scalar_force, vr, grad_kernel, xdist, dvx, artificial_viscosity

          thermal_energy_sum = thermal_energy_sum + vr * scalar_force
          
       end do ! end of loop over neighbours

!       temperature_change_tmp( part_indices( actual_particle ) ) = real( pelabel( part_indices( actual_particle) ), 8 ) 

       
       energy_factor = 1._8 * ( kappa - 1._8 ) /const /2._8 
       
       particles(local_particle_index)%results%sph_force(2) = energy_factor * thermal_energy_sum
       ! TODO: add art_visc term to temp_change
       
!write(50+me, *) temperature_change_tmp( part_indices( actual_particle ) )
       
    end do
!close(50)
    
  end subroutine sph_sum_force


!   subroutine sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

!     use treetypes, only: &
!          t_particle

!     use treevars, only: &
!          tree_nodes


!     implicit none
!     include 'mpif.h'

    
!     integer, intent(in) :: np_local    !< # particles on this CPU
!     type(t_particle), intent(inout) :: particles(:)
!     integer, intent(in) :: itime  ! timestep
!     integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
!     integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
    
!     integer :: ierr


! !  subroutine sph_force(nppm_ori, ex_tmp, ey_tmp, ez_tmp, n_nn_tmp, nshortm, itime_tmp, periodic_x, periodic_y, periodic_z, &
! !       boxlength_x, boxlength_y, boxlength_z, idim_tmp, max_npass_tmp, nshort_tmp, pstart_tmp, thermal_constant, art_vis_alpha, art_vis_beta!, temperature_change_tmp, kappa, sph_factor )


! !    use treevars, &
! !         ONLY: &
! !         nlev, &
! !         boxsize, &
! !         xmin, ymin, zmin, &
! !         ntwig, &
! !         ntwigp, &
! !         nleaf, &
! !         x, &     ! \bug ab: for testing
! !         xcoc     ! \bug ab: for testing


!     implicit none
!     include 'mpif.h'


!     INTEGER, INTENT(in) :: n_nn_tmp
!     INTEGER, INTENT(in) :: nshortm
!     INTEGER, INTENT(in) :: itime_tmp
!     LOGICAL, INTENT(in) :: periodic_x
!     LOGICAL, INTENT(in) :: periodic_y
!     LOGICAL, INTENT(in) :: periodic_z
!     REAL*8, INTENT(in) :: boxlength_x
!     REAL*8, INTENT(in) :: boxlength_y
!     REAL*8, INTENT(in) :: boxlength_z
!     INTEGER, INTENT(in) :: idim_tmp
!     INTEGER, INTENT(in) :: max_npass_tmp
!     INTEGER, INTENT(in), DIMENSION(max_npass_tmp) :: nshort_tmp
!     INTEGER, INTENT(in), DIMENSION(max_npass_tmp) :: pstart_tmp
!     REAL, INTENT(in) :: sph_factor
    
! !    INTEGER, INTENT(in) :: npshort                          !< number of particles in current chunk
! !    INTEGER, DIMENSION(npshort), INTENT(in) :: pshort       !< array with particle indices
! !    INTEGER, INTENT(in) :: pass
!     INTEGER, INTENT(in) :: nppm_ori
!     REAL*8, DIMENSION(nppm_ori), INTENT(inout) :: ex_tmp,ey_tmp,ez_tmp
!     REAL*8, INTENT(in) :: thermal_constant
!     REAL, INTENT(in) :: art_vis_alpha
!     REAL, INTENT(in) :: art_vis_beta
!     REAL, INTENT(in) :: kappa
!     REAL*8, DIMENSION(nppm_ori), INTENT(out) :: temperature_change_tmp
    
!     INTEGER :: nps
!     INTEGER :: ip1
!     REAL :: eps
!     INTEGER :: i
!     INTEGER :: jpass
!     REAL*8 :: ttrav_nn, tfetch_nn            ! timing for search_nn

!     REAL*8 :: s
!     INTEGER*8 :: ix, iy, iz
!     INTEGER :: nbits
!     INTEGER*8 :: local_key

!     REAL*8 kernel
!     REAL*8 dist

!     CHARACTER(100) :: nn_filename



!     call sph_density(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

!     call update_particle_props(np_local, particles)
  
!     call sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)
    
!   end subroutine sph


  subroutine update_particle_props(np_local, particles)
    
    use treetypes, only: &
         t_particle
    
    use treevars, only: &
         tree_nodes, &
         nleaf, &
         nleaf_me
    
    use module_htable
    
    ! only for sort test
    use tree_utils
    
    
    use physvars, only: &
         my_rank, &
         n_cpu
    
    implicit none
    include 'mpif.h'
    
    
    integer, intent(in) :: np_local    !< # particles on this CPU
    type(t_particle), intent(inout) :: particles(:)
    
    integer :: nleaf_non_local
    integer*8, allocatable :: non_local_node_keys(:)
    integer*8, allocatable :: key_arr_cp(:)
    integer*8, allocatable :: non_local_node_owner(:)
    integer, allocatable :: int_arr(:)
    integer :: num_request
    integer :: i
    integer :: ierr
    integer, allocatable :: requests_per_process(:)
 
    nleaf_non_local = nleaf ! bigger than necessary, TODO: find a better estimation for this

    allocate( non_local_node_keys(nleaf_non_local), non_local_node_owner(nleaf_non_local), requests_per_process(n_cpu), &
         key_arr_cp(nleaf_non_local), int_arr(nleaf_non_local), STAT=ierr )
    ! TODO: remove key_arr_cp and int_arr after sort test
    ! TODO: test STAT

    num_request = 0

    ! get leafs from hashtabel with owner .ne. my_rank
    do i = 1, maxaddress
       if( htable_entry_is_valid(i) ) then
          if( (htable(i)%owner .ne. my_rank) .and. htable(i)%node>0 ) then
             if( htable(i)%owner .ne. mod(my_rank+1,2)) write(*,*) 'strange owner:', my_rank, htable(i)%owner, htable(i)%key

             num_request = num_request + 1
             non_local_node_keys(num_request) = htable(i)%key
             non_local_node_owner(num_request) = htable(i)%owner
          end if
       end if
    end do
    
    if( nleaf_non_local > num_request) write (*,*) 'on rank', my_rank, 'nleaf_non_local:', nleaf_non_local, 'num_request:', num_request


    ! write (*,*) my_rank, num_request
    
    ! ! test sort with keys
    ! write(*,*) ' starting sort test'
    ! key_arr_cp = non_local_node_keys
    
    ! call sort(non_local_node_keys, int_arr)
    
    ! do i= 1, nleaf_non_local
    !    if( key_arr_cp(int_arr(i)) .ne. non_local_node_keys(i) ) write(*,*) 'Error in sort.'
    ! end do
    
    ! do i= 1, nleaf_non_local-1
    !    if( non_local_node_keys(i) > non_local_node_keys(i+1)) write(*,*) 'Error in sort.'
    ! end do
  
    ! write(*,*) 'sort test finished'


    ! sort keys accorting to owner

    do i = 1, num_request
       if(non_local_node_owner(i) > n_cpu -1 ) write(*,*) 'owner > n_cpu:', i, non_local_node_owner(i), non_local_node_keys(i)
    end do


    
    call sort(non_local_node_owner(1:num_request), int_arr(1:num_request))
    
    key_arr_cp(1:num_request) = non_local_node_keys(1:num_request)

    do i= 1, num_request
       non_local_node_keys(i) = key_arr_cp(int_arr(i))
    end do

    requests_per_process = 0

    do i = 1, num_request
       ! use owner + 1, because owner is from 0 and array index from 1
       requests_per_process(non_local_node_owner(i)+1) = requests_per_process(non_local_node_owner(i)+1) + 1
    end do

    if( requests_per_process(my_rank+1) .ne. 0) write (*,*) 'on rank', my_rank, 'requests for self is non-zero:', & 
         requests_per_process(my_rank+1)


!    call MPI_ALLTOALL(






    deallocate( non_local_node_keys, non_local_node_owner, requests_per_process, key_arr_cp, int_arr, STAT=ierr )

    
    
    
    
    
    
  end subroutine update_particle_props
  


end module module_sph

