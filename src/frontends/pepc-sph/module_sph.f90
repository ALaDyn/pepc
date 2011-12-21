!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_sph

  use module_interaction_specific_types, only: &
       num_neighbour_particles

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
  
  integer :: idim = 3

  real*8 :: thermal_constant             
  real*8 :: kappa                        
  real*8 :: art_vis_alpha                
  real*8 :: art_vis_beta                 
  logical :: use_artificial_viscosity

  real*8, parameter :: pi = 3.1415926535897932384626433832795028842   ! 64-bit
  
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



  subroutine sph_initialize(idim_tmp, thermal_constant_tmp, kappa_tmp, art_vis_alpha_tmp, art_vis_beta_tmp)
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: idim_tmp
    real*8, intent(in) :: thermal_constant_tmp
    real*8, intent(in) :: kappa_tmp
    real*8, intent(in) :: art_vis_alpha_tmp
    real*8, intent(in) :: art_vis_beta_tmp

    integer :: ierr
    
    
    ! set module variables
    idim             = idim_tmp
    thermal_constant = thermal_constant_tmp
    kappa            = kappa_tmp
    art_vis_alpha    = art_vis_alpha_tmp
    art_vis_beta     = art_vis_beta_tmp

    use_artificial_viscosity = .false.


  end subroutine sph_initialize






  subroutine sph_kernel_tests(idim_)
    
    implicit none
    
    integer, intent(in) :: idim_

    real*8 :: r
    real*8 :: h
    real*8 :: kernel
    real*8 :: grad_kernel
    integer :: i, steps


    idim = idim_

    steps = 100
    h = 1.

    do i = 0, steps
       
       r = -2. + 4./steps*i
       call sph_kernel(abs(r),h,kernel)
       call sph_grad_kernel(abs(r), h, grad_kernel)
       
       grad_kernel = grad_kernel * r

       write(77, *) i, r, kernel, grad_kernel

    end do

    close(77)

  end subroutine sph_kernel_tests










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

  subroutine sph_kernel(r, h, kernel)
    
    implicit none
    include 'mpif.h'
    
    ! \bug ab: using idim for dimensions, be carefull with idim .ne. 3 because of pepc
    real*8, intent(in) :: r                                 !< scalar distance between two particles, should allways be >= 0
    real*8, intent(in) :: h                                 !< sph smoothing length h
    real*8, intent(out) :: kernel                           !< scalar sph kernel, W
    
    real*8 :: q                                             !< r/h
    real*8 :: sph_kernel_factor
    
    integer :: ierr
    
    q = r/h

    if( idim .eq. 3 ) then 
      sph_kernel_factor = 1._8/pi/h**3                                     ! see monaghan 1992, S. 554, 555
    else if( idim .eq. 2 ) then
       sph_kernel_factor = 10._8/7._8/pi/h**2
    else if( idim .eq. 1 ) then
       sph_kernel_factor = 2._8/3._8/h
    else 
       write (*,*) "idim not one of 1, 2, 3. Terminating:", idim
       
       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if


    
    ! see monaghan 1992, S. 554, 555
    if (q >= 0._8 .and. q <= 1._8) then
       kernel       = sph_kernel_factor * (1._8 - 3._8/2._8 *(r/h)**2 + 3._8/4._8 * (r/h)**3)
    else if ( q <= 2._8 ) then
       kernel       = sph_kernel_factor * (1._8/4._8 * (2._8 - (r/h) )**3)
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

    use module_pepc_types, only: &
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
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index)
    do local_particle_index = 1, np_local
       particles(local_particle_index)%results%h = sqrt(particles(local_particle_index)%results%maxdist2)/2._8
    end do
    !$OMP END PARALLEL DO

    
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index, h, kernel, actual_neighbour, actual_node)
    do local_particle_index = 1, np_local
       
       h = particles(local_particle_index)%results%h

       call sph_kernel( 0._8, h, kernel) ! particle self
       particles(local_particle_index)%results%rho = particles(local_particle_index)%data%q *kernel

       do actual_neighbour = 1, num_neighbour_particles
          actual_node = particles(local_particle_index)%results%neighbour_nodes(actual_neighbour)

          call sph_kernel(  sqrt( particles(local_particle_index)%results%dist2(actual_neighbour) ), h, kernel)

          particles(local_particle_index)%results%rho = particles(local_particle_index)%results%rho + &
               tree_nodes(actual_node)%q *kernel

       end do

    end do
    !$OMP END PARALLEL DO


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

  subroutine sph_grad_kernel(r, h, grad_kernel)

    implicit none
    include 'mpif.h'
    
    ! \bug ab: using idim for dimensions, be carefull with idim .ne. 3 because of pepc
    real*8, intent(in) :: r                                 !< scalar distance between two particles, should be >= 0
    real*8, intent(in) :: h                                 !< sph smoothing length h
    real*8, intent(out) :: grad_kernel                      !< scalar part of gradient of sph kernel, Nabla W,
    ! has to be multiplied by the distance vector by the calling function

    real*8 :: q                                             !< r/h
    real*8 :: sph_kernel_factor

    integer :: ierr

    q = r/h


    if( idim .eq. 3 ) then
       sph_kernel_factor = 1._8/pi/h**3                                     ! see monaghan 1992, S. 554, 555
    else if( idim .eq. 2 ) then
       sph_kernel_factor = 10._8/7._8/pi/h**2
    else if( idim .eq. 1 ) then
       sph_kernel_factor = 2._8/3._8/h
    else 
       write (*,*) "idim not one of 1, 2, 3. Terminating."
       
       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if

    
    if (q >= 0._8 .and. q <= 1._8) then
       grad_kernel = sph_kernel_factor * (9._8 * r - 12._8 * h )/(4._8 * h**3 )
    else if ( q > 1._8 .and. q <= 2._8 ) then
       grad_kernel = sph_kernel_factor * (-3._8) * (2._8 -  r/h )**2 / (4._8 * h * r )
    else ! q > 2, or negative
       write (*,*) "SPH grad kernel: q not in [0,2]. Should never happen:", q
       
       ! ab: contribution of last particle is always 0, because h is defined as half of the distance to this particle, so we are really using n-1 neighbours
       
       call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
    end if

  end subroutine sph_grad_kernel




  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! subroutine to sum the sph acceleration, not the force, becaue the tree_walk returns the field,
  ! which, in case of gravity, equals the acceleration.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

    
    use module_pepc_types, only: &
         t_particle

    use treevars, only: &
         tree_nodes

    use physvars, only: &
         my_rank

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
    real*8 :: h1, h2                          !< smoothing length of particle and current interaction partner
    real*8 :: grad_kernel_1, grad_kernel_2    !< kernel gradient for particle and current interaction partner with h1 and h2
    integer*8 :: actual_node

    real*8 :: const
    real*8, dimension(3) :: dist              !< distance vector from particle to current interaction partner
    real*8 :: distance                        !< scalar distance from particle to current interaction partner
    real*8 :: artificial_viscosity
    real*8 :: vr, rr, mu, eta
    real*8 :: sound_speed
    real*8 :: scalar_force
    real*8 :: thermal_energy_sum
    real*8 :: thermal_energy_factor
    real*8, dimension(3) :: dv                !< velocity difference between particle and current interaction partner
    real*8 :: energy_factor
    integer :: dim


    CHARACTER(100) :: forcefile


    do local_particle_index = 1, np_local
       particles(local_particle_index)%results%sph_force = [ 0._8, 0._8, 0._8 ]
    end do
    
    
    const = thermal_constant  ! m**2 s**-2 K

    ! p = rho *const*T , quick n dirty
    
    
    energy_factor = 1._8 * ( kappa - 1._8 ) /const /2._8 


!write(forcefile,'(a,i6.6,a,i6.6,a)') "sph_", itime, "_", me, ".list"
!open(50, FILE=forcefile)


    ! thermal energy change:
    ! du_i/dt = 1/2 SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) (v_i -v_j ) grad_kernel ) (Monaghan 1992 S. 549)
    ! dT/dt = (gamma -1) /k_B du/dt
    ! gamma = kappa (adiabatic exponent)
    ! \todo ab: correction factor f_i (siehe Springel 2010, SEREN 2011) ??
    
    ! acceleration:
    ! dv_i/dt = -SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) grad_kernel ) (Monaghan 1992 S. 546)
    
    ! for artificial viscosity (art_vis) see Monaghan 1992, S. 550
    
    

    
    
    ! _$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index, h1, h2, thermal_energy_sum, actual_neighbour, actual_node, grad_kernel_1, grad_kernel_2, &
    ! _$OMP& rr, vr, dim, dist, dv, thermal_energy_factor, eta, sound_speed, mu, artificial_viscosity, scalar_force, distance)
    do local_particle_index = 1, np_local
       
       h1 = particles(local_particle_index)%results%h
       
!write(50+me, *) "new part:", part_indices( actual_particle ), pelabel( part_indices( actual_particle ) ), tree_x( part_indices( actual_particle) ), h
       
       thermal_energy_sum = 0._8


       do actual_neighbour = 1, num_neighbour_particles
          
          actual_node = particles(local_particle_index)%results%neighbour_nodes(actual_neighbour)
          
          ! scalar distance
          distance = sqrt( particles(local_particle_index)%results%dist2(actual_neighbour) )

          h2 = tree_nodes(actual_node)%h
          
          call sph_grad_kernel( distance , h1, grad_kernel_1 )
!          call sph_grad_kernel( distance , h2, grad_kernel_2 )
          ! grad_kernel is negative, so the minus in the force turns to a plus
          
          
          !write(50+me, *) '  ', sqrt( dist2_list(actual_neighbour, actual_particle))
          
          
          ! PEPCs sum_force returns the electric field components ex, ey, ez. In velocities the accelerations are then computed as
          ! a = charge / mass * e(x|y|z). For gravity charge/mass = 1, so the fields and accelerations are equal. 
          ! Because of this the accelerations caused by sph_force can be added to the fields.
          
          ! distance vector
          dist = particles(local_particle_index)%results%dist_vector(:, actual_neighbour)

          if( abs(sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2) - distance ) > 1E-7 ) then
             write(*,*) '|dist| - distance > 1E-7:', my_rank, local_particle_index, actual_neighbour, distance, dist
          end if






          ! rr: art_vis, distance squared
          rr =  particles(local_particle_index)%results%dist2(actual_neighbour)
          
          vr = 0._8

          ! for all dimension
          do dim = 1, idim
             dv(dim) = tree_nodes(actual_node)%v(dim) - particles(local_particle_index)%data%v(dim)
             
             ! vr: art_vis, scalar product of velocity difference and distance
             vr =  vr + dist(dim) * dv(dim)
          end do
          
          thermal_energy_factor = vr
          
          
          ! TODO: make eta parameter???
          eta = 0.1_8 * h1                                                ! art_vis
          sound_speed = ( sqrt( const * particles(local_particle_index)%data%temperature )  &
               + sqrt( const * tree_nodes(actual_node)%temperature ) )/ 2. ! mean sound_speed 
          
          if ( use_artificial_viscosity .and. (vr < 0._8) ) then
             
             mu = ( h1 * vr ) / ( rr + eta * eta )                        ! art_vis
             
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
               ) * grad_kernel_1
! + grad_kernel_2 ) / 2._8
          

          particles(local_particle_index)%results%sph_force = particles(local_particle_index)%results%sph_force + scalar_force * dist

          write(78, *) my_rank, local_particle_index, h1, particles(local_particle_index)%data%temperature, particles(local_particle_index)%results%rho, actual_neighbour, distance, grad_kernel_1, tree_nodes(actual_node)%q, tree_nodes(actual_node)%temperature, tree_nodes(actual_node)%h, tree_nodes(actual_node)%rho, scalar_force*dist(1), particles(local_particle_index)%results%sph_force(1)

          
!write(50+me, *) scalar_force, vr, grad_kernel, xdist, dvx, artificial_viscosity

          thermal_energy_sum = thermal_energy_sum + vr * scalar_force
          
       end do ! end of loop over neighbours
       
       
       particles(local_particle_index)%results%temperature_change = energy_factor * thermal_energy_sum
       ! TODO: add art_visc term to temp_change
       
       !write(50+me, *) temperature_change_tmp( part_indices( actual_particle ) )
       
    end do
    ! _$OMP END PARALLEL DO
      
    !close(50)
    
  end subroutine sph_sum_force
  
  


  subroutine sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes, idim)
    
    use module_pepc_types, only: &
         t_particle
    
    use treevars, only: &
         tree_nodes

    use physvars, only: &
         thermal_constant

    implicit none
    include 'mpif.h'

    
    integer, intent(in) :: np_local    !< # particles on this CPU
    type(t_particle), intent(inout) :: particles(:)
    integer, intent(in) :: itime  ! timestep
    integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
    integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
    integer, intent(in) :: idim    ! dimension

    integer :: ierr

    call sph_initialize(idim, thermal_constant, 1.4_8, 2._8, 1._8)

    call sph_density(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)
    
    call update_particle_props(np_local, particles)
    
    do ierr = 1, np_local
       particles(ierr)%results%sph_force = 0.
    end do

    call sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)


  end subroutine sph







  subroutine update_particle_props(np_local, particles)
    
    use module_pepc_types, only: &
         t_particle
    
    use treevars, only: &
         tree_nodes, &
         nleaf, &
         nleaf_me
    
    use module_htable
    
    ! only for sort test
    use module_utils

    use physvars, only: &
         my_rank, &
         n_cpu

    implicit none
    include 'mpif.h'

    ! Data structure for shipping updated sph properties
    type t_property_update
       integer*8 :: key                                                  !< key
       integer   :: owner                                                !< owner
       real*8    :: smoothing_length                                     !< \bug ab: comments needed
       real*8    :: rho                                                  !<
       real*8    :: v(1:3)                                               !< velocity
       REAL*8    :: temperature                                          !< SPH temperature
    end type t_property_update

    integer, parameter :: nprops_property_update = 6
    
        
    integer, intent(in) :: np_local    !< # particles on this CPU
    type(t_particle), intent(inout) :: particles(:)
    
    integer :: nleaf_non_local
    integer*8, allocatable :: non_local_node_keys(:)
    integer*8, allocatable :: key_arr_cp(:)
    integer*8, allocatable :: non_local_node_owner(:)
    integer, allocatable :: non_local_node_node(:)
    integer, allocatable :: node_arr_cp(:)
    integer, allocatable :: int_arr(:)
    integer :: num_request
    integer :: i
    integer :: ierr
    integer, allocatable :: requests_per_process(:)
    integer, allocatable :: requests_from_process(:)
    integer*8, allocatable :: requested_keys(:)
    integer, allocatable :: sdispls(:)
    integer, allocatable :: rdispls(:)
    integer :: total_num_requests_from_others
    integer :: disp
    integer :: actual_address
    integer :: actual_node

    type(t_property_update), allocatable :: packed_updates(:)
    type(t_property_update), allocatable :: received_updates(:)

    type(t_property_update)  :: dummy_property_update
    integer :: mpi_type_property_update
    integer, parameter :: max_props = nprops_property_update
    ! address calculation
    integer, dimension(1:max_props) :: blocklengths, displacements, types
    integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
        
   
    nleaf_non_local = nleaf - nleaf_me ! bigger than necessary, TODO: find a better estimation for this

    allocate( non_local_node_node(nleaf_non_local), non_local_node_keys(nleaf_non_local), non_local_node_owner(nleaf_non_local), requests_per_process(n_cpu), &
         key_arr_cp(nleaf_non_local), node_arr_cp(nleaf_non_local), int_arr(nleaf_non_local), requests_from_process(n_cpu), &
         sdispls(n_cpu), rdispls(n_cpu), STAT=ierr )
    ! TODO: remove key_arr_cp and int_arr after sort test
    ! TODO: test STAT

    num_request = 0

    ! get leafs from hashtabel with owner .ne. my_rank
    ! TODO: do not search in htable, but in array containing nodes ? (Idee von Lukas)
    do i = 1, maxaddress
       if( htable_entry_is_valid(i) ) then
          if( (htable(i)%owner .ne. my_rank) .and. htable(i)%node>0 ) then
             if( htable(i)%owner > n_cpu-1) write(*,*) 'strange owner:', my_rank, htable(i)%owner, htable(i)%key

             num_request = num_request + 1
             non_local_node_keys(num_request) = htable(i)%key
             non_local_node_owner(num_request) = htable(i)%owner
             non_local_node_node(num_request) = htable(i)%node
          end if
       end if
    end do
    
!    if( nleaf_non_local > num_request) write (*,*) 'on rank', my_rank, 'nleaf_non_local:', nleaf_non_local, 'num_request:', num_request


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



    ! TODO: remove this test or change it to debug
    do i = 1, num_request
       if(non_local_node_owner(i) > n_cpu -1 ) write(*,*) 'owner > n_cpu:', i, non_local_node_owner(i), non_local_node_keys(i)
    end do
    
    
    ! sort keys accorting to owner    
    call sort(non_local_node_owner(1:num_request), int_arr(1:num_request))

    do i = 1, num_request
       if(non_local_node_owner(i) > n_cpu -1 ) write(*,*) 'after sort, owner > n_cpu:', i, non_local_node_owner(i), non_local_node_keys(i)
    end do


    ! sort keys according to owners
    key_arr_cp(1:num_request) = non_local_node_keys(1:num_request)
    do i= 1, num_request
       non_local_node_keys(i) = key_arr_cp(int_arr(i))
    end do

    ! sort nodes according to owners
    node_arr_cp(1:num_request) = non_local_node_node(1:num_request)
    do i= 1, num_request
       non_local_node_node(i) = node_arr_cp(int_arr(i))
    end do



    requests_per_process = 0

    do i = 1, num_request
       ! use owner + 1, because owner is from 0 and array index from 1
       requests_per_process(non_local_node_owner(i)+1) = requests_per_process(non_local_node_owner(i)+1) + 1
    end do


    if( requests_per_process(my_rank+1) .ne. 0) write (*,*) 'on rank', my_rank, 'requests for self is non-zero:', & 
         requests_per_process(my_rank+1)

!    write(*,*) 'requests 1:', my_rank, requests_per_process

    call MPI_ALLTOALL(requests_per_process, 1, MPI_INTEGER, requests_from_process, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    
    total_num_requests_from_others = sum(requests_from_process)

!    write(*,*) 'requests 2:', my_rank, requests_from_process


    allocate(requested_keys(total_num_requests_from_others))

    disp = 0
    do i = 1, n_cpu
       sdispls(i) = disp
       disp = disp + requests_per_process(i)
    end do
       
    disp = 0
    do i = 1, n_cpu
       rdispls(i) = disp
       disp = disp + requests_from_process(i)
    end do

!    write (*,*) 'sdispls:', my_rank, sdispls
!    write (*,*) 'rdispls:', my_rank, rdispls


    call MPI_ALLTOALLV(non_local_node_keys, requests_per_process, sdispls, MPI_INTEGER8, &
         requested_keys, requests_from_process, rdispls, MPI_INTEGER8, MPI_COMM_WORLD, ierr)

    ! test whether requested keys are locally known
    do i = 1, total_num_requests_from_others
       actual_address = key2addr(requested_keys(i), 'update properties: test requested keys')
    end do

    
    ! register propertyupdate data type
    blocklengths(1:nprops_property_update)  = [1, 1, 1, 1, 3, 1]
    types(1:nprops_property_update)         = [MPI_INTEGER8, MPI_INTEGER, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
    call MPI_GET_ADDRESS( dummy_property_update,                  address(0), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%key,              address(1), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%owner,            address(2), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%smoothing_length, address(3), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%rho,              address(4), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%v,                address(5), ierr )
    call MPI_GET_ADDRESS( dummy_property_update%temperature,      address(6), ierr )
    displacements(1:nprops_property_update) = int(address(1:nprops_property_update) - address(0))
    call MPI_TYPE_STRUCT( nprops_property_update, blocklengths, displacements, types, mpi_type_property_update, ierr )
    call MPI_TYPE_COMMIT( mpi_type_property_update, ierr)


    allocate( packed_updates(total_num_requests_from_others), received_updates(num_request) )

    ! pack everything together
    do i = 1, total_num_requests_from_others
       actual_address = key2addr(requested_keys(i), 'update properties: packaging particles')

       actual_node = htable(actual_address)%node

       packed_updates(i) = t_property_update( htable(actual_address)%key, htable(actual_address)%owner, particles(actual_node)%results%h, &
            particles(actual_node)%results%rho, tree_nodes(actual_node)%v, tree_nodes(actual_node)%temperature )

       ! test whether requested key is parent of particle key
       ! TODO: write a test funciton for this?
       ! write(*,'(i3,a,O30,O30)') my_rank, 'packing:', htable(actual_address)%key, particles(actual_node)%key

    end do
    
    disp = 0
    do i = 1, n_cpu
       sdispls(i) = disp
       disp = disp + requests_from_process(i)
    end do
       
    disp = 0
    do i = 1, n_cpu
       rdispls(i) = disp
       disp = disp + requests_per_process(i)
    end do


    call MPI_ALLTOALLV(packed_updates, requests_from_process, sdispls, mpi_type_property_update, &
         received_updates, requests_per_process, rdispls, mpi_type_property_update, MPI_COMM_WORLD, ierr)


    do i= 1, num_request
       if( received_updates(i)%key .ne. non_local_node_keys(i) ) write(*,*) 'Error in update on', my_rank, 'key mismatch', received_updates(i)%key, non_local_node_keys(i)
    end do

    do i= 1, num_request
       actual_address = key2addr( received_updates(i)%key, 'update properties: updating properties of remote nodes' )
       actual_node = htable( actual_address )%node
       
       if(actual_node .ne. non_local_node_node(i) ) then 
          write(76, *) my_rank, i, received_updates(i)%key, actual_node, non_local_node_node(i)
       end if


       tree_nodes(actual_node)%rho         = received_updates(i)%rho
       tree_nodes(actual_node)%temperature = received_updates(i)%temperature
       tree_nodes(actual_node)%v           = received_updates(i)%v
       tree_nodes(actual_node)%h           = received_updates(i)%smoothing_length
       
    end do


    do i = 1, np_local
       
       tree_nodes(i)%rho         = particles(i)%results%rho
       tree_nodes(i)%temperature = particles(i)%data%temperature
       tree_nodes(i)%h           = particles(i)%results%h

    end do






    ! TODO: update tree_nodes%h for all parents


    deallocate( non_local_node_node, non_local_node_keys, non_local_node_owner, requests_per_process, key_arr_cp, int_arr, node_arr_cp, STAT=ierr )

    
    deallocate( requested_keys, STAT=ierr ) 
    
    deallocate( packed_updates, received_updates, STAT=ierr )
    
    
    
  end subroutine update_particle_props
  


end module module_sph

