!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Interaction specific wrappers to support frontends with old interfaces
!> nothing is obligatory for the treecode here
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_pepc_wrappers
     use treetypes, only : t_particle, t_calc_force_params
     use module_interaction_specific, only : t_particle_results, t_particle_data
     implicit none
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      type(t_particle), public, dimension(:), allocatable :: particles

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public pepc_fields_coulomb_wrapper
      public pepc_grid_fields_coulomb_wrapper

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
    !>   Calculate fields and potential for supplied particle coordinates p_x, p_y, p_z and charges p_q
    !>
    !>   Returns fields Ex, Ey, Ez and potential pot excluding external terms
    !>   @param[in] np_local local number of particles
    !>   @param[in] npart_total total particle number
    !>   @param[in] p_x dimension(1:np_local) - x-component of particle coordinates
    !>   @param[in] p_y dimension(1:np_local) - y-component of particle coordinates
    !>   @param[in] p_z dimension(1:np_local) - z-component of particle coordinates
    !>   @param[in] p_q dimension(1:np_local) - particle charge
    !>   @param[inout] p_w dimension(1:np_local) - particle workload from previous iteration (should be set to 1.0 for the first timestep)
    !>   @param[in] p_label dimension(1:np_local) - particle label (may any number except zero)
    !>   @param[out] p_Ex dimension(1:np_local) - x-component of electric field
    !>   @param[out] p_Ey dimension(1:np_local) - y-component of electric field
    !>   @param[out] p_Ez dimension(1:np_local) - z-component of electric field
    !>   @param[out] p_pot dimension(1:np_local) - electric potential
    !>   @param[in] np_mult_ memory allocation parameter
    !>   @param[in] cf_par parameters for force summation
    !>   @param[in] itime current simulation timestep number
    !>   @param[in] weighted selector for load balancing
    !>   @param[in] curve_type selector for type of space filling curve
    !>   @param[in] num_neighbours number of neighbour boxes to be considered during tree walk
    !>   @param[in] neighbours shift vectors to neighbour boxes
    !>   @param[in] no_dealloc if set to .true., deallocation of tree-structures is prevented to allow for front-end triggered diagnostics
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_fields_coulomb_wrapper(np_local,npart_total,p_x, p_y, p_z, p_q, p_w, p_label, &
                    p_Ex, p_Ey, p_Ez, p_pot, np_mult_, cf_par, itime, weighted, curve_type, &
                    num_neighbours, neighbours, no_dealloc, no_restore)
        use treevars
        use module_pepcfields
        implicit none
        integer, intent(inout) :: np_local  ! # particles on this CPU
        integer, intent(in) :: npart_total ! total # simulation particles
        real, intent(in) :: np_mult_       ! multipole opening angle
        type(t_calc_force_params), intent(in) :: cf_par
        integer, intent(in) :: itime  ! timestep
        integer, intent(in) :: weighted
        real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc
        real*8, intent(in), dimension(np_local) :: p_q ! charges, masses
        integer, intent(in), dimension(np_local) :: p_label  ! particle label
        real*8, intent(out), dimension(np_local) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
        integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
        integer, intent(in) :: neighbours(3, num_neighbours) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
        real*8, dimension(np_local) :: p_w ! work loads
        integer, intent(in) :: curve_type ! type of space-filling curve
        logical, intent(in) :: no_dealloc, no_restore

        integer :: i

        type(t_particle_results), dimension(:), allocatable :: particle_results

        if (allocated(particles))        deallocate(particles)
        if (allocated(particle_results)) deallocate(particle_results)

        allocate(particles(1:np_local))
        allocate(particle_results(1:np_local))

        if (force_debug) then
            write (*,'(a7,a50/2i5,4f15.2)') 'PEPC | ','Params: itime, mac, theta, eps, force_const:', &
            itime, cf_par%mac, cf_par%theta, cf_par%eps, cf_par%force_const
            write (*,'(a7,a20/(i16,4f15.3,i8))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_q(i), &
            p_label(i),i=1,npp)
        endif

        do i=1,np_local
            particles(i) = t_particle( [p_x(i), p_y(i), p_z(i)],       &  ! position
                                              max(p_w(i), 1._8),       &  ! workload from last step
                                                           -1_8,       &  ! key - will be assigned later
                                                     p_label(i),       &  ! particle label for tracking purposes
                                                             me,       &  ! particle owner
                                         t_particle_data( p_q(i) )  )     ! charge etc
        end do

        call pepc_fields(np_local, npart_total, particles, particle_results, &
                             np_mult_, cf_par, itime, weighted, curve_type, num_neighbours, neighbours, no_dealloc, no_restore)

        ! read data from particle_coordinates, particle_results, particle_properties
        do i=1,np_local
          p_ex(i)  = particle_results(i)%e(1)
          p_ey(i)  = particle_results(i)%e(2)
          p_ez(i)  = particle_results(i)%e(3)
          p_pot(i) = particle_results(i)%pot
          p_w(i)   =  particles(i)%work
        end do

        if (force_debug) then
            write (ipefile,'("Tree forces:"/"   p    q   m   pot  ",f8.2)')
            write (*,'("Tree forces:"/"   p    q   m   ux   pot  ",f8.2)')
            write (ipefile,'("Tree forces:"/"   p    q   m   pot  ",f8.2)') cf_par%force_const

            do i=1,np_local
                write (ipefile,'(1x,i7,3(1pe14.5))') particles(i)%label, particles(i)%data%q, p_pot(i), p_ex(i)
                write (*,'(1x,i7,3(1pe14.5))') particles(i)%label, particles(i)%x(1), particles(i)%data%q, p_pot(i)
            end do

        endif

    end subroutine



    subroutine pepc_grid_fields_coulomb_wrapper(ngp,p_x, p_y, p_z, p_label, p_Ex, p_Ey, p_Ez, p_pot, &
                              cf_par, num_neighbour_boxes, neighbour_boxes)
      use treevars, only : me
      use module_pepcfields
      implicit none
      integer, intent(in) :: ngp
      real*8, intent(in) :: p_x(ngp), p_y(ngp), p_z(ngp)
      integer, intent(in) :: p_label(ngp)
      real*8, intent(out) :: p_Ex(ngp), p_Ey(ngp), p_Ez(ngp), p_pot(ngp)
      type(t_calc_force_params), intent(in) :: cf_par
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

      type(t_particle),         dimension(:), allocatable :: grid_particles
      type(t_particle_results), dimension(:), allocatable :: grid_particle_results

      integer :: i

      allocate(grid_particles(ngp), grid_particle_results(ngp))

      do i=1,ngp
        grid_particles(i) = t_particle( [p_x(i), p_y(i), p_z(i)],       &  ! position
                                                              1.,       &  ! workload from last step
                                                            -1_8,       &  ! key - will be assigned later
                                                      p_label(i),       &  ! particle label for tracking purposes
                                                              me,       &  ! particle owner
                                          t_particle_data( 0.0 )  )        ! charge etc
      end do


      call pepc_grid_fields(ngp, grid_particles, grid_particle_results, cf_par, num_neighbour_boxes, neighbour_boxes)

        do i=1,ngp
          p_ex(i)  = grid_particle_results(i)%e(1)
          p_ey(i)  = grid_particle_results(i)%e(2)
          p_ez(i)  = grid_particle_results(i)%e(3)
          p_pot(i) = grid_particle_results(i)%pot
        end do

      deallocate(grid_particles, grid_particle_results)

    end subroutine


end module module_pepc_wrappers
