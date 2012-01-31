module particle_pusher

  implicit none
  private


    public velocities
    public push

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Calculate velocities from sph force
    !>   without any constraints
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine velocities(p_start,p_finish,delta_t)
      use physvars, only: &
           particles

      use module_interaction_specific_types, only: &
           PARTICLE_TYPE_FIXED
           
      implicit none

      real, intent(in) :: delta_t
      integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

      integer :: p

      ! unconstrained motion by default
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
      do p = p_start, p_finish
         if( btest(particles(p)%data%type, PARTICLE_TYPE_FIXED) ) then
            ! don't move this particle
            particles(p)%data%v          = 0
            particles(p)%data%v_minus_half = 0
         else
            ! for gravity, mass and charge are equal so q * e / m = e
            ! because the velocity at the same timestep as the coordinate is needed for sph ( v(t + dt) ), and v_minus_half is dt/2 behind:
            particles(p)%data%v          = particles(p)%data%v_minus_half + 3._8/2._8 * delta_t * particles(p)%results%sph_force
            ! v(t - dt/2) for leap frog
            particles(p)%data%v_minus_half = particles(p)%data%v_minus_half + delta_t * particles(p)%results%sph_force
         end if
      end do
      !$OMP END PARALLEL DO
      

    end subroutine velocities


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Update particle positions - used with leap-frog scheme
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push(ips,ipf,delt)

      use physvars, only: &
           particles
      
      implicit none
      
      integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
      real, intent(in) :: delt
      integer :: p

      
      !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
      do p=ips,ipf
         particles(p)%x = particles(p)%x + particles(p)%data%v_minus_half * delt
      end do
      !$OMP END PARALLEL DO

    end subroutine push

end module particle_pusher




