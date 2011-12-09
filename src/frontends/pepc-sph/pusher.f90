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
           
      implicit none

      real, intent(in) :: delta_t
      integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

      integer :: p

      ! unconstrained motion by default
      do p = p_start, p_finish
         ! for gravity, mass and charge are equal so q * e / m = e
         ! because the velocity at the same timestep as the coordinate is needed for sph ( v(t + dt) ):
         particles(p)%data%v          = particles(p)%data%v_and_half + delta_t * particles(p)%results%sph_force / 2._8
         ! v(t + dt *(1 + 1/2)) for leap frog
         particles(p)%data%v_and_half = particles(p)%data%v_and_half + delta_t * particles(p)%results%sph_force
      end do

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

      do p=ips,ipf
         particles(p)%x = particles(p)%x + particles(p)%data%v_and_half * delt
      end do


    end subroutine push

end module particle_pusher




