module particle_pusher
  implicit none
  private


    public velocities
    public push

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Calculate velocities from fields(dorces)
    !>   without any constraints
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine velocities(p_start,p_finish,delta_t)
      use physvars
      implicit none

      real, intent(in) :: delta_t
      integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

      integer :: p

      !  Available ensemble modes
      !      1 = NVE - total energy conserved

      ! unconstrained motion by default (scheme=1)
      do p = p_start, p_finish
         u(1:3,p) = u(1:3,p) + delta_t * particles(p)%data%q * particle_Results(p)%e(1:3) / m(p)
      end do

    end subroutine velocities


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Update particle positions - used with leap-frog scheme
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push(ips,ipf,delt)

      use physvars
      integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
      real, intent(in) :: delt
      integer :: p

      do p=ips,ipf
         particles(p)%x(1:3) = particles(p)%x(1:3) + u(1:3,p) * delt
      end do

    end subroutine push

end module particle_pusher




