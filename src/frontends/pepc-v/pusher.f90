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
         ux(p) = ux(p) + delta_t * q(p)*ex(p)/m(p)
         uy(p) = uy(p) + delta_t * q(p)*ey(p)/m(p)
         uz(p) = uz(p) + delta_t * q(p)*ez(p)/m(p)
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
         x(p)=x(p)+ux(p)*delt
         y(p)=y(p)+uy(p)*delt
         z(p)=z(p)+uz(p)*delt
      end do

    end subroutine push

end module particle_pusher




