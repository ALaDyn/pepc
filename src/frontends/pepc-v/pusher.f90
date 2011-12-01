module particle_pusher
  implicit none
  private


    public push_rk2

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Update particle positions - using 2nd order RK
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_rk2(stage)

      use physvars
      use module_remesh
      integer, intent(in) :: stage   ! In which RK stage are we?
      integer :: i

      do i=1,np

        if (stage == 1) then
            ! Euler predictor
            vortex_particles(i)%data%u_rk(1:3)  = vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%af_rk(1:3) = vortex_particles(i)%results%af(1:3)
            vortex_particles(i)%x(1:3)     = vortex_particles(i)%x(1:3) + dt*vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha(1:3) + dt*vortex_particles(i)%results%af(1:3)
        else
            ! Trapezoidal corrector
            vortex_particles(i)%x(1:3) = vortex_particles(i)%x(1:3)-dt*vortex_particles(i)%data%u_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%u_rk(1:3)+vortex_particles(i)%results%u(1:3))
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha(1:3)-dt*vortex_particles(i)%data%af_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%af_rk(1:3)+vortex_particles(i)%results%af(1:3))
        end if

      end do

      call kick_out_particles()

    end subroutine push_rk2

end module particle_pusher




