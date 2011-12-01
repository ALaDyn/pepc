!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is directly involved in force calculation
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_calc_force

     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8 !< TODO: currently unused
      real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear
      private calc_2nd_algebraic_condensed
      !private calc_2nd_algebraic_decomposed
      !private calc_6th_algebraic_decomposed
      private G_core
      private G_decomp


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
        !> generic Multipole Acceptance Criterion
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function mac(particle, node, cf_par, dist2, boxlength2)
            use treetypes
            implicit none

            logical :: mac
            integer, intent(in) :: node
            type(t_calc_force_params), intent(in) :: cf_par
            real*8, intent(in) :: dist2
            real*8, intent(in) :: boxlength2
            type(t_particle), intent(in) :: particle

            select case (cf_par%mac)
                case (0)
                    ! Barnes-Hut-MAC
                    mac = (cf_par%theta2 * dist2 > boxlength2)
                case default
                    ! N^2 code
                    mac = .false.
            end select

        end function

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> clears result in t_particle datatype - usually, this function does not need to be touched
        !> due to dependency on treetypes and(!) on module_interaction_specific, the
        !> function cannot reside in module_interaction_specific that may not include treetypes
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particleresults_clear(particles, nparticles)
          use treetypes
          implicit none
          type(t_particle), intent(inout) :: particles(nparticles)
          integer, intent(in) :: nparticles

          particles(1:nparticles)%results = EMPTY_PARTICLE_RESULTS

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, inode, delta, dist2, vbox, cf_par)
          use treetypes
          implicit none

          integer, intent(in) :: inode
          type(t_particle), intent(inout) :: particle
          real*8, intent(in) :: vbox(3), delta(3), dist2
          !> Force law struct has following content (defined in module treetypes)
          !> These need to be included/defined in call to fields from frontend
          !>    real    :: eps
          !>    real    :: force_const
          !>    integer :: force_law   0= no interaction (default); 2=2nd order condensed algebraic kernel; 3=2nd order decomposed algebraic kernel
          type(t_calc_force_params), intent(in) :: cf_par

          real*8 :: u(3), af(3)

          select case (cf_par%force_law)
            case (2)  !  use 2nd order algebraic kernel, condensed
                call calc_2nd_algebraic_condensed(particle, inode, delta, dist2, cf_par, u, af)

            case (3)  !  TODO: use 2nd order algebraic kernel, decomposed
                !call calc_2nd_algebraic_decomposed(particle, inode, delta, dist2, cf_par, u, af)
                u = 0.
                af = 0.

            case (4)  !  TODO: use 6th order algebraic kernel, decomposed
                !call calc_6th_algebraic_decomposed(inode, delta, dist2, cf_par, u, af)
                u = 0.
                af = 0.

            case default
              u = 0.
              af = 0.
          end select

          ! TODO: factor out multiplication of force_const, does not depend on actual interaction-pair
          particle%results%u(1:3)    = particle%results%u(1:3)     - cf_par%force_const * u(1:3)
          particle%results%af(1:3)   = particle%results%af(1:3)    + cf_par%force_const * af(1:3)

          particle%work = particle%work + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle (not required in vortex bubu, yet)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles, cf_par)
          use treetypes
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)
          type(t_calc_force_params), intent(in) :: cf_par

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D 2nd order condensed algebraic kernel interaction
        !> of particle p with tree node inode, results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine calc_2nd_algebraic_condensed(particle, inode, d, dist2, cf_par, u, af)
            use treetypes
            use treevars
            use module_interaction_specific
            implicit none

            type(t_particle), intent(in) :: particle
            integer, intent(in) :: inode !< index of particle to interact with
            real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
            type(t_calc_force_params), intent(in) :: cf_par !< Force parameters - see module_treetypes
            real*8, intent(out) ::  u(1:3), af(1:3)

            type(t_multipole_data), pointer :: t

            integer :: i1, i2, i3 !< helper variables for the tensor structures

            real*8 :: dx, dy, dz, sig2 !< temp variables for distance and smoothing radius
            real*8 :: Gc25,Gc35,Gc45,Gc55,MPa1,DPa1,DPa2,QPa1,QPa2 !< prefactors for the multipole expansion
            real*8, dimension(3) :: vort !< temp variables for vorticity (or better: alpha)

            ! tensors allow nice and short code and better comparison with (my!) theory
            real*8, dimension(3) :: m0, CP0 !< data structures for the monopole moments
            real*8, dimension(3,3) :: m1, CP1 !< data structures for the dipole moments
            real*8, dimension(3,3,3) :: m2, CP2 !< data structures for the quadrupole moments

            t=>tree_nodes(inode)

            dx = d(1)
            dy = d(2)
            dz = d(3)

            sig2 = cf_par%eps**2

            vort = [particle%data%alpha(1),particle%data%alpha(2),particle%data%alpha(3)]  ! need particle's vorticity for cross-product here

            m0 = [t%chargex,t%chargex,t%chargex]       ! monopole moment tensor
            CP0 = cross_prod(m0,vort)                  ! cross-product for 1st expansion term

            m1 = reshape([t%xdip1,t%xdip2,t%xdip3, &    ! dipole moment tensor
                          t%ydip1,t%ydip2,t%ydip3, &
                          t%zdip1,t%zdip2,t%zdip3],[3,3])
            CP1 = reshape([cross_prod(m1(:,1),vort), &                ! cross-product for 2nd expansion term
                           cross_prod(m1(:,2),vort), &
                           cross_prod(m1(:,3),vort)],[3,3])

            m2 = reshape([t%xxquad1,t%xxquad2,t%xxquad3, &                                                  ! quadrupole moment tensor
                          t%xyquad1,t%xyquad2,t%xyquad3, &
                          t%xzquad1,t%xzquad2,t%xzquad3, &
                          t%xyquad1,t%xyquad2,t%xyquad3, &
                          t%yyquad1,t%yyquad2,t%yyquad3, &
                          t%yzquad1,t%yzquad2,t%yzquad3, &
                          t%xzquad1,t%xzquad2,t%xzquad3, &
                          t%yzquad1,t%yzquad2,t%yzquad3, &
                          t%zzquad1,t%zzquad2,t%zzquad3],[3,3,3])
            CP2 = reshape([cross_prod(m2(:,1,1),vort),cross_prod(m2(:,2,1),vort),cross_prod(m2(:,3,1),vort), &          ! cross-product for 3rd expansion term
                           cross_prod(m2(:,1,2),vort),cross_prod(m2(:,2,2),vort),cross_prod(m2(:,3,2),vort), &
                           cross_prod(m2(:,1,3),vort),cross_prod(m2(:,2,3),vort),cross_prod(m2(:,3,3),vort)],[3,3,3])

            ! precompute kernel function evaluations of various order
            Gc25 = G_core(dist2,sig2,2.5D00)
            Gc35 = G_core(dist2,sig2,3.5D00)
            Gc45 = G_core(dist2,sig2,4.5D00)
            Gc55 = G_core(dist2,sig2,5.5D00)

            MPa1 = 3.0D00*Gc35*dot_product(d,CP0)   ! monopole prefactor for af
            DPa1 = 15.0D00*G_core(dist2,sig2,4.5D00)*sum( (/ (sum( (/ (CP1(i2,i1)*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )  ! dipole prefators for af
            DPa2 = (CP1(1,1)+CP1(2,2)+CP1(3,3))
            QPa1 = 52.5D00*G_core(dist2,sig2,5.5D00)*sum( (/ (sum( (/ (sum( (/ (CP2(i3,i2,i1)*d(i3),i3=1,3) /) )*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) ! quadrupole prefactors for af
            QPa2 = dot_product(d,CP2(1,1,:)+CP2(2,2,:)+CP2(3,3,:)+CP2(1,:,1)+CP2(2,:,2)+CP2(3,:,3)+CP2(:,1,1)+CP2(:,2,2)+CP2(:,3,3))

            u(1) = Gc25* (dy*m0(3)-dz*m0(2)) &                                                                                       ! MONOPOLE

                    + 3.0D00*Gc35* sum( (/ ((m1(3,i1)*dy-m1(2,i1)*dz)*d(i1),i1=1,3) /) ) - Gc25* (m1(3,2)-m1(2,3)) &                 ! DIPOLE

                    - 1.5D00*Gc35* ( sum( (/ (m2(3,i1,i1)*dy-m2(2,i1,i1)*dz,i1=1,3) /) ) + 2.0*sum( (/ ((m2(3,i1,2)-m2(2,i1,3))*d(i1),i1=1,3) /) ) ) + &
                      7.5D00*Gc45* sum( (/ (sum( (/ ((dy*m2(3,i2,i1)-dz*m2(2,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            u(2) = Gc25* (dz*m0(1)-dx*m0(3)) &                                                                                       ! MONOPOLE

                    + 3.0D00*Gc35* sum( (/ ((m1(1,i1)*dz-m1(3,i1)*dx)*d(i1),i1=1,3) /) ) - Gc25* (m1(1,3)-m1(3,1)) &                 ! DIPOLE

                    - 1.5D00*Gc35* ( sum( (/ (m2(1,i1,i1)*dz-m2(3,i1,i1)*dx,i1=1,3) /) ) + 2.0*sum( (/ ((m2(1,i1,3)-m2(3,i1,1))*d(i1),i1=1,3) /) ) ) + &
                      7.5D00*Gc45* sum( (/ (sum( (/ ((dz*m2(1,i2,i1)-dx*m2(3,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            u(3) = Gc25* (dx*m0(2)-dy*m0(1)) &                                                                                       ! MONOPOLE

                    + 3.0D00*Gc35* sum( (/ ((m1(2,i1)*dx-m1(1,i1)*dy)*d(i1),i1=1,3) /) ) - Gc25* (m1(2,1)-m1(1,2)) &                 ! DIPOLE

                    - 1.5D00*Gc35* ( sum( (/ (m2(2,i1,i1)*dx-m2(1,i1,i1)*dy,i1=1,3) /) ) + 2.0*sum( (/ ((m2(2,i1,1)-m2(1,i1,2))*d(i1),i1=1,3) /) ) ) + &
                      7.5D00*Gc45* sum( (/ (sum( (/ ((dx*m2(2,i2,i1)-dy*m2(1,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            af(1) = Mpa1*dx - Gc25*CP0(1)  &                                                                                               ! MONOPOLE

                    + DPa1*dx - 3.0D00*Gc35* ( dot_product(CP1(:,1)+CP1(1,:),d) + dx*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dx - &
                        7.5D00*Gc45*( sum( (/ (sum( (/ ((CP2(i2,i1,1)+CP2(i2,1,i1)+CP2(1,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dx*QPa2) + &
                        1.5D00*Gc35*( CP2(1,1,1)+CP2(2,2,1)+CP2(3,3,1)+CP2(1,1,1)+CP2(2,1,2)+CP2(3,1,3)+CP2(1,1,1)+CP2(1,2,2)+CP2(1,3,3) ) ! QUADRUPOLE


            af(2) = Mpa1*dy - Gc25*CP0(2)  &                                                                                               ! MONOPOLE

                    + DPa1*dy - 3.0D00*Gc35* ( dot_product(CP1(:,2)+CP1(2,:),d) + dy*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dy - &
                        7.5D00*Gc45*( sum( (/ (sum( (/ ((CP2(i2,i1,2)+CP2(i2,2,i1)+CP2(2,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dy*QPa2) + &
                        1.5D00*Gc35*( CP2(1,1,2)+CP2(2,2,2)+CP2(3,3,2)+CP2(1,2,1)+CP2(2,2,2)+CP2(3,2,3)+CP2(2,1,1)+CP2(2,2,2)+CP2(2,3,3) ) ! QUADRUPOLE


            af(3) = Mpa1*dz - Gc25*CP0(3)  &                                                                                               ! MONOPOLE

                    + DPa1*dz - 3.0D00*Gc35* ( dot_product(CP1(:,3)+CP1(3,:),d) + dz*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dz - &
                        7.5D00*Gc45*( sum( (/ (sum( (/ ((CP2(i2,i1,3)+CP2(i2,3,i1)+CP2(3,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dz*QPa2) + &
                        1.5D00*Gc35*( CP2(1,1,3)+CP2(2,2,3)+CP2(3,3,3)+CP2(1,3,1)+CP2(2,3,2)+CP2(3,3,3)+CP2(3,1,1)+CP2(3,2,2)+CP2(3,3,3) ) ! QUADRUPOLE

        end subroutine calc_2nd_algebraic_condensed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D 2nd order decomposed algebraic kernel interaction
        !> of particle p with tree node inode that is shifted by the lattice
        !> vector vbox results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        subroutine calc_2nd_algebraic_decomposed(inode, d, dist2, cf_par, u, af)
!
!        end subroutine calc_2nd_algebraic_decomposed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Helper functions for multipole expansion
        !> of particle p with tree node inode that is shifted by the lattice
        !> vector vbox results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        function G_core(r,s,factor)
           implicit none

           real*8 :: G_core
           real*8, intent(in) :: r,s,factor

           G_core = (r+factor*s)/((r+s)**factor)

        end function


        function G_decomp(r,s,tau)
           implicit none

           real*8 :: G_decomp
           real*8, intent(in) :: r,s,tau

           G_decomp = 1.0/((r+s)**tau)

        end function


        function cross_prod(vec_a, vec_b)
            implicit none

            real*8, dimension(3) :: cross_prod
            real*8, dimension(3), intent(in) :: vec_a, vec_b

            cross_prod(1) = vec_a(2)*vec_b(3) - vec_a(3)*vec_b(2)
            cross_prod(2) = vec_a(3)*vec_b(1) - vec_a(1)*vec_b(3)
            cross_prod(3) = vec_a(1)*vec_b(2) - vec_a(2)*vec_b(1)

        end function

  end module module_calc_force
