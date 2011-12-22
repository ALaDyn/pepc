!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_interaction_specific
     use module_pepc_types
     use module_interaction_specific_types
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

      integer, public :: force_law    = 2      !< 2 = calc_2nd_algebraic_condensed
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      real*8, public  :: theta2       = 0.6**2.  !< square of multipole opening angle
      real*8, public  :: sig2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)

      namelist /calc_force_vortex/ force_law, mac_select, theta2, sig2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public multipole_from_particle
      public shift_multipoles_up
      public results_add
      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear
      public calc_force_init
      public calc_force_finalize
      public calc_force_prepare
      public get_number_of_interactions_per_particle
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
      !> Computes multipole properties of a single particle
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine multipole_from_particle(particle_pos, particle, multipole)
        implicit none
        real*8, intent(in) :: particle_pos(3)
        type(t_particle_data), intent(in) :: particle
        type(t_tree_node_interaction_data), intent(out) :: multipole

        multipole = t_tree_node_interaction_data(particle_pos, sqrt(dot_product(particle%alpha, particle%alpha)), particle%alpha(1), particle%alpha(2), particle%alpha(3), &
                                     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> Accumulates multipole properties of child nodes to parent node
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine shift_multipoles_up(parent, children)
        implicit none
        type(t_tree_node_interaction_data), intent(out) :: parent
        type(t_tree_node_interaction_data), intent(in) :: children(:)

        integer :: nchild, j

        real*8 :: shift(1:3)

        nchild = size(children)

        parent%chargex     = SUM( children(1:nchild)%chargex )
        parent%chargey     = SUM( children(1:nchild)%chargey )
        parent%chargez     = SUM( children(1:nchild)%chargez )
        parent%abs_charge  = SUM( children(1:nchild)%abs_charge )

        ! centre of charge
        parent%coc        = [0., 0., 0.]
        do j=1,nchild
          parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%abs_charge )
        end do
        parent%coc(1:3) = parent%coc(1:3) / parent%abs_charge

        ! multipole properties
        parent%xdip1    = 0.
        parent%ydip1    = 0.
        parent%zdip1    = 0.
        parent%xdip2    = 0.
        parent%ydip2    = 0.
        parent%zdip2    = 0.
        parent%xdip3    = 0.
        parent%ydip3    = 0.
        parent%zdip3    = 0.
        parent%xxquad1  = 0.
        parent%xyquad1  = 0.
        parent%xzquad1  = 0.
        parent%yzquad1  = 0.
        parent%yyquad1  = 0.
        parent%zzquad1  = 0.
        parent%xxquad2  = 0.
        parent%xyquad2  = 0.
        parent%xzquad2  = 0.
        parent%yzquad2  = 0.
        parent%yyquad2  = 0.
        parent%zzquad2  = 0.
        parent%xxquad3  = 0.
        parent%xyquad3  = 0.
        parent%xzquad3  = 0.
        parent%yzquad3  = 0.
        parent%yyquad3  = 0.
        parent%zzquad3  = 0.

        do j=1,nchild
          shift(1:3) = parent%coc(1:3) - children(j)%coc

          ! dipole moment
          parent%xdip1 = parent%xdip1 + children(j)%xdip1 - children(j)%chargex*shift(1)
          parent%ydip1 = parent%ydip1 + children(j)%ydip1 - children(j)%chargex*shift(2)
          parent%zdip1 = parent%zdip1 + children(j)%zdip1 - children(j)%chargex*shift(3)

          parent%xdip2 = parent%xdip2 + children(j)%xdip2 - children(j)%chargey*shift(1)
          parent%ydip2 = parent%ydip2 + children(j)%ydip2 - children(j)%chargey*shift(2)
          parent%zdip2 = parent%zdip2 + children(j)%zdip2 - children(j)%chargey*shift(3)

          parent%xdip3 = parent%xdip3 + children(j)%xdip3 - children(j)%chargez*shift(1)
          parent%ydip3 = parent%ydip3 + children(j)%ydip3 - children(j)%chargez*shift(2)
          parent%zdip3 = parent%zdip3 + children(j)%zdip3 - children(j)%chargez*shift(3)

          ! quadrupole moment
          parent%xxquad1 = parent%xxquad1 + children(j)%xxquad1 - children(j)%xdip1*shift(1) - children(j)%xdip1*shift(1) + children(j)%chargex*shift(1)*shift(1)
          parent%xyquad1 = parent%xyquad1 + children(j)%xyquad1 - children(j)%xdip1*shift(2) - children(j)%ydip1*shift(1) + children(j)%chargex*shift(1)*shift(2)
          parent%xzquad1 = parent%xzquad1 + children(j)%xzquad1 - children(j)%xdip1*shift(3) - children(j)%zdip1*shift(1) + children(j)%chargex*shift(1)*shift(3)
          parent%yzquad1 = parent%yzquad1 + children(j)%yzquad1 - children(j)%ydip1*shift(3) - children(j)%zdip1*shift(2) + children(j)%chargex*shift(2)*shift(3)
          parent%yyquad1 = parent%yyquad1 + children(j)%yyquad1 - children(j)%ydip1*shift(2) - children(j)%ydip1*shift(2) + children(j)%chargex*shift(2)*shift(2)
          parent%zzquad1 = parent%zzquad1 + children(j)%zzquad1 - children(j)%zdip1*shift(3) - children(j)%zdip1*shift(3) + children(j)%chargex*shift(3)*shift(3)

          parent%xxquad2 = parent%xxquad2 + children(j)%xxquad2 - children(j)%xdip2*shift(1) - children(j)%xdip2*shift(1) + children(j)%chargey*shift(1)*shift(1)
          parent%xyquad2 = parent%xyquad2 + children(j)%xyquad2 - children(j)%xdip2*shift(2) - children(j)%ydip2*shift(1) + children(j)%chargey*shift(1)*shift(2)
          parent%xzquad2 = parent%xzquad2 + children(j)%xzquad2 - children(j)%xdip2*shift(3) - children(j)%zdip2*shift(1) + children(j)%chargey*shift(1)*shift(3)
          parent%yzquad2 = parent%yzquad2 + children(j)%yzquad2 - children(j)%ydip2*shift(3) - children(j)%zdip2*shift(2) + children(j)%chargey*shift(2)*shift(3)
          parent%yyquad2 = parent%yyquad2 + children(j)%yyquad2 - children(j)%ydip2*shift(2) - children(j)%ydip2*shift(2) + children(j)%chargey*shift(2)*shift(2)
          parent%zzquad2 = parent%zzquad2 + children(j)%zzquad2 - children(j)%zdip2*shift(3) - children(j)%zdip2*shift(3) + children(j)%chargey*shift(3)*shift(3)

          parent%xxquad3 = parent%xxquad3 + children(j)%xxquad3 - children(j)%xdip3*shift(1) - children(j)%xdip3*shift(1) + children(j)%chargez*shift(1)*shift(1)
          parent%xyquad3 = parent%xyquad3 + children(j)%xyquad3 - children(j)%xdip3*shift(2) - children(j)%ydip3*shift(1) + children(j)%chargez*shift(1)*shift(2)
          parent%xzquad3 = parent%xzquad3 + children(j)%xzquad3 - children(j)%xdip3*shift(3) - children(j)%zdip3*shift(1) + children(j)%chargez*shift(1)*shift(3)
          parent%yzquad3 = parent%yzquad3 + children(j)%yzquad3 - children(j)%ydip3*shift(3) - children(j)%zdip3*shift(2) + children(j)%chargez*shift(2)*shift(3)
          parent%yyquad3 = parent%yyquad3 + children(j)%yyquad3 - children(j)%ydip3*shift(2) - children(j)%ydip3*shift(2) + children(j)%chargez*shift(2)*shift(2)
          parent%zzquad3 = parent%zzquad3 + children(j)%zzquad3 - children(j)%zdip3*shift(3) - children(j)%zdip3*shift(3) + children(j)%chargez*shift(3)*shift(3)

        end do

        parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2+(parent%coc(3)-children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> adds res2 to res1
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine results_add(res1, res2)
        implicit none
        type(t_particle_results), intent(inout) :: res1
        type(t_particle_results), intent(in) :: res2

        res1%u    = res1%u    + res2%u
        res1%af   = res1%af   + res2%af
      end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> initializes interaction specific parameters, redas them from file
        !> if optional argument para_file_name is given
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_init(para_file_available, para_file_name, my_rank)
          implicit none
          logical, intent(in) :: para_file_available
          character(*), intent(in) :: para_file_name
          integer, intent(in) :: my_rank
          integer, parameter :: para_file_id = 47

          if (para_file_available) then
            open(para_file_id,file=para_file_name)

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section calc_force_vortex: ", para_file_name
            read(para_file_id,NML=calc_force_vortex)

            close(para_file_id)
          endif

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> computes derived parameters for calc force module
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_prepare()
          implicit none
          ! nothing to do here
        end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> subroutine must return the estimated number of iteractions per particle
      !> for the current mac and/or parameters and the supplied total number of particles
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
        implicit none
        integer, intent(in) :: npart_total !< total number of particles
        integer, intent(out) :: nintmax !< maximum number of interactions per particle

        real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

        ! Estimate of interaction list length - Hernquist expression
        ! applies for BH-MAC
        invnintmax = max(theta2 / (35.*log(1.*npart_total)) , 1._8/npart_total)
        nintmax    = int(1._8/invnintmax)

      end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> finalizes the calc force module at end of simulation
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_finalize()
          implicit none
          ! nothing to do here
        end subroutine calc_force_finalize


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> generic Multipole Acceptance Criterion
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function mac(particle, node, dist2, boxlength2)
            use module_pepc_types
            use treevars, only : tree_nodes
            implicit none

            logical :: mac
            integer, intent(in) :: node
            real*8, intent(in) :: dist2
            real*8, intent(in) :: boxlength2
            type(t_particle), intent(in) :: particle

            select case (mac_select)
                case (0)
                    ! Barnes-Hut-MAC
                    mac = (theta2 * dist2 > boxlength2)
                case (1)
                    ! Bmax-MAC
                    mac = (theta2 * dist2 > min(tree_nodes(node)%bmax**2,3.0*boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
                case default
                    ! N^2 code
                    mac = .false.
            end select

        end function

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> clears result in t_particle datatype - usually, this function does not need to be touched
        !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
        !> function cannot reside in module_interaction_specific that may not include module_pepc_types
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particleresults_clear(particles, nparticles)
          use module_pepc_types
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
        subroutine calc_force_per_interaction(particle, inode, delta, dist2, vbox)
          use module_pepc_types
          implicit none

          integer, intent(in) :: inode
          type(t_particle), intent(inout) :: particle
          real*8, intent(in) :: vbox(3), delta(3), dist2
          !> Force law struct has following content (defined in module module_pepc_types)
          !> These need to be included/defined in call to fields from frontend
          !>    real    :: eps
          !>    real    :: force_const
          !>    integer :: force_law   0= no interaction (default); 2=2nd order condensed algebraic kernel; 3=2nd order decomposed algebraic kernel

          real*8 :: u(3), af(3)

          select case (force_law)
            case (2)  !  use 2nd order algebraic kernel, condensed

                if (inode > 0) then
                    ! It's a leaf, do direct summation without ME stuff
                    call calc_2nd_algebraic_condensed_direct(particle, inode, delta, dist2, u, af)
                else
                    ! It's not a leaf, do ME
                    call calc_2nd_algebraic_condensed(particle, inode, delta, dist2, u, af)
                end if

            case (3)  !  TODO: use 2nd order algebraic kernel, decomposed
                !call calc_2nd_algebraic_decomposed(particle, inode, delta, dist2, u, af)
                u = 0.
                af = 0.

            case (6)  ! use 6th order algebraic kernel, condensed


                if (inode > 0) then
                    ! It's a leaf, do direct summation without ME stuff
                    call calc_6th_algebraic_condensed_direct(particle, inode, delta, dist2, u, af) !TODO: 6xth order direct summation
                else
                    ! It's not a leaf, do ME
                    call calc_6th_algebraic_condensed(particle, inode, delta, dist2, u, af)
                end if


            case default
              u = 0.
              af = 0.
          end select

          particle%results%u(1:3)    = particle%results%u(1:3)     -  u(1:3)
          particle%results%af(1:3)   = particle%results%af(1:3)    + af(1:3)
          particle%work = particle%work + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle (not required in vortex bubu, yet)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles)
          use module_pepc_types
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D 2nd order condensed algebraic kernel interaction
        !> of particle p with tree node inode, results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine calc_2nd_algebraic_condensed(particle, inode, d, dist2, u, af)
            use module_pepc_types
            use treevars
            use module_interaction_specific_types
            implicit none

            type(t_particle), intent(in) :: particle
            integer, intent(in) :: inode !< index of particle to interact with
            real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
            real*8, intent(out) ::  u(1:3), af(1:3)

            type(t_tree_node_interaction_data), pointer :: t

            integer :: i1, i2, i3 !< helper variables for the tensor structures

            real*8 :: dx, dy, dz !< temp variables for distance
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

            vort = [particle%data%alpha(1),particle%data%alpha(2),particle%data%alpha(3)]  ! need particle's vorticity for cross-product here

            m0 = [t%chargex,t%chargey,t%chargez]       ! monopole moment tensor
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
        !> Calculates 3D 2nd order condensed algebraic kernel interaction
        !> of particle p with tree node inode, results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_6th_algebraic_condensed(particle, inode, d, dist2, u, af)
            use module_pepc_types
            use treevars
            use module_interaction_specific_types
            implicit none

            type(t_particle), intent(in) :: particle
            integer, intent(in) :: inode !< index of particle to interact with
            real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
            real*8, intent(out) ::  u(1:3), af(1:3)

            type(t_tree_node_interaction_data), pointer :: t

            integer :: i1, i2, i3 !< helper variables for the tensor structures

            real*8 :: dx, dy, dz !< temp variables for distance
            real*8 :: MPa1,DPa1,DPa2,QPa1,QPa2 !< prefactors for the multipole expansion
            real*8 :: pre1, pre2, pre3, pre4, pre5, pre6
            real*8, dimension(3) :: vort !< temp variables for vorticity (or better: alpha)

            ! tensors allow nice and short code and better comparison with (my!) theory
            real*8, dimension(3) :: m0, CP0 !< data structures for the monopole moments
            real*8, dimension(3,3) :: m1, CP1 !< data structures for the dipole moments
            real*8, dimension(3,3,3) :: m2, CP2 !< data structures for the quadrupole moments

            t=>tree_nodes(inode)

            dx = d(1)
            dy = d(2)
            dz = d(3)

            vort = [particle%data%alpha(1),particle%data%alpha(2),particle%data%alpha(3)]  ! need particle's vorticity for cross-product here

            m0 = [t%chargex,t%chargey,t%chargez]       ! monopole moment tensor
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

            pre1 = 0.0078125*(128.0*1.0*G_decomp(dist2,sig2,1.5D00) + 64.0*sig2*3.0*G_decomp(dist2,sig2,2.5D00) + 48.0*sig2**2*5.0*G_decomp(dist2,sig2,3.5D00) + &
                       40.0*sig2**3*7.0*G_decomp(dist2,sig2,4.5D00) - 70.0*sig2**4*9.0*G_decomp(dist2,sig2,5.5D00) + 315.0*sig2**5*11.0*G_decomp(dist2,sig2,6.5D00)) ! only Monopole u prefactor
                                                                                                                                                                     ! 2nd Dipole u prefactor
                                                                                                                                                                     ! 1st Monopole af prefactor

            pre2 = 0.0078125*(128.0*3.0*G_decomp(dist2,sig2,2.5D00) + 64.0*sig2*15.0*G_decomp(dist2,sig2,3.5D00) + 48.0*sig2**2*35.0*G_decomp(dist2,sig2,4.5D00) + &
                       40.0*sig2**3*63.0*G_decomp(dist2,sig2,5.5D00) - 70.0*sig2**4*99.0*G_decomp(dist2,sig2,6.5D00) + 315.0*sig2**5*143.0*G_decomp(dist2,sig2,7.5D00)) ! 1st Dipole u prefactor
                                                                                                                                                                        ! 2nd Monopole af prefactor (+dot_product)
                                                                                                                                                                        ! 2nd Dipole af prefactor

            pre3 = 0.5*0.0078125*(128.0*3.0*G_decomp(dist2,sig2,2.5D00) + 64.0*sig2*15.0*G_decomp(dist2,sig2,3.5D00) + 48.0*sig2**2*35.0*G_decomp(dist2,sig2,4.5D00) + &
                       40.0*sig2**3*63.0*G_decomp(dist2,sig2,5.5D00) - 70.0*sig2**4*99.0*G_decomp(dist2,sig2,6.5D00) + 315.0*sig2**5*143.0*G_decomp(dist2,sig2,7.5D00)) ! 1st Quadrupole u prefactor

            pre4 = 0.5*0.0078125*(128.0*15.0*G_decomp(dist2,sig2,3.5D00) + 64.0*sig2*105.0*G_decomp(dist2,sig2,4.5D00) + 48.0*sig2**2*315.0*G_decomp(dist2,sig2,5.5D00) + &
                       40.0*sig2**3*693.0*G_decomp(dist2,sig2,6.5D00) - 70.0*sig2**4*1287.0*G_decomp(dist2,sig2,7.5D00) + 315.0*sig2**5*2145.0*G_decomp(dist2,sig2,8.5D00)) ! 2nd Quadrupole u prefactor
                                                                                                                                                                           ! 2nd Quadrupole af prefactor (+...)

            pre5 = 0.0078125*(128.0*3.0*G_decomp(dist2,sig2,2.5D00) + 64.0*sig2*15.0*G_decomp(dist2,sig2,3.5D00) + 48.0*sig2**2*35.0*G_decomp(dist2,sig2,4.5D00) + &
                       40.0*sig2**3*63.0*G_decomp(dist2,sig2,5.5D00) - 70.0*sig2**4*99.0*G_decomp(dist2,sig2,6.5D00) + 315.0*sig2**5*143.0*G_decomp(dist2,sig2,7.5D00)) ! 1st Dipole af prefactor (+sumsum)

            pre6 = 0.5*0.0078125*(128.0*105.0*G_decomp(dist2,sig2,4.5D00) + 64.0*sig2*945.0*G_decomp(dist2,sig2,5.5D00) + 48.0*sig2**2*3465.0*G_decomp(dist2,sig2,6.5D00) + &
                       40.0*sig2**3*9009.0*G_decomp(dist2,sig2,7.5D00) - 70.0*sig2**4*19305.0*G_decomp(dist2,sig2,8.5D00) + 315.0*sig2**5*36465.0*G_decomp(dist2,sig2,9.5D00)) ! 3rd Quadrupole af prefactor (+...)


            MPa1 = pre2 * dot_product(d,CP0)   ! monopole prefactor for af
            DPa1 = pre5 * sum( (/ (sum( (/ (CP1(i2,i1)*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )  ! dipole prefators for af
            DPa2 = (CP1(1,1)+CP1(2,2)+CP1(3,3))
            QPa1 = pre6 * sum( (/ (sum( (/ (sum( (/ (CP2(i3,i2,i1)*d(i3),i3=1,3) /) )*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) ! quadrupole prefactors for af
            QPa2 = dot_product(d,CP2(1,1,:)+CP2(2,2,:)+CP2(3,3,:)+CP2(1,:,1)+CP2(2,:,2)+CP2(3,:,3)+CP2(:,1,1)+CP2(:,2,2)+CP2(:,3,3))

            u(1) = pre1 * (dy*m0(3)-dz*m0(2)) &                                                                                       ! MONOPOLE

                    + pre2 * sum( (/ ((m1(3,i1)*dy-m1(2,i1)*dz)*d(i1),i1=1,3) /) ) - pre1 * (m1(3,2)-m1(2,3)) &                 ! DIPOLE

                    - pre3 * ( sum( (/ (m2(3,i1,i1)*dy-m2(2,i1,i1)*dz,i1=1,3) /) ) + 2.0*sum( (/ ((m2(3,i1,2)-m2(2,i1,3))*d(i1),i1=1,3) /) ) ) + &
                      pre4 * sum( (/ (sum( (/ ((dy*m2(3,i2,i1)-dz*m2(2,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            u(2) = pre1 * (dz*m0(1)-dx*m0(3)) &                                                                                       ! MONOPOLE

                    + pre2 * sum( (/ ((m1(1,i1)*dz-m1(3,i1)*dx)*d(i1),i1=1,3) /) ) - pre1 * (m1(1,3)-m1(3,1)) &                 ! DIPOLE

                    - pre3 * ( sum( (/ (m2(1,i1,i1)*dz-m2(3,i1,i1)*dx,i1=1,3) /) ) + 2.0*sum( (/ ((m2(1,i1,3)-m2(3,i1,1))*d(i1),i1=1,3) /) ) ) + &
                      pre4 * sum( (/ (sum( (/ ((dz*m2(1,i2,i1)-dx*m2(3,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            u(3) = pre1 * (dx*m0(2)-dy*m0(1)) &                                                                                       ! MONOPOLE

                    + pre2 * sum( (/ ((m1(2,i1)*dx-m1(1,i1)*dy)*d(i1),i1=1,3) /) ) - pre1 * (m1(2,1)-m1(1,2)) &                 ! DIPOLE

                    - pre3 * ( sum( (/ (m2(2,i1,i1)*dx-m2(1,i1,i1)*dy,i1=1,3) /) ) + 2.0*sum( (/ ((m2(2,i1,1)-m2(1,i1,2))*d(i1),i1=1,3) /) ) ) + &
                      pre4 * sum( (/ (sum( (/ ((dx*m2(2,i2,i1)-dy*m2(1,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) )           ! QUADRUPOLE


            af(1) = Mpa1*dx - pre1 * CP0(1)  &                                                                                               ! MONOPOLE

                    + DPa1*dx - pre2 * ( dot_product(CP1(:,1)+CP1(1,:),d) + dx*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dx - &
                        pre4 * ( sum( (/ (sum( (/ ((CP2(i2,i1,1)+CP2(i2,1,i1)+CP2(1,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dx*QPa2) + &
                        pre3 * ( CP2(1,1,1)+CP2(2,2,1)+CP2(3,3,1)+CP2(1,1,1)+CP2(2,1,2)+CP2(3,1,3)+CP2(1,1,1)+CP2(1,2,2)+CP2(1,3,3) ) ! QUADRUPOLE


            af(2) = Mpa1*dy - pre1 * CP0(2)  &                                                                                               ! MONOPOLE

                    + DPa1*dy - pre2 * ( dot_product(CP1(:,2)+CP1(2,:),d) + dy*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dy - &
                        pre4 * ( sum( (/ (sum( (/ ((CP2(i2,i1,2)+CP2(i2,2,i1)+CP2(2,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dy*QPa2) + &
                        pre3 * ( CP2(1,1,2)+CP2(2,2,2)+CP2(3,3,2)+CP2(1,2,1)+CP2(2,2,2)+CP2(3,2,3)+CP2(2,1,1)+CP2(2,2,2)+CP2(2,3,3) ) ! QUADRUPOLE


            af(3) = Mpa1*dz - pre1 * CP0(3)  &                                                                                               ! MONOPOLE

                    + DPa1*dz - pre2 * ( dot_product(CP1(:,3)+CP1(3,:),d) + dz*DPa2 ) &                                              ! DIPOLE

                    + QPa1*dz - &
                        pre4 * ( sum( (/ (sum( (/ ((CP2(i2,i1,3)+CP2(i2,3,i1)+CP2(3,i2,i1))*d(i2),i2=1,3) /) )*d(i1),i1=1,3) /) ) + dz*QPa2) + &
                        pre3 * ( CP2(1,1,3)+CP2(2,2,3)+CP2(3,3,3)+CP2(1,3,1)+CP2(2,3,2)+CP2(3,3,3)+CP2(3,1,1)+CP2(3,2,2)+CP2(3,3,3) ) ! QUADRUPOLE

        end subroutine calc_6th_algebraic_condensed

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D 2nd order condensed algebraic kernel interaction
        !> of particle p with tree *particle*, results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_2nd_algebraic_condensed_direct(particle, inode, d, dist2, u, af)

            use module_pepc_types
            use treevars, only : tree_nodes
            implicit none

            integer, intent(in) :: inode
            type(t_particle), intent(inout) :: particle
            real*8, intent(in) :: d(3), dist2
            real*8, intent(out) :: u(1:3), af(1:3)

            type(t_tree_node_interaction_data), pointer :: t
            real*8, dimension(3) :: m0, CP0 !< data structures for the monopole moments
            real*8 :: dx, dy, dz, Gc25, MPa1
            real*8, dimension(3) :: vort !< temp variables for vorticity (or better: alpha)

            t=>tree_nodes(inode)

            dx = d(1)
            dy = d(2)
            dz = d(3)

            m0 = [t%chargex,t%chargey,t%chargez]       ! monopole moment tensor
            vort = [particle%data%alpha(1),particle%data%alpha(2),particle%data%alpha(3)]  ! need particle's vorticity for cross-product here
            CP0 = cross_prod(m0,vort)                  ! cross-product for 1st expansion term

            Gc25 = G_core(dist2,sig2,2.5D00)
            MPa1 = 3.0D00*G_core(dist2,sig2,3.5D00)*dot_product(d,CP0)

            u(1) = Gc25* (dy*m0(3)-dz*m0(2))
            u(2) = Gc25* (dz*m0(1)-dx*m0(3))
            u(3) = Gc25* (dx*m0(2)-dy*m0(1))

            af(1) = Mpa1*dx - Gc25*CP0(1)
            af(2) = Mpa1*dy - Gc25*CP0(2)
            af(3) = Mpa1*dz - Gc25*CP0(3)

        end subroutine calc_2nd_algebraic_condensed_direct

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D 2nd order condensed algebraic kernel interaction
        !> of particle p with tree node inode, results are returned in u and af
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_6th_algebraic_condensed_direct(particle, inode, d, dist2, u, af)
            use module_pepc_types
            use treevars
            use module_interaction_specific_types
            implicit none

            type(t_particle), intent(in) :: particle
            integer, intent(in) :: inode !< index of particle to interact with
            real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
            real*8, intent(out) ::  u(1:3), af(1:3)

            type(t_tree_node_interaction_data), pointer :: t

            real*8 :: dx, dy, dz !< temp variables for distance
            real*8 :: sig4, sig6, sig8, sig10, Gd15, Gd25, Gd35, Gd45, Gd55, Gd65, Gd75
            real*8 :: pre1, pre2, MPa1
            real*8, dimension(3) :: vort !< temp variables for vorticity (or better: alpha)
            real*8, dimension(3) :: m0, CP0 !< data structures for the monopole moments


            t=>tree_nodes(inode)

            dx = d(1)
            dy = d(2)
            dz = d(3)

            vort = [particle%data%alpha(1),particle%data%alpha(2),particle%data%alpha(3)]  ! need particle's vorticity for cross-product here

            m0 = [t%chargex,t%chargey,t%chargez]       ! monopole moment tensor
            CP0 = cross_prod(m0,vort)                  ! cross-product for 1st expansion term


            ! precompute kernel function evaluations of various order
            Gd15 = G_decomp(dist2,sig2,1.5D00)
            Gd25 = G_decomp(dist2,sig2,2.5D00)
            Gd35 = G_decomp(dist2,sig2,3.5D00)
            Gd45 = G_decomp(dist2,sig2,4.5D00)
            Gd55 = G_decomp(dist2,sig2,5.5D00)
            Gd65 = G_decomp(dist2,sig2,6.5D00)
            Gd75 = G_decomp(dist2,sig2,7.5D00)

            sig4  = sig2*sig2
            sig6  = sig4*sig2
            sig8  = sig6*sig2
            sig10 = sig8*sig2

            pre1 = 0.0078125*(128.0*Gd15 + 192.0*sig2*Gd25 + 240.0*sig4*Gd35 + 280.0*sig6*Gd45 - 630.0*sig8*Gd55 + 3465.0*sig10*Gd65)
            pre2 = 0.0078125*(384.0*Gd25 + 960.0**sig2*Gd35 + 1680.0*sig4*Gd45 + 2520.0*sig6*Gd55 - 6930.0*sig8*Gd65 + 45045.0*sig10*Gd75)

            MPa1 = pre2 * dot_product(d,CP0)   ! monopole prefactor for af

            u(1) = pre1 * (dy*m0(3)-dz*m0(2))                                                                                        ! MONOPOLE
            u(2) = pre1 * (dz*m0(1)-dx*m0(3))                                                                                        ! MONOPOLE
            u(3) = pre1 * (dx*m0(2)-dy*m0(1))                                                                                        ! MONOPOLE

            af(1) = Mpa1*dx - pre1 * CP0(1)                                                                                                 ! MONOPOLE
            af(2) = Mpa1*dy - pre1 * CP0(2)                                                                                                 ! MONOPOLE
            af(3) = Mpa1*dz - pre1 * CP0(3)                                                                                                 ! MONOPOLE

        end subroutine calc_6th_algebraic_condensed_direct




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

  end module module_interaction_specific
