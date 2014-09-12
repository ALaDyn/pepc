! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2012 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains types for building the sim domain form walls
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_geometry
    use module_cmdline
    use variables
    use helper
    use module_mirror_boxes, only: mirror_box_layers,periodicity


    implicit none

    contains

!======================================================================================

    SUBROUTINE init_boundaries()
        implicit none

        real(KIND=8),allocatable :: x0(:,:)
        real(KIND=8),allocatable :: e1(:,:)
        real(KIND=8),allocatable :: e2(:,:)
        real(KIND=8),allocatable :: n(:,:)

        integer,allocatable :: type(:)

        integer,allocatable :: opposite_bnd(:)
        logical,allocatable :: reflux_particles(:)
        integer,allocatable :: nwp(:)

        integer :: nbnd,nbnd_max
        integer :: rc,ib,fid=12

        namelist /geometry/ x0,e1,e2,n,type,opposite_bnd,reflux_particles,nwp,nbnd

        nbnd_max=1000
        allocate(x0(nbnd_max,3),stat=rc)
        allocate(e1(nbnd_max,3),stat=rc)
        allocate(e2(nbnd_max,3),stat=rc)
        allocate(n(nbnd_max,3),stat=rc)
        allocate(type(nbnd_max),stat=rc)
        allocate(opposite_bnd(nbnd_max),stat=rc)
        allocate(reflux_particles(nbnd_max),stat=rc)
        allocate(nwp(nbnd_max),stat=rc)

        IF(root) write(*,'(a,a)') " == reading parameter file, section geometry: ", trim(input_file)
        open(fid,file=trim(input_file))
        read(fid,NML=geometry)
        rewind(fid)

        deallocate(x0)
        deallocate(e1)
        deallocate(e2)
        deallocate(n)
        deallocate(type)
        deallocate(opposite_bnd)
        deallocate(reflux_particles)
        deallocate(nwp)

        allocate(x0(nbnd,3),stat=rc)
        allocate(e1(nbnd,3),stat=rc)
        allocate(e2(nbnd,3),stat=rc)
        allocate(n(nbnd,3),stat=rc)
        allocate(type(nbnd),stat=rc)
        allocate(opposite_bnd(nbnd),stat=rc)
        allocate(reflux_particles(nbnd),stat=rc)
        allocate(nwp(nbnd),stat=rc)

        nwp=0
        reflux_particles=.false.
        opposite_bnd=0
        type=0
        x0=0.
        e1=0.
        e2=0.
        n=0.

        read(fid,NML=geometry)
        close(fid)

        IF(nbnd<=0) THEN
            IF (root) write(*,'(a)') " number of boundaries not set or set to invalid value "
            STOP
        ELSE
            IF (root) write(*,'(a,i3,a)') " == initializing ",nbnd," boundaries"
        END IF

        nb=nbnd
        allocate(boundaries(nb),stat=rc)

        DO ib=1,nb
            CALL init_boundary(x0(ib,:),e1(ib,:),e2(ib,:),n(ib,:),type(ib),ib,boundaries(ib))
            IF (nwp(ib)>0) CALL init_wall(nwp(ib),boundaries(ib))
            CALL set_refluxing(reflux_particles(ib),boundaries(ib))
        END DO

        DO ib=1,nb
            IF (opposite_bnd(ib)/=0) THEN
                CALL set_periodic_bc(boundaries(ib),boundaries(opposite_bnd(ib)))
            END IF
        END DO

        CALL check_boundaries()

        deallocate(x0)
        deallocate(e1)
        deallocate(e2)
        deallocate(n)
        deallocate(type)
        deallocate(opposite_bnd)
        deallocate(reflux_particles)
        deallocate(nwp)

    END SUBROUTINE init_boundaries

!======================================================================================

    SUBROUTINE check_boundaries()
        implicit none

        integer :: ib
        real*8 :: test
        real*8 :: eps=1.0e-12

        IF (root) THEN
            DO ib=1,nb
                CALL check_boundary(boundaries(ib))

                IF (boundaries(ib)%type==2) THEN
                    IF(boundaries(ib)%opp_bnd==0) THEN
                        write(*,'(a,i3,a)')"No opposing boundary set for boundary ",boundaries(ib)%indx,"."
                        STOP
                    END IF
                    test=dotproduct(boundaries(ib)%n,boundaries(boundaries(ib)%opp_bnd)%n)
                    IF (real_unequal(test,-1._8,eps)) THEN
                        write(*,*)
                        write(*,'(a,i3,a,i3,a)')"Boundaries ",boundaries(ib)%indx," and ",boundaries(ib)%opp_bnd," are set as a pair of opposing boundaries, but are not parallel."
                        write(*,*)boundaries(ib)%n
                        write(*,*)boundaries(boundaries(ib)%opp_bnd)%n
                        write(*,*)
                        STOP
                    END IF
                END IF
            END DO
        END IF

    END SUBROUTINE check_boundaries

!======================================================================================

    integer function count_wallparticles()
        implicit none

        type(t_boundary)::wall
        integer :: ib,tnwp_aux

        tnwp_aux=0
        DO ib=1,nb
            wall=boundaries(ib)
            tnwp_aux=tnwp_aux+wall%nwp
        END DO

        count_wallparticles=tnwp_aux

   END function count_wallparticles

 !======================================================================================

     SUBROUTINE hit_wall(p,wall,hit,x_hit,x_hit_rel)
        use module_pepc_types
        use variables
        implicit none

        type(t_particle), intent(in) :: p
        type(t_boundary), intent(in) :: wall
        logical,intent(out) :: hit
        real*8 :: xold(3),x(3),n(3),deltax(3)
        real(KIND=8), intent(out) :: x_hit(3),x_hit_rel(2)

        hit=.false.
        x_hit=0.0_8
        x_hit_rel=0.0_8

        IF (species(p%data%species)%physical_particle .eqv. .false.) THEN
            RETURN !no hit
        END IF

        n=wall%n
        deltax=p%x-wall%x0

        IF (dotproduct(n,deltax)>0) THEN
            RETURN !no hit
        END IF

        x = p%x
        xold = p%x - dt*p%data%v

        CALL get_intersect(xold,x,wall,x_hit)
        CALL check_hit(x_hit(1),x_hit(2),x_hit(3),wall,hit,x_hit_rel)


        RETURN

    END SUBROUTINE hit_wall

!======================================================================================

    SUBROUTINE init_wall(nwp,boundary)
        implicit none

        type(t_boundary), intent(inout) :: boundary
        integer, intent(in) :: nwp

        IF (boundary%type/=0) THEN
            write(*,*) "Problem with boundary",boundary%indx
            write(*,*) "Wallparticles can only be initialized for absorbing walls (type=0)."
            STOP
        END IF

        boundary%nwp=nwp

    END SUBROUTINE init_wall

 !======================================================================================

    SUBROUTINE init_boundary(x0,e1,e2,n,typ,indx,wall)
        implicit none

        type(t_boundary), intent(inout) :: wall
        real*8, intent(in), dimension(3) :: x0,e1,e2,n
        integer, intent(in) :: typ,indx

        wall%x0=x0
        wall%e1=e1
        wall%e2=e2
        wall%n=n/sqrt(dotproduct(n,n))
        wall%type=typ
        wall%indx=indx
        wall%reflux_particles=.false.

    END SUBROUTINE init_boundary

 !======================================================================================

    SUBROUTINE set_refluxing(reflux_particles,boundary)
        implicit none

        type(t_boundary), intent(inout) :: boundary
        logical, intent(in) :: reflux_particles

        IF (reflux_particles .eqv. .false.) THEN
            boundary%reflux_particles=reflux_particles
        ELSE
            IF ((boundary%type==0) .OR. (boundary%type==3) .OR. (boundary%type==4)) THEN
                boundary%reflux_particles=reflux_particles
            ELSE
                write(*,*) "Problem with boundary",boundary%indx
                write(*,*) "Refluxing conditions can only be set for absorbing boundaries (type=0, type=4 or type=3)."
                STOP
            END IF
        END IF
    END SUBROUTINE set_refluxing

!======================================================================================

    SUBROUTINE set_periodic_bc(wall,opp_wall)
        implicit none

        type(t_boundary), intent(inout) :: wall,opp_wall


        IF ((wall%type/=2) .OR. (opp_wall%type/=2)) THEN
            write(*,'(a,i3,a,i3,a)') "Error while trying to set periodic boundaries for boundaries ",wall%indx," and ",opp_wall%indx,": Both walls have to have type 2"
            STOP
        END IF

        wall%opp_bnd=opp_wall%indx
        opp_wall%opp_bnd=wall%indx
        wall%dist=dotproduct(wall%n,(opp_wall%x0-wall%x0))
        opp_wall%dist=wall%dist

    END SUBROUTINE set_periodic_bc

!======================================================================================

    SUBROUTINE check_boundary(wall)
        type(t_boundary), intent(inout) :: wall
        real*8 :: eps=1.0e-12
        real*8 :: test(2)
        logical :: ok

        test(1)=dotproduct(wall%e1,wall%n)
        test(2)=dotproduct(wall%e2,wall%n)

        ok = (real_equal_zero(test(1),eps)) .AND. (real_equal_zero(test(2),eps))

        IF (ok .eqv. .false.) THEN
            write(*,*) "Boundary set incorrectly:"
            write(*,*) "x0:", wall%x0
            write(*,*) "e1:", wall%e1
            write(*,*) "e2:", wall%e2
            write(*,*) "n:", wall%n
            write(*,*) "The surface normal is not perpENDicular to the plane."
            STOP
        END IF

    END SUBROUTINE check_boundary

!======================================================================================

    SUBROUTINE check_hit(px,py,pz,wall,hit,x_hit_rel)
        type(t_boundary), intent(in) :: wall
        logical,intent(out) :: hit
        real*8,intent(in) :: px,py,pz
        real*8 :: eps=1.0e-12
        real*8 :: y(3),a(3),b(3),lambda,mu,dist
        real(KIND=8), intent(out) :: x_hit_rel(2)

        x_hit_rel=0.0_8
        hit = .false.

        a=wall%e1
        b=wall%e2

        y(1)=px-wall%x0(1)
        y(2)=py-wall%x0(2)
        y(3)=pz-wall%x0(3)

        lambda=0.
        mu=0.

        dist=dotproduct(wall%n,y)

        IF (real_unequal_zero(dist,eps)) THEN !Punkt liegt nicht in der Ebene
            RETURN
        ELSE !Punkt liegt in der Ebene. Check ob er auch in der Flaeche liegt:
            IF (real_unequal_zero(a(1),eps) .AND. (real_unequal_zero(b(2),eps) .OR. real_unequal_zero(b(3),eps)))  THEN
                IF (real_unequal_zero(b(2),eps)) THEN
                    mu = (y(2)/b(2) - (a(2)*y(1))/(b(2)*a(1))) / (1 - (b(1)*a(2) / (b(2)*a(1))))
                    lambda = (y(1) - (mu*b(1))) / a(1)
                ELSE IF (real_unequal_zero(b(3),eps)) THEN
                    mu = (y(3)/b(3) - (a(3)*y(1))/(b(3)*a(1))) / (1 - (b(1)*a(3) / (b(3)*a(1))))
                    lambda = (y(1) - (mu*b(1))) / a(1)
                END IF
            ELSE IF (real_unequal_zero(a(2),eps) .AND. (real_unequal_zero(b(1),eps) .OR. real_unequal_zero(b(3),eps)))  THEN
                IF (real_unequal_zero(b(1),eps)) THEN
                    mu = (y(1)/b(1) - (a(1)*y(2))/(b(1)*a(2))) / (1 - (b(2)*a(1) / (b(1)*a(2))))
                    lambda = (y(2) - (mu*b(2))) / a(2)
                ELSE IF (real_unequal_zero(b(3),eps)) THEN
                    mu = (y(3)/b(3) - (a(3)*y(2))/(b(3)*a(2))) / (1 - (b(2)*a(3) / (b(3)*a(2))))
                    lambda = (y(2) - (mu*b(2))) / a(2)
                END IF
            ELSE IF (real_unequal_zero(a(3),eps) .AND. (real_unequal_zero(b(2),eps) .OR. real_unequal_zero(b(1),eps)))  THEN
                IF (real_unequal_zero(b(2),eps)) THEN
                    mu = (y(2)/b(2) - (a(2)*y(3))/(b(2)*a(3))) / (1 - (b(3)*a(2) / (b(2)*a(3))))
                    lambda = (y(3) - (mu*b(3))) / a(3)
                ELSE IF (real_unequal_zero(b(1),eps)) THEN
                    mu = (y(1)/b(1) - (a(1)*y(3))/(b(1)*a(3))) / (1 - (b(3)*a(1) / (b(1)*a(3))))
                    lambda = (y(3) - (mu*b(3))) / a(3)
                END IF
            END IF

            IF ((lambda<1.).AND.(mu<1.) .AND. (mu>0.) .AND. (lambda>0.)) THEN ! Der Punkt liegt in der Flaeche
                hit=.true.
                x_hit_rel(1) = lambda
                x_hit_rel(2) = mu
            END IF
        END IF
    END SUBROUTINE check_hit

!======================================================================================


    SUBROUTINE get_intersect(xold,xnew,wall,intersect)
        implicit none

        type(t_boundary), intent(in) :: wall
        real*8, intent(in), dimension(3) :: xold,xnew
        real*8 :: n(3),x0(3),lambda,a,b,dx(3)
        real*8, intent(out), dimension(3) :: intersect

        x0= wall%x0
        n = wall%n
        dx = (xold-x0)
        a = dotproduct(dx,n)
        dx = (xnew-x0)
        b = dotproduct(dx,n)
        lambda = abs(a) / (abs(a)+abs(b))
        intersect = xold + (xnew-xold)*lambda

    END SUBROUTINE get_intersect

END module module_geometry
