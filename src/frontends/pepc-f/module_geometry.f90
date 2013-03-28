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

    implicit none

    contains

!======================================================================================

    subroutine init_boundaries()
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
            call init_boundary(x0(ib,:),e1(ib,:),e2(ib,:),n(ib,:),type(ib),ib,boundaries(ib))
            IF (nwp(ib)>0) call init_wall(nwp(ib),boundaries(ib))
            call set_refluxing(reflux_particles(ib),boundaries(ib))
        END DO

        DO ib=1,nb
            IF (opposite_bnd(ib)/=0) THEN
                call set_periodic_bc(boundaries(ib),boundaries(opposite_bnd(ib)))
            END IF
        END DO

        call check_boundaries()
        tnpps(0)=count_wallparticles()
        call add_wallparticles_to_boundaries()

        deallocate(x0)
        deallocate(e1)
        deallocate(e2)
        deallocate(n)
        deallocate(type)
        deallocate(opposite_bnd)
        deallocate(reflux_particles)
        deallocate(nwp)

    end subroutine init_boundaries

!======================================================================================

    subroutine check_boundaries()
        implicit none

        integer :: ib
        real*8 :: test
        real*8 :: eps=1.0e-19

        if (root) then
            do ib=1,nb
                call check_boundary(boundaries(ib))

                if (boundaries(ib)%type==2) then
                    if(boundaries(ib)%opp_bnd==0) then
                        write(*,'(a,i3,a)')"No opposing boundary set for boundary ",boundaries(ib)%indx,"."
                        STOP
                    end if
                    test=dotproduct(boundaries(ib)%n,boundaries(boundaries(ib)%opp_bnd)%n)
                    if ((test<-1.-eps).or.(test>-1.+eps)) then
                        write(*,*)
                        write(*,'(a,i3,a,i3,a)')"Boundaries ",boundaries(ib)%indx," and ",boundaries(ib)%opp_bnd," are set as a pair of opposing boundaries, but are not parallel."
                        write(*,*)boundaries(ib)%n
                        write(*,*)boundaries(boundaries(ib)%opp_bnd)%n
                        write(*,*)
                        STOP
                    end if
                end if
            end do
        end if

    end subroutine check_boundaries

!======================================================================================

    integer function count_wallparticles()
        implicit none

        type(t_boundary)::wall
        integer :: ib,tnwp_aux

        tnwp_aux=0
        do ib=1,nb
            wall=boundaries(ib)
            tnwp_aux=tnwp_aux+wall%nwp
        end do

        count_wallparticles=tnwp_aux

   end function count_wallparticles

 !======================================================================================

    subroutine add_wallparticles_to_boundaries()
        implicit none

        integer :: ib,label,i

        label=-1
        do ib=1,nb
            if (boundaries(ib)%nwp/=0) then
                do i=1,boundaries(ib)%nwp
                    boundaries(ib)%wp_labels(i)=label
                    label=label-1
                end do
            end if
        end do

    end subroutine add_wallparticles_to_boundaries

!======================================================================================

    subroutine hit_wall(p,wall,hit)
        use module_pepc_types
        use variables
        implicit none

        type(t_particle), intent(in) :: p
        type(t_boundary), intent(in) :: wall
        logical,intent(out) :: hit
        real*8 :: xold(3),x(3),intersect(3),n(3),deltax(3)

        hit=.false.

        if (p%data%species == 0) then
             hit = .false.
            return
        end if

        n=wall%n
        deltax=p%x-wall%x0
        if (dotproduct(n,deltax)>0) then
            hit=.false.
            return
        end if

        x=p%x
        xold=p%x - dt*p%data%v
        call get_intersect(xold,x,wall,intersect)
        call check_hit(intersect(1),intersect(2),intersect(3),wall,hit)

        return

    end subroutine hit_wall

!======================================================================================

    subroutine init_wall(nwp,boundary)
        implicit none

        type(t_boundary), intent(inout) :: boundary
        integer, intent(in) :: nwp
        integer :: rc

        IF (boundary%type/=0) THEN
            write(*,*) "Problem with boundary",boundary%indx
            write(*,*) "Wallparticles can only be initialized for absorbing walls (type=0)."
            STOP
        END IF

        boundary%nwp=nwp
        deallocate(boundary%wp_labels)
        allocate(boundary%wp_labels(nwp),stat=rc)
        boundary%wp_labels=0

    end subroutine init_wall

 !======================================================================================

    subroutine init_boundary(x0,e1,e2,n,typ,indx,wall)
        implicit none

        type(t_boundary), intent(inout) :: wall
        real*8, intent(in), dimension(3) :: x0,e1,e2,n
        integer, intent(in) :: typ,indx
        integer :: rc

        wall%x0=x0
        wall%e1=e1
        wall%e2=e2
        wall%n=n/dotproduct(n,n)
        wall%type=typ
        wall%indx=indx
        wall%reflux_particles=.false.
        allocate(wall%wp_labels(0),stat=rc)

    end subroutine init_boundary

 !======================================================================================

    subroutine set_refluxing(reflux_particles,boundary)
        implicit none

        type(t_boundary), intent(inout) :: boundary
        logical, intent(in) :: reflux_particles

        IF (reflux_particles .eqv. .false.) THEN
            boundary%reflux_particles=reflux_particles
        ELSE
            IF ((boundary%type==0) .or. (boundary%type==3) .or. (boundary%type==4)) THEN
                boundary%reflux_particles=reflux_particles
            ELSE
                write(*,*) "Problem with boundary",boundary%indx
                write(*,*) "Refluxing conditions can only be set for absorbing boundaries (type=0, type=4 or type=3)."
                STOP
            END IF
        END IF
    end subroutine set_refluxing

!======================================================================================

    subroutine set_periodic_bc(wall,opp_wall)
        implicit none

        type(t_boundary), intent(inout) :: wall,opp_wall


        if ((wall%type/=2) .or. (opp_wall%type/=2)) then
            write(*,'(a,i3,a,i3,a)') "Error while trying to set periodic boundaries for boundaries ",wall%indx," and ",opp_wall%indx,": Both walls have to have type 2"
            STOP
        end if

        wall%opp_bnd=opp_wall%indx
        opp_wall%opp_bnd=wall%indx
        wall%dist=dotproduct(wall%n,(opp_wall%x0-wall%x0))
        opp_wall%dist=wall%dist

    end subroutine set_periodic_bc

!======================================================================================

    subroutine check_boundary(wall)
        type(t_boundary), intent(inout) :: wall

        real*8 :: test(2)
        logical :: ok

        test(1)=dotproduct(wall%e1,wall%n)
        test(2)=dotproduct(wall%e2,wall%n)

        ok= (test(1)==0) .and. (test(2)==0)

        if (ok .eqv. .false.) then
            write(*,*) "Boundary set incorrectly:"
            write(*,*) "x0:", wall%x0
            write(*,*) "e1:", wall%e1
            write(*,*) "e2:", wall%e2
            write(*,*) "n:", wall%n
            write(*,*) "The surface normal is not perpendicular to the plane."
            STOP
        end if

    end subroutine check_boundary

!======================================================================================

    subroutine check_hit(px,py,pz,wall,hit)
        type(t_boundary), intent(in) :: wall
        logical,intent(out) :: hit
        real*8,intent(in) :: px,py,pz
        real*8 :: eps=1.0e-19
        real*8 :: y(3),a(3),b(3),lambda,mu,dist

        a=wall%e1
        b=wall%e2

        y(1)=px-wall%x0(1)
        y(2)=py-wall%x0(2)
        y(3)=pz-wall%x0(3)

        dist=dotproduct(wall%n,y)

        if ((a(1)/=0) .and. ((b(2)/=0) .or. (b(3)/=0)))  then
            if (b(2)/=0) then
                mu = (y(2)/b(2) - (a(2)*y(1))/(b(2)*a(1))) / (1 - (b(1)*a(2) / (b(2)*a(1))))
                lambda = (y(1) - (mu*b(1))) / a(1)
            else if (b(3)/=0) then
                mu = (y(3)/b(3) - (a(3)*y(1))/(b(3)*a(1))) / (1 - (b(1)*a(3) / (b(3)*a(1))))
                lambda = (y(1) - (mu*b(1))) / a(1)
            end if
        else if ((a(2)/=0) .and. ((b(1)/=0) .or. (b(3)/=0)))  then
            if (b(1)/=0) then
                mu = (y(1)/b(1) - (a(1)*y(2))/(b(1)*a(2))) / (1 - (b(2)*a(1) / (b(1)*a(2))))
                lambda = (y(2) - (mu*b(2))) / a(2)
            else if (b(3)/=0) then
                mu = (y(3)/b(3) - (a(3)*y(2))/(b(3)*a(2))) / (1 - (b(2)*a(3) / (b(3)*a(2))))
                lambda = (y(2) - (mu*b(2))) / a(2)
            end if
        else if ((a(3)/=0) .and. ((b(2)/=0) .or. (b(1)/=0)))  then
            if (b(2)/=0) then
                 mu = (y(2)/b(2) - (a(2)*y(3))/(b(2)*a(3))) / (1 - (b(3)*a(2) / (b(2)*a(3))))
                 lambda = (y(3) - (mu*b(3))) / a(3)
            else if (b(1)/=0) then
                 mu = (y(1)/b(1) - (a(1)*y(3))/(b(1)*a(3))) / (1 - (b(3)*a(1) / (b(1)*a(3))))
                lambda = (y(3) - (mu*b(3))) / a(3)
            end if
        end if

        hit=.false.
        if ((dist<eps).and.(dist>-eps)) then
            !write(*,*) "Punkt liegt in der Ebene"
            if ((lambda<1.).and.(mu<1.) .and. (mu>0.) .and. (lambda>0.)) then
                !write(*,*) "Der Punkt liegt in der Flaeche"
                hit=.true.
            else
                !write(*,*) "Der Punkt liegt nicht in der Flaeche",lambda,mu,wall%indx
                !write(*,*) px,py,pz
            end if
        else
            !write(*,*) "Punkt liegt nicht in der Ebene",dist
            !write(*,*) px,py,pz
        end if
    end subroutine check_hit

!======================================================================================

    subroutine get_intersect(xold,xnew,wall,intersect)
        implicit none

        type(t_boundary), intent(in) :: wall
        real*8, intent(in), dimension(3) :: xold,xnew
        real*8 :: lambda,n(3),x0(3)
        real*8, intent(out), dimension(3) :: intersect

        x0=wall%x0
        n=wall%n
        lambda=(dotproduct(n,x0)-dotproduct(n,xold)) / dotproduct(n, xnew-xold)
        intersect=xold + lambda*(xnew-xold)

    end subroutine get_intersect

end module module_geometry
