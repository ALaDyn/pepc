! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2023 Juelich Supercomputing Centre, 
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

module manipulate_particles

    use module_pepc_kinds
    implicit none

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Initialize particles with different setup (choose via ispecial)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine special_start()

        use physvars
        use files
        use mpi
        implicit none

        integer :: ierr
        integer(kind_particle) :: j, k, ind, ind0, i, m, l
        real :: par_rand_res
        real(kind_physics) :: part_2d, rc, xi1, xi2, xi, part_3d, eta1, eta2, eta, v(3), xt, yt, zt
        real(kind_physics), dimension(3,3) :: D1, D2, D3, D4   !< rotation matrices for ring setup
        real(kind_physics), allocatable :: xp(:), yp(:), zp(:), volp(:), wxp(:), wyp(:), wzp(:)  !< helper arrays for ring setups

        ! weird helper variables for sphere setup (ask Holger Dachsel)
        real(kind_physics) ::  a,b,c,cth,sth,cphi,sphi,s
        real(kind_physics), parameter :: zero=0.d0, one=1.d0, mone=-one, two=2.d0, five=5.d0, nine=9.d0

        !     interface
        !         subroutine par_rand(res, iseed)
        !            real, intent(inout) :: res
        !            integer, intent(in), optional :: iseed
        !         end subroutine
        !     end interface


        ! Set up particle data
        config: select case(ispecial)
        case(1)                               ! Vortex ring setup, side-by-side

            allocate(xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

            j = 0

            do k = 1, nc
                part_2d = 2*pi/(8*k)
                rc = (1+12*k**2)/(6*k)*rl

                do l = 1,8*k
                    j = j+1
                    xi1 = part_2d*(l-1)
                    xi2 = part_2d*l
                    xi = (xi2-xi1)/2+xi1
                    xp(j) = rc*cos(xi)
                    yp(j) = rc*sin(xi)
                    zp(j) = 0
                    volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
                    wxp(j) = 0.
                    wyp(j) = 0.
                    wzp(j) = g*exp(-(rc/rmax)**2)
                end do
            end do

            xp(ns) = 0.
            yp(ns) = 0.
            zp(ns) = 0.
            wxp(ns) = 0.
            wyp(ns) = 0.
            wzp(ns) = g
            volp(ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

            j = 0
            ind0 = 0
            ind = 0
            part_3d = 2*pi/Nphi
            do m = 1, Nphi
                eta1 = part_3d*(m-1)
                eta2 = part_3d*m
                eta = (eta2 - eta1)/2+eta1
                v(1) = cos(eta2)
                v(2) = sin(eta2)
                v(3) = 0
                D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                v(1) = cos(eta)
                v(2) = sin(eta)
                v(3) = 0
                D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
                D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                D4 = matmul(D1,D3)
                do i = 1, Ns
                    if (m==1) then
                        v(1) = xp(i) + (r_torus+rmax)*cos(eta)
                        v(2) = yp(i) + (r_torus+rmax)*sin(eta)
                        v(3) = zp(i);
                        xp(i) = dot_product(v,D2(1:3,1))
                        yp(i) = dot_product(v,D2(1:3,2))
                        zp(i) = dot_product(v,D2(1:3,3))
                        v(1) = wxp(i);
                        v(2) = wyp(i);
                        v(3) = wzp(i);
                        wxp(i) = dot_product(v,D2(1:3,1))
                        wyp(i) = dot_product(v,D2(1:3,2))
                        wzp(i) = dot_product(v,D2(1:3,3))
                    else
                        v(1) = xp(i)
                        v(2) = yp(i)
                        v(3) = zp(i)
                        xp(i) = dot_product(v,D4(1:3,1))
                        yp(i) = dot_product(v,D4(1:3,2))
                        zp(i) = dot_product(v,D4(1:3,3))
                        v(1) = wxp(i)
                        v(2) = wyp(i)
                        v(3) = wzp(i)
                        wxp(i) = dot_product(v,D4(1:3,1))
                        wyp(i) = dot_product(v,D4(1:3,2))
                        wzp(i) = dot_product(v,D4(1:3,3))
                    end if
                    ind0 = ind0 + 1
                    if (mod(ind0-1,n_cpu) == my_rank) then
                        ind = ind + 1
                        if (ind .gt.np) then
                            write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                        end if
                        vortex_particles(ind)%x(1) = xp(i) - (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(2) = yp(i)
                        vortex_particles(ind)%x(3) = zp(i)
                        vortex_particles(ind)%data%alpha(1) = wxp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(2) = wyp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(3) = wzp(i)*volp(i)
                        ind = ind + 1
                        vortex_particles(ind)%x(1) = xp(i) + (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(2) = yp(i)
                        vortex_particles(ind)%x(3) = zp(i)
                        vortex_particles(ind)%data%alpha(1) = wxp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(2) = wyp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(3) = wzp(i)*volp(i)
                    end if
                end do
            end do
            np = ind
            deallocate(xp,yp,zp,volp,wxp,wyp,wzp)

        case(2) ! Vortex ring setup, offset collision

            allocate(xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

            j = 0

            do k = 1, nc
                part_2d = 2*pi/(8*k)
                rc = (1+12*k**2)/(6*k)*rl

                do l = 1,8*k
                    j = j+1
                    xi1 = part_2d*(l-1)
                    xi2 = part_2d*l
                    xi = (xi2-xi1)/2+xi1
                    xp(j) = rc*cos(xi)
                    yp(j) = rc*sin(xi)
                    zp(j) = 0
                    volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
                    wxp(j) = 0.
                    wyp(j) = 0.
                    wzp(j) = g*exp(-(rc/rmax)**2)
                end do
            end do

            xp(ns) = 0.
            yp(ns) = 0.
            zp(ns) = 0.
            wxp(ns) = 0.
            wyp(ns) = 0.
            wzp(ns) = g
            volp(ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

            j = 0
            ind0 = 0
            ind = 0
            part_3d = 2*pi/Nphi
            do m = 1, Nphi
                eta1 = part_3d*(m-1)
                eta2 = part_3d*m
                eta = (eta2 - eta1)/2+eta1
                v(1) = cos(eta2)
                v(2) = sin(eta2)
                v(3) = 0
                D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                v(1) = cos(eta)
                v(2) = sin(eta)
                v(3) = 0
                D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
                D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                D4 = matmul(D1,D3)
                do i = 1, Ns
                    if (m==1) then
                        v(1) = xp(i) + (r_torus+rmax)*cos(eta)
                        v(2) = yp(i) + (r_torus+rmax)*sin(eta)
                        v(3) = zp(i);
                        xp(i) = dot_product(v,D2(1:3,1))
                        yp(i) = dot_product(v,D2(1:3,2))
                        zp(i) = dot_product(v,D2(1:3,3))
                        v(1) = wxp(i);
                        v(2) = wyp(i);
                        v(3) = wzp(i);
                        wxp(i) = dot_product(v,D2(1:3,1))
                        wyp(i) = dot_product(v,D2(1:3,2))
                        wzp(i) = dot_product(v,D2(1:3,3))
                    else
                        v(1) = xp(i)
                        v(2) = yp(i)
                        v(3) = zp(i)
                        xp(i) = dot_product(v,D4(1:3,1))
                        yp(i) = dot_product(v,D4(1:3,2))
                        zp(i) = dot_product(v,D4(1:3,3))
                        v(1) = wxp(i)
                        v(2) = wyp(i)
                        v(3) = wzp(i)
                        wxp(i) = dot_product(v,D4(1:3,1))
                        wyp(i) = dot_product(v,D4(1:3,2))
                        wzp(i) = dot_product(v,D4(1:3,3))
                    end if
                    ind0 = ind0 + 1
                    if (mod(ind0-1,n_cpu) == my_rank) then
                        ind = ind + 1
                        if (ind .gt. np) then
                            write(*,*) 'something is wrong here: to many particles in init of first ring',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                        end if
                        vortex_particles(ind)%x(1) = xp(i) - (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(2) = yp(i) - (torus_offset(2)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(3) = zp(i) - (torus_offset(3)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%data%alpha(1) = wxp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(2) = wyp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(3) = wzp(i)*volp(i)
                    end if
                end do
            end do

            j = 0

            do k = 1, nc
                part_2d = 2*pi/(8*k)
                rc = (1+12*k**2)/(6*k)*rl

                do l = 1,8*k
                    j = j+1
                    xi1 = part_2d*(l-1)
                    xi2 = part_2d*l
                    xi = (xi2-xi1)/2+xi1
                    xp(j) = rc*cos(xi)
                    yp(j) = rc*sin(xi)
                    zp(j) = 0
                    volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
                    wxp(j) = 0.
                    wyp(j) = 0.
                    wzp(j) = G*exp(-(rc/rmax)**2)
                end do
            end do

            xp(ns) = 0.
            yp(ns) = 0.
            zp(ns) = 0.
            wxp(ns) = 0.
            wyp(ns) = 0.
            wzp(ns) = g
            volp(Ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

            j = 0
            ind0 = 0
            part_3d = 2*pi/Nphi
            do m = 1, Nphi
                eta1 = part_3d*(m-1)
                eta2 = part_3d*m
                eta = (eta2 - eta1)/2+eta1
                v(1) = cos(eta2)
                v(2) = sin(eta2)
                v(3) = 0
                D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                v(1) = cos(eta)
                v(2) = sin(eta)
                v(3) = 0
                D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
                D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
                D4 = matmul(D1,D3)
                do i = 1, Ns
                    if (m==1) then
                        v(1) = xp(i) + (r_torus+rmax)*cos(eta)
                        v(2) = yp(i) + (r_torus+rmax)*sin(eta)
                        v(3) = zp(i);
                        xp(i) = dot_product(v,D2(1:3,1))
                        yp(i) = dot_product(v,D2(1:3,2))
                        zp(i) = dot_product(v,D2(1:3,3))
                        v(1) = wxp(i);
                        v(2) = wyp(i);
                        v(3) = -wzp(i);
                        wxp(i) = dot_product(v,D2(1:3,1))
                        wyp(i) = dot_product(v,D2(1:3,2))
                        wzp(i) = dot_product(v,D2(1:3,3))
                    else
                        v(1) = xp(i)
                        v(2) = yp(i)
                        v(3) = zp(i)
                        xp(i) = dot_product(v,D4(1:3,1))
                        yp(i) = dot_product(v,D4(1:3,2))
                        zp(i) = dot_product(v,D4(1:3,3))
                        v(1) = wxp(i)
                        v(2) = wyp(i)
                        v(3) = -wzp(i)
                        wxp(i) = dot_product(v,D4(1:3,1))
                        wyp(i) = dot_product(v,D4(1:3,2))
                        wzp(i) = dot_product(v,D4(1:3,3))
                    end if
                    ind0 = ind0 + 1
                    if (mod(ind0-1,n_cpu) == my_rank) then
                        ind = ind + 1
                        if (ind .gt. np) then
                            write(*,*) 'something is wrong here: to many particles in init of second ring',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                        end if
                        vortex_particles(ind)%x(1) = xp(i) + (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(2) = yp(i) + (torus_offset(2)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%x(3) = zp(i) + (torus_offset(3)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
                        vortex_particles(ind)%data%alpha(1) = wxp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(2) = wyp(i)*volp(i)
                        vortex_particles(ind)%data%alpha(3) = wzp(i)*volp(i)
                    end if
                end do
            end do
            np = ind

            deallocate(xp,yp,zp,volp,wxp,wyp,wzp)

        case(3)  ! Sphere setup

            a = nine/five*sqrt(dble(n-1)/dble(n))
            j = 0
            do i = 1,n

                if(i.eq.1) then
                    cth = mone
                    sth = zero
                    cphi = one
                    sphi = zero
                    b = cphi
                    c = sphi
                elseif(i.eq.n) then
                    cth = one
                    sth = zero
                    cphi = one
                    sphi = zero
                else
                    cth = dble(2*i-n-1)/dble(n-1)
                    sth = two*(sqrt(dble(i-1)/dble(n-1))*sqrt(dble(n-i)/dble(n-1)))
                    s = a*sqrt(dble(n-1)/(dble(i-1)*dble(n-i)))
                    cphi = b*cos(s)-c*sin(s)
                    sphi = c*cos(s)+b*sin(s)
                    b = cphi
                    c = sphi
                endif
                if (mod(i+my_rank,n_cpu) == 0) then
                    if (j .gt. np-1) then
                        write(*,*) 'something is wrong here: to many particles in init',my_rank,j,np,i
                        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                    end if
                    j = j+1
                    vortex_particles(j)%x(1) = sth*cphi
                    vortex_particles(j)%x(2) = sth*sphi
                    vortex_particles(j)%x(3) = cth
                    vortex_particles(j)%data%alpha(1) = 0.3D01/(0.8D01*pi)*sth*sphi*h**2
                    vortex_particles(j)%data%alpha(2) = 0.3D01/(0.8D01*pi)*sth*(-cphi)*h**2
                    vortex_particles(j)%data%alpha(3) = 0.
                end if
            end do
            np = j

        case(4) ! Vortex wakes

            allocate(xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

            j = 0

            do k = 1, nc
                part_2d = 2*pi/(8*k)
                rc = (1+12*k**2)/(6*k)*rl

                do l = 1,8*k
                    j = j+1
                    xi1 = part_2d*(l-1)
                    xi2 = part_2d*l
                    xi = (xi2-xi1)/2+xi1
                    xp(j) = rc*cos(xi)
                    yp(j) = rc*sin(xi)
                    zp(j) = 0
                    volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
                    wxp(j) = 0.
                    wyp(j) = 0.
                    wzp(j) = g*exp(-(rc/rmax)**2)
                end do
            end do

            xp(ns) = 0.
            yp(ns) = 0.
            zp(ns) = 0.
            wxp(ns) = 0.
            wyp(ns) = 0.
            wzp(ns) = g
            volp(ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

            part_3d = 4*pi/nphi
            ind0 = 0
            ind = 0
            zp(1:ns) = zp(1:ns) - 2*pi - part_3d
            do k = 1,nphi

                zp(1:ns) = zp(1:ns)+part_3d

                do i = 1, ns
                    ind0 = ind0 + 1
                    if (mod(ind0-1,n_cpu) == my_rank) then
                        ind = ind + 1
                        if (ind .gt.np) then
                            write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                        end if
                        vortex_particles(ind)%x(1) = xp(i)
                        vortex_particles(ind)%x(2) = yp(i) - torus_offset(2)
                        vortex_particles(ind)%x(3) = zp(i)
                        vortex_particles(ind)%data%alpha(1) = 0.
                        vortex_particles(ind)%data%alpha(2) = 0.
                        vortex_particles(ind)%data%alpha(3) = -(exp(-(zp(i)-pi)**2)+exp(-(zp(i)+pi)**2))*wzp(i)*volp(i)
                        ind = ind + 1
                        vortex_particles(ind)%x(1) = xp(i)
                        vortex_particles(ind)%x(2) = yp(i) + torus_offset(2)
                        vortex_particles(ind)%x(3) = zp(i)
                        vortex_particles(ind)%data%alpha(1) = 0.
                        vortex_particles(ind)%data%alpha(2) = 0.
                        vortex_particles(ind)%data%alpha(3) = +(exp(-(zp(i)-pi)**2)+exp(-zp(i)**2)+exp(-(zp(i)+pi)**2))*wzp(i)*volp(i)
                    end if
                end do

            end do

            np = ind

        case(5)  ! Different wakes

            ind0 = 0
            ind = 0
            do i=1,ceiling(1.0/h)
                do j=1,ceiling(2*pi/h)
                    do k=1,nc
                        ind0 = ind0 + 1
                        if (mod(ind0+my_rank,n_cpu) == 0) then
                            ind = ind+1
                            if (ind .gt. np-1) then
                                write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n,n_cpu
                                call MPI_ABORT(MPI_COMM_WORLD,1, ierr)
                            end if
                            xt = (i-1)*h !+ h/2
                            yt = -pi + (j-1)*h !+ h/2
                            zt = -pi + (k-1)*h

                            vortex_particles(ind)%x(1) = xt + torus_offset(1)
                            vortex_particles(ind)%x(2) = yt + torus_offset(2)
                            vortex_particles(ind)%x(3) = zt + torus_offset(3)
                            vortex_particles(ind)%data%alpha(1) = 0.
                            vortex_particles(ind)%data%alpha(2) = 0.
                            vortex_particles(ind)%data%alpha(3) = -g/2*(1-tanh(yt)**2)*h**3 !* (exp(-zt**2/2)+exp(-(zt-pi/2)**2/2)+exp(-(zt+pi/2)**2/2))
                            ind = ind + 1
                            vortex_particles(ind)%x(1) = xt - torus_offset(1)
                            vortex_particles(ind)%x(2) = yt - torus_offset(2)
                            vortex_particles(ind)%x(3) = zt - torus_offset(3)
                            vortex_particles(ind)%data%alpha(1) = 0.
                            vortex_particles(ind)%data%alpha(2) = 0.
                            vortex_particles(ind)%data%alpha(3) = +g/2*(1-tanh(yt)**2)*h**3 !* (exp(-zt**2/2)+exp(-(zt-pi/2)**2/2)+exp(-(zt+pi/2)**2/2))
                        end if
                    end do
                end do
            end do
            np = ind


        case(98) ! Random cubic setup (for testing purpose only)

            j = 0

            do i = 1, n

                xt = 0.
                yt = 0.
                zt = 0.

                call par_rand(par_rand_res)
                xt = par_rand_res
                call par_rand(par_rand_res)
                yt = par_rand_res
                call par_rand(par_rand_res)
                zt = par_rand_res

                if (mod(i+my_rank,n_cpu) == 0) then
                    j = j+1
                    if (j .gt. np) then
                        write(*,*) 'something is wrong here: to many particles in init',my_rank,j,n
                        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
                    end if
                    vortex_particles(j)%x(1) = xt
                    vortex_particles(j)%x(2) = yt
                    vortex_particles(j)%x(3) = zt
                    call par_rand(par_rand_res)
                    vortex_particles(j)%data%alpha(1) = par_rand_res*h**3
                    call par_rand(par_rand_res)
                    vortex_particles(j)%data%alpha(2) = -par_rand_res*h**3
                    call par_rand(par_rand_res)
                    vortex_particles(j)%data%alpha(3) = par_rand_res*h**3
                end if
            end do
            np = j

        case(99) ! Read-in MPI checkpoints

            call read_in_checkpoint()

        end select config

        ! shrink vortex_particles to correct size
        ! reasoning: vortex_particles may have been allocated larger than required
        ! but may be queried for its size at a later stage
        if (np .ne. size(vortex_particles)) then
           block
              type(t_particle), allocatable :: temp_particles(:)
              call move_alloc(vortex_particles, temp_particles)
              vortex_particles = temp_particles(1:np)
           end block
        end if

        ! initial dump if we did not read in via MPI
        if (ispecial .ne. 99) then
            call dump(0,ts)
            call kick_out_particles()
            call reset_labels() ! works on vortex_particles
            vortex_particles(1:np)%work = 1.
        end if

    end subroutine special_start

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Update particle positions - using 2nd order RK
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_rk2(stage)

      use physvars
      integer, intent(in) :: stage   ! In which RK stage are we?
      integer(kind_particle) :: i

      do i=1,np

        if (stage == 1) then
            ! Euler predictor
            vortex_particles(i)%data%x_rk(1:3)      = vortex_particles(i)%x(1:3)
            vortex_particles(i)%data%alpha_rk(1:3)  = vortex_particles(i)%data%alpha(1:3)
            vortex_particles(i)%data%u_rk(1:3)      = vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%af_rk(1:3)     = vortex_particles(i)%results%af(1:3)
            vortex_particles(i)%x(1:3) = vortex_particles(i)%x(1:3) + dt*vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha(1:3) + dt*vortex_particles(i)%results%af(1:3)
        else
            ! Trapezoidal corrector
            vortex_particles(i)%x(1:3) = vortex_particles(i)%data%x_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%u_rk(1:3)+vortex_particles(i)%results%u(1:3))
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%af_rk(1:3)+vortex_particles(i)%results%af(1:3))
        end if

      end do

      call kick_out_particles()

    end subroutine push_rk2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   4th order remeshing using parallel sorting
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine remeshing()

        use physvars
        use module_timings
        use omp_lib
        use treevars, only: num_threads
        use mpi
        implicit none

        integer :: mesh_supp, ierr, i1, i2, i3, xtn, ytn, ztn, omp_thread_num
        integer(kind_particle) :: i, k
        integer(kind_particle) :: m_np, m_nppm, m_n, tmp
        real(kind_physics) :: frac, xt, yt, zt, axt, ayt, azt, wt
        real(kind_physics), allocatable :: mesh_offset(:)
        real(kind_physics), dimension(3) :: total_vort, total_vort_full_pre, total_vort_full_mid, total_vort_full_post
        type(t_particle_short), allocatable :: m_part(:)
        integer, parameter :: t_remesh_interpol = t_userdefined_first + 2
        integer, parameter :: t_remesh_sort = t_userdefined_first + 3

        total_vort = 0.
        do i = 1,np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
        end do
        call MPI_ALLREDUCE(total_vort,total_vort_full_pre,3,MPI_KIND_PHYSICS,MPI_SUM,MPI_COMM_WORLD,ierr)

        !! Define mesh structure and support, currently only none or 4 are possible
        mesh_supp = 4 !TODO: make this a parameter (and adapt remeshing itself to it)
        allocate(mesh_offset(mesh_supp))
        mesh_offset = (/-0.15D+01, -0.5D+00, 0.5D+00, 0.15D+01/)

        ! Define array limits for new particles on the mesh
        m_np = np*mesh_supp**3
        m_n  = n *mesh_supp**3
        m_nppm = ceiling(1.05*max(1.0*m_n/n_cpu,1.0*m_np)) ! allow 5% fluctuation, just a safety factor, since we`ll kick out particles beforehand

        ! TODO: Define global max. #particles (do we need this?)
        call MPI_ALLREDUCE(m_nppm,tmp,1,MPI_KIND_PARTICLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        m_nppm = tmp

        allocate(m_part(m_nppm),STAT=ierr)

        if (ierr .ne. 0) then
            write(*,*) 'something is wrong here: remeshing allocation failed', my_rank, ierr
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
        end if

        call timer_start(t_remesh_interpol)

        ! Set number of openmp threads to the same number as pthreads used in the walk
        !$ call omp_set_num_threads(num_threads)

        ! Inform the user that openmp is used, and with how many threads
        !$OMP PARALLEL PRIVATE(omp_thread_num)
        !$ omp_thread_num = OMP_GET_THREAD_NUM()
        !$ if( (my_rank .eq. 0) .and. (omp_thread_num .eq. 0) ) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads.'
        !$OMP END PARALLEL

        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,xt,yt,zt,xtn,ytn,ztn,axt,ayt,azt,wt,k,frac)
        do i=1,np
            xt  = vortex_particles(i)%x(1)
            yt  = vortex_particles(i)%x(2)
            zt  = vortex_particles(i)%x(3)
            xtn = nint(xt/m_h)
            ytn = nint(yt/m_h)
            ztn = nint(zt/m_h)
            axt = vortex_particles(i)%data%alpha(1)
            ayt = vortex_particles(i)%data%alpha(2)
            azt = vortex_particles(i)%data%alpha(3)
            wt  = vortex_particles(i)%work
            ! For each particle i define its neighbor mesh points and project onto them
            k = (i-1)*mesh_supp**3
            do i1 = 1,mesh_supp
                do i2 = 1,mesh_supp
                    do i3 = 1,mesh_supp
                        k=k+1
                        m_part(k)%x(1) = m_h*(xtn + mesh_offset(i1))
                        m_part(k)%x(2) = m_h*(ytn + mesh_offset(i2))
                        m_part(k)%x(3) = m_h*(ztn + mesh_offset(i3))
                        frac = ip_kernel(dabs(xt-m_part(k)%x(1))/m_h,kernel_c)*&
                               ip_kernel(dabs(yt-m_part(k)%x(2))/m_h,kernel_c)*&
                               ip_kernel(dabs(zt-m_part(k)%x(3))/m_h,kernel_c)
                        m_part(k)%data%alpha(1) = frac*axt
                        m_part(k)%data%alpha(2) = frac*ayt
                        m_part(k)%data%alpha(3) = frac*azt
                        m_part(k)%work = frac*wt
                        !m_part(k)%data%x_rk = 0.
                        !m_part(k)%data%alpha_rk = 0.
                        !m_part(k)%data%u_rk(1:3) = 0.
                        !m_part(k)%data%af_rk(1:3) = 0.
                        !m_part(k)%results%u(1:3) = 0.
                        !m_part(k)%results%af(1:3) = 0.
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        ! Reset the number of openmp threads to 1.
        !$ call omp_set_num_threads(1)

        call timer_stop(t_remesh_interpol)

        deallocate(vortex_particles)

        call timer_start(t_remesh_sort)
        call sort_remesh(m_part, m_np, m_nppm)
        call timer_stop(t_remesh_sort)

        np = m_np
        allocate(vortex_particles(np))
        vortex_particles(1:np) = m_part(1:m_np)

        deallocate(mesh_offset, m_part)

        total_vort = 0.
        do i = 1,np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
        end do
        call MPI_ALLREDUCE(total_vort,total_vort_full_mid,3,MPI_KIND_PHYSICS,MPI_SUM,MPI_COMM_WORLD,ierr)

        !call kick_out_particles()
        call MPI_ALLREDUCE(np,n,1,MPI_KIND_PARTICLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        if (1.25*n/n_cpu .lt. np) then
            write(*,*) 'warning, rank',my_rank,' appears to be heavily imbalanced:',1.0*np/(1.0*n/n_cpu)
        end if
        
        call reset_labels() ! works on vortex_particles

        total_vort = 0.
        do i = 1,np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
        end do
        call MPI_ALLREDUCE(total_vort,total_vort_full_post,3,MPI_KIND_PHYSICS,MPI_SUM,MPI_COMM_WORLD,ierr)

        if (my_rank == 0) then
            write(*,*) '   Vorticity before remeshing (x,y,z,norm2):', total_vort_full_pre, sqrt(dot_product(total_vort_full_pre,total_vort_full_pre))
            write(*,*) 'Vorticity before pop. control (x,y,z,norm2):', total_vort_full_mid, sqrt(dot_product(total_vort_full_mid,total_vort_full_mid))
            write(*,*) '    Vorticity after remeshing (x,y,z,norm2):', total_vort_full_post, sqrt(dot_product(total_vort_full_post,total_vort_full_post))
        end if

    end subroutine remeshing


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Interpolation function for remeshing
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ip_kernel(dist,c)
       use physvars, only: pi

        implicit none

        real(kind_physics), intent(in) :: dist, c
        real(kind_physics) :: ip_kernel

        ip_kernel = 0.0D+00
        if (dist .le. 1) then
            ip_kernel = 0.1D01 - (0.25D+01)*dist**2 + (0.15D+01)*dist**3 - c**2*(0.2D+01-(0.9D+01)*dist**2+(0.6D+01)*dist**3)
        end if
        if ((dist .gt. 1) .and. (dist .le. 2)) then
            ip_kernel = (0.5D+00)*((0.1D+01)-dist-(0.2D+01)*c**2+(0.4D+01)*c**2*dist)*((0.2D+01)-dist)**2
        end if

    end function ip_kernel


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Sorting function for remeshing
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sort_remesh(particles, m_np, m_nppm)

        use physvars
        use treevars, only : maxlevel, idim
        use module_sort, only : sort
        use mpi
        implicit none

        type(t_particle_short), intent(inout) :: particles(*)
        integer(kind_particle), intent(inout) :: m_np
        integer(kind_particle), intent(in) :: m_nppm

        integer :: ierr
        integer(kind_particle) :: i, j, k
        type(t_particle_short) :: bound_parts_loc(2), bound_parts(0:2*n_cpu-1), ship_parts(m_np), get_parts(m_nppm)
        integer :: prev, next, nbits
!TODO should this be kind_default or really kind_particle?!?!?
        integer :: irnkl2(m_nppm), indxl(m_np), irnkl(m_nppm)
        real(kind_physics) :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local, s, local_work(m_nppm)
        real(kind_physics) :: xboxsize, yboxsize, zboxsize, boxsize, xmax, xmin, ymax, ymin, zmax, zmin, thresh2
        integer(kind_particle) :: ix(m_np), iy(m_np), iz(m_np)
        integer(kind_key) :: sorted_keys(m_nppm), local_keys(m_nppm)
        integer :: fposts(n_cpu+1),gposts(n_cpu+1),islen(n_cpu),irlen(n_cpu)
        integer(kind_particle) :: npnew, npold
        
        integer(kind_key) :: iplace
        
        interface
           subroutine slsort_keys(nin, nmax, keys, workload, balance_weight, max_imbalance, nout, indxl, &
                                  irnkl, scounts, rcounts, sdispls, rdispls, keys2, irnkl2, size, rank, comm)
              use module_pepc_kinds
              integer(kind_particle), intent(in) :: nin
              integer(kind_particle), intent(in) :: nmax
              integer(kind_key), intent(inout) :: keys(*)
              real*8, intent(inout) :: workload(*)
              integer(kind_default), intent(in) :: balance_weight
              real*8, intent(in) :: max_imbalance
              integer(kind_particle), intent(out) :: nout
              integer(kind_default), intent(out) :: indxl(*), irnkl(*), scounts(*), rcounts(*), sdispls(*), rdispls(*)
              integer(kind_key), intent(out) :: keys2(*)
              integer(kind_default), intent(out) :: irnkl2(*)
              integer(kind_pe), intent(in) :: size, rank
              integer(kind_default), intent(in) :: comm
           end subroutine slsort_keys
        end interface
        
        iplace = 2_8**(idim * maxlevel)

        indxl = 0
        irnkl = 0
        islen = 0
        irlen = 0

        ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
        if (my_rank == 0) then
            prev = n_cpu - 1
        else
            prev = my_rank-1
        endif

        if (my_rank == n_cpu-1 ) then
            next = 0
        else
            next = my_rank+1
        endif

        xmin_local = minval(particles(1:m_np)%x(1))
        xmax_local = maxval(particles(1:m_np)%x(1))
        ymin_local = minval(particles(1:m_np)%x(2))
        ymax_local = maxval(particles(1:m_np)%x(2))
        zmin_local = minval(particles(1:m_np)%x(3))
        zmax_local = maxval(particles(1:m_np)%x(3))

        ! Find global limits
        call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_KIND_PHYSICS, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_KIND_PHYSICS, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymin_local, ymin, 1, MPI_KIND_PHYSICS, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymax_local, ymax, 1, MPI_KIND_PHYSICS, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_KIND_PHYSICS, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_KIND_PHYSICS, MPI_MAX,  MPI_COMM_WORLD, ierr )

        xboxsize = xmax-xmin
        yboxsize = ymax-ymin
        zboxsize = zmax-zmin

        ! Safety margin - put buffer region around particles
        xmax = xmax + xboxsize/10000
        xmin = xmin - xboxsize/10000
        ymax = ymax + yboxsize/10000
        ymin = ymin - yboxsize/10000
        zmax = zmax + zboxsize/10000
        zmin = zmin - zboxsize/10000

        boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

        s=boxsize/2**maxlevel       ! refinement length

        !! Start key generation
        ! TODO: Use module_spacefilling here, problem: this module uses treevars variables, but we do not (e.g. npp vs. m_np)

        ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
        ix(1:m_np) = ( particles(1:m_np)%x(1) - xmin )/s           ! partial keys
        iy(1:m_np) = ( particles(1:m_np)%x(2) - ymin )/s           !
        iz(1:m_np) = ( particles(1:m_np)%x(3) - zmin )/s

        ! construct keys by interleaving coord bits and add placeholder bit
        ! - note use of 64-bit constants to ensure correct arithmetic
        nbits = maxlevel+1
        do j = 1,m_np
            local_keys(j) = iplace
            do i=0,nbits-1
                local_keys(j) = local_keys(j) &
                + 8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) )
            end do
        end do

        ! Sort particles locally according to local_keys
        particles(1:m_np)%key = local_keys(1:m_np)
        call sort(local_keys(1:m_np), indxl(1:m_np))
        ship_parts(1:m_np) = particles(indxl(1:m_np))
        particles(1:m_np) = ship_parts(1:m_np)

        ! Check if sort finished locally and find inner doublets
        k = 0
        do i=2,m_np
            if (particles(i)%key .lt. particles(i-1)%key) then
                write (*,'(a,i5,2z20)') 'Sorting locally failed: i,particles(i)%key,particles(i-1)%key=',i,particles(i)%key,particles(i-1)%key
                call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
            else if (particles(i)%key == particles(i-1)%key) then
                particles(i)%data%alpha(1) = particles(i)%data%alpha(1) + particles(i-1)%data%alpha(1)
                particles(i)%data%alpha(2) = particles(i)%data%alpha(2) + particles(i-1)%data%alpha(2)
                particles(i)%data%alpha(3) = particles(i)%data%alpha(3) + particles(i-1)%data%alpha(3)
                particles(i)%work = particles(i)%work + particles(i-1)%work
            else
                k = k+1
                particles(k) = particles(i-1)
            end if
        end do
        ! Last original particle will always be transfered to next/final active particle
        k = k + 1
        particles(k) = particles(m_np)

        m_np = k

        ! Use Parallel Sort by Parallel Search (SLSORT by M. Hofmann, Chemnitz)
        indxl = 0
        npold = m_np
        npnew = m_np
        local_keys(1:npold) = particles(1:npold)%key
        local_work(1:npold) = particles(1:npold)%work

        call slsort_keys(npold, m_nppm, local_keys, local_work, 0, 0.05D0, npnew, indxl, irnkl, islen, irlen, fposts, gposts, &
                         sorted_keys, irnkl2, n_cpu, my_rank, MPI_COMM_WORLD)

        ! Permute particles according to arrays from slsort
        m_np = npnew
        ship_parts(1:npold) = particles(indxl(1:npold))
        call MPI_ALLTOALLV(ship_parts, islen, fposts, MPI_TYPE_PARTICLE_SHORT_sca, &
        get_parts, irlen, gposts, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD,ierr)
        particles(irnkl(1:m_np)) = get_parts(1:m_np)
        particles(1:m_np)%key = sorted_keys(1:m_np)

        ! Check if sort finished and find inner doublets
        k = 0
        do i=2,m_np
            if (particles(i)%key .lt. particles(i-1)%key) then
                write (*,'(a,i5,2z20)') 'Sorting globally failed: i,particles(i)%key,particles(i-1)%key=',i,particles(i)%key,particles(i-1)%key
                call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
            else if (particles(i)%key == particles(i-1)%key) then
                particles(i)%data%alpha(1) = particles(i)%data%alpha(1) + particles(i-1)%data%alpha(1)
                particles(i)%data%alpha(2) = particles(i)%data%alpha(2) + particles(i-1)%data%alpha(2)
                particles(i)%data%alpha(3) = particles(i)%data%alpha(3) + particles(i-1)%data%alpha(3)
                particles(i)%work = particles(i)%work + particles(i-1)%work
            else
                k = k+1
                particles(k) = particles(i-1)
            end if
        end do
        ! Last original particle will always be transfered to next/final active particle
        k = k + 1
        particles(k) = particles(m_np)

        if (k == 0) then
            write(*,*) 'something is wrong here: no particles left after compression',my_rank
            call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
        end if
        if (k == 1) then
            write(*,*) 'Warning: Only one particle left on this process, candidate for m_np=0',my_rank
        end if

        bound_parts_loc(1) = particles(1)
        bound_parts_loc(2) = particles(k)
        call MPI_ALLGATHER(bound_parts_loc, 2, MPI_TYPE_PARTICLE_SHORT_sca, &
        bound_parts, 2, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD, ierr)

        ! Eliminate right boundary, if doublet with at least on neighbor
        if (bound_parts(2*my_rank+1)%key .eq. bound_parts(2*next)%key) then
            write(*,*) my_rank, 'is looking towards', next, 'for elimination'
            k = k-1
        end if

        ! Accumulate boundary particles, looking left for doublets
        i = prev
        do while (bound_parts(2*my_rank)%key .eq. bound_parts(2*i+1)%key)
            write(*,*) my_rank, 'is looking towards', i, 'for accumulation'
            particles(1)%data%alpha(1) = particles(1)%data%alpha(1) + bound_parts(2*i+1)%data%alpha(1)
            particles(1)%data%alpha(2) = particles(1)%data%alpha(2) + bound_parts(2*i+1)%data%alpha(2)
            particles(1)%data%alpha(3) = particles(1)%data%alpha(3) + bound_parts(2*i+1)%data%alpha(3)
            particles(1)%work = particles(1)%work + bound_parts(2*i+1)%work
            i = i-1
            if (i.lt.0) exit
        end do

        m_np = k

        ! Kick out particles (cannot use subroutinee here, since we work on a temp-array) 
        thresh2 = thresh**2
        k = 0
        do i = 1,m_np
            if (dot_product(particles(i)%data%alpha,particles(i)%data%alpha) .gt. thresh2) then
                k = k+1
                particles(k) = particles(i)
            end if
        end do
        m_np = k
        
        ! Rebalance: Use Parallel Sort by Parallel Search (SLSORT by M. Hofmann, Chemnitz)
        indxl = 0
        npold = m_np
        npnew = m_np
        local_keys(1:npold) = particles(1:npold)%key
        local_work(1:npold) = particles(1:npold)%work

        call slsort_keys(npold, m_nppm, local_keys, local_work, 0, 0.05D0, npnew, indxl, irnkl, islen, irlen, fposts, gposts, &
                         sorted_keys, irnkl2, n_cpu, my_rank, MPI_COMM_WORLD)

        ! Permute particles according to arrays from slsort
        m_np = npnew
        ship_parts(1:npold) = particles(indxl(1:npold))
        call MPI_ALLTOALLV(ship_parts, islen, fposts, MPI_TYPE_PARTICLE_SHORT_sca, &
        get_parts, irlen, gposts, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD,ierr)
        particles(irnkl(1:m_np)) = get_parts(1:m_np)
        particles(1:m_np)%key = sorted_keys(1:m_np)
        particles(1:m_np)%work = 1. !TODO: is this elegant? 

    end subroutine sort_remesh


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Kick out particles according to threshold (compare with |alpha|)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kick_out_particles()

        use physvars
        use mpi
        implicit none

        integer :: ierr
        integer(kind_particle) :: i, k
        real(kind_physics) :: thresh2
        type(t_particle), allocatable :: temp_particles(:)

        ! kick out particles with vorticity magnitude below threshold
        thresh2 = thresh**2
        k = 0
        do i = 1,np
            if (dot_product(vortex_particles(i)%data%alpha,vortex_particles(i)%data%alpha) .gt. thresh2) then
                k = k+1
                vortex_particles(k) = vortex_particles(i)
            end if
        end do
        np = k
        call MPI_ALLREDUCE(np,n,1,MPI_KIND_PARTICLE,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! shrink vortex_particles to new size
        call move_alloc(vortex_particles, temp_particles)
        vortex_particles = temp_particles(1:np)

        if (1.25*n/n_cpu .lt. np) then
            write(*,*) 'warning, rank',my_rank,' appears to be heavily imbalanced:',1.0*np/(1.0*n/n_cpu)
        end if

    end subroutine kick_out_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Initialize/reset labels after manipulation, esp. after remeshing
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine reset_labels()

        use physvars
        use mpi
        implicit none

        integer(kind_particle) :: i, nscan
        integer :: ierr

        ! Define valid labels for all local particles
        nscan = 0
        call MPI_SCAN(np,nscan,1,MPI_KIND_PARTICLE,MPI_SUM,MPI_COMM_WORLD,ierr)
        nscan = nscan-np
        do i=1,np
            vortex_particles(i)%label = nscan+i
        end do

    end subroutine reset_labels

    subroutine direct_sum(np_local,particles,results,my_rank,n_cpu)

        use module_pepc_types
        use module_interaction_specific_types, only: t_particle_results
        use module_directsum
        use mpi
        implicit none

        integer(kind_particle), intent(in) :: np_local
        integer, intent(in) :: my_rank, n_cpu
        type(t_particle), intent(in) :: particles(1:np_local)
        type(t_particle_results), intent(out) :: results(1:np_local)

        integer(kind_particle) :: i
        type(t_particle_results), allocatable :: directresults(:)
        integer(kind_particle) :: indices(1:np_local)

        if (my_rank==0) write(*,'("PEPC-V | ", a)') 'Starting direct summation ...'

        do i = 1,np_local
            indices(i) = i
        end do

        allocate(directresults(int(np_local, kind=kind_particle)))
        call directforce(particles, indices, int(np_local, kind=kind_particle), directresults, MPI_COMM_WORLD)
        results(1:np_local) = directresults(1:np_local)

        deallocate(directresults)

        if (my_rank==0) write(*,'("PEPC-V | ", a)') '                          ... done!'

    end subroutine direct_sum

    ! portable random number generator, see numerical recipes
    ! check for the random numbers:
    ! the first numbers should be 0.2853809, 0.2533582 and 0.2533582
    subroutine par_rand(res, iseed)

        !use physvars
        implicit none

        integer :: idum, idum2, iy, j, k
        integer :: iv(32)

        real, intent(inout) :: res
        integer, intent(in), optional :: iseed

        integer :: IM1, IM2, IMM1, IA1, IA2, IQ1, IQ2, IR1, IR2, NTAB, NDIV
        real    :: AM, RNMX

        save

        data idum, idum2 /-1, 123456789/

        IM1 = 2147483563
        IM2 = 2147483399
        AM  = 1.0/IM1
        IMM1 = IM1-1
        IA1 = 40014
        IA2 = 40692
        IQ1 = 53668
        IQ2 = 52774
        IR1 = 12211
        IR2 = 3791
        NTAB = 32
        NDIV = 1+IMM1/NTAB
        RNMX = 1.0 - 1.2e-7

        if (idum < 0) then

            if (present(iseed)) then
                idum = iseed
            else
                idum = 1
            endif

            idum2 = idum

            do j = NTAB+7,0,-1
                k = idum/IQ1
                idum = IA1 * (idum-k*IQ1) - k*IR1
                if (idum < 0 ) idum = idum + IM1

                if (j<NTAB) iv(j+1) = idum

            end do
            iy = iv(1)

        end if

        k = idum/IQ1
        idum = IA1 * (idum-k*IQ1) - k*IR1
        if (idum < 0) idum = idum + IM1

        k = idum2/IQ2
        idum2 = IA2 * (idum2-k*IQ2) - k*IR2
        if (idum2 < 0) idum2 = idum2 + IM2

        j = iy/NDIV + 1
        iy = iv(j)-idum2
        iv(j) = idum

        if (iy < 1) iy = iy + IMM1
        res = AM*iy
        if (res > RNMX) res = RNMX

    end subroutine par_rand


end module manipulate_particles
