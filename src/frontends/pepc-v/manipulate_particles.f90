module manipulate_particles

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
        implicit none
        include 'mpif.h'

        integer :: j, k, ind, ind0, i, m, nscan, ierr, l
        real*8 :: part_2d, rc, xi1, xi2, xi, part_3d, eta1, eta2, eta, v(3), thresh2
        real*8, dimension(3,3) :: D1, D2, D3, D4   !< rotation matrices for ring setup
        real*8, allocatable :: xp(:), yp(:), zp(:), volp(:), wxp(:), wyp(:), wzp(:)  !< helper arrays for ring setups

        ! weird helper variables for sphere setup (ask Holger Dachsel)
        real*8 a,b,c,cth,sth,cphi,sphi,s
        real*8 zero
        parameter(zero=0.d0)
        real*8 one
        parameter(one=1.d0)
        real*8 mone
        parameter(mone=-one)
        real*8 two
        parameter(two=2.d0)
        real*8 five
        parameter(five=5.d0)
        real*8 nine
        parameter(nine=9.d0)


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
                        if (ind .gt.np-1) then
                            write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,ierr)
                            stop
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
                        if (ind .gt. np-1) then
                            write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,ierr)
                            stop
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
                        if (ind .gt. np-1) then
                            write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,np,n
                            call MPI_ABORT(MPI_COMM_WORLD,ierr)
                            stop
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
                        call MPI_ABORT(MPI_COMM_WORLD,ierr)
                        stop
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

        case(99) ! Read-in MPI checkpoints

            call read_in_checkpoint()

        end select config

        call kick_out_particles()
        call reset_labels()

        vortex_particles(1:np)%work = 0.

        ! initial dump if we did not read in via MPI
        if (ispecial .ne. 99) call dump(0,ts)

    end subroutine special_start

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Update particle positions - using 2nd order RK
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine push_rk2(stage)

      use physvars
      integer, intent(in) :: stage   ! In which RK stage are we?
      integer :: i

      do i=1,np

        if (stage == 1) then
            ! Euler predictor
            vortex_particles(i)%data%u_rk(1:3)  = vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%af_rk(1:3) = vortex_particles(i)%results%af(1:3)
            vortex_particles(i)%x(1:3) = vortex_particles(i)%x(1:3) + dt*vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha(1:3) + dt*vortex_particles(i)%results%af(1:3)
        else
            ! Trapezoidal corrector
            vortex_particles(i)%x(1:3) = vortex_particles(i)%x(1:3)-dt*vortex_particles(i)%data%u_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%u_rk(1:3)+vortex_particles(i)%results%u(1:3))
            vortex_particles(i)%data%alpha(1:3) = vortex_particles(i)%data%alpha(1:3)-dt*vortex_particles(i)%data%af_rk(1:3) + 0.5*dt*(vortex_particles(i)%data%af_rk(1:3)+vortex_particles(i)%results%af(1:3))
        end if

      end do

      call kick_out_particles()

    end subroutine push_rk2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Kick out particles according to threshold (compare with |alpha|)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kick_out_particles()

        use physvars
        implicit none
        include 'mpif.h'

        integer :: i, k, ierr, nscan
        real*8 :: thresh2

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
        call MPI_ALLREDUCE(np,n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    end subroutine kick_out_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Initialize/reset labels after manipulation, esp. after remeshing
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine reset_labels()

        use physvars
        implicit none
        include 'mpif.h'

        integer :: i, ierr, nscan

        ! Define valid labels for all local particles
        nscan = 0
        call MPI_SCAN(np,nscan,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        nscan = nscan-np
        do i=1,np
            vortex_particles(i)%label = nscan+i
        end do

    end subroutine reset_labels


end module manipulate_particles
