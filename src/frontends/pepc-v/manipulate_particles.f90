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

        vortex_particles(1:np)%work = 1.

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
    !>   4th order remeshing using parallel sorting
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine remeshing(step)

        use physvars
        implicit none
        include 'mpif.h'

        integer, intent(in) :: step
        integer :: mesh_supp, m_np, m_nppm, tmp, ierr, k, i, i1, i2, i3
        real*8 :: frac
        real*8, allocatable :: mesh_offset(:)
        real*8, dimension(3) :: total_vort, total_vort_full_pre, total_vort_full_post

        type(t_particle), allocatable :: m_part(:)

        if ((rem_freq .gt. 0) .and. (mod(step,rem_freq)==0)) then

            total_vort = 0.
            do i = 1,np
                total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
            end do
            call MPI_ALLREDUCE(total_vort,total_vort_full_pre,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

            !! Define mesh structure and support, currently only none or 4 are possible
            mesh_supp = 4 !TODO: make this a parameter (and adapt remeshing itself to it)
            allocate(mesh_offset(mesh_supp))
            mesh_offset = (/-0.15D+01, -0.5D+00, 0.5D+00, 0.15D+01/)

            ! Define array limits for new particles on the mesh
            m_np = np*mesh_supp**3
            m_nppm = int(1.1*max(m_np,1000)) ! allow 10% fluctuation

            ! TODO: Define global max. #particles (do we need this?)
            call MPI_ALLREDUCE(m_nppm,tmp,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
            m_nppm = tmp

            allocate(m_part(m_nppm),STAT=ierr)

            if (ierr .ne. 0) then
                write(*,*) 'something is wrong here: remeshing allocation failed', my_rank, ierr
                call MPI_ABORT(MPI_COMM_WORLD,ierr)
                stop
            end if

            !! Start required remeshing, i.e. project particles onto mesh
            k=0
            do i=1,np
                ! For each particle i define its neighbor mesh points and project onto them
                do i1 = 1,mesh_supp
                    do i2 = 1,mesh_supp
                        do i3 = 1,mesh_supp
                            k=k+1
                            m_part(k)%x(1) = m_h*(nint(vortex_particles(i)%x(1)/m_h) + mesh_offset(i1))
                            m_part(k)%x(2) = m_h*(nint(vortex_particles(i)%x(2)/m_h) + mesh_offset(i2))
                            m_part(k)%x(3) = m_h*(nint(vortex_particles(i)%x(3)/m_h) + mesh_offset(i3))
                            frac = ip_kernel(dabs(vortex_particles(i)%x(1)-m_part(k)%x(1))/m_h,kernel_c)*&
                                   ip_kernel(dabs(vortex_particles(i)%x(2)-m_part(k)%x(2))/m_h,kernel_c)*&
                                   ip_kernel(dabs(vortex_particles(i)%x(3)-m_part(k)%x(3))/m_h,kernel_c)
                            m_part(k)%data%alpha(1) = frac*vortex_particles(i)%data%alpha(1)
                            m_part(k)%data%alpha(2) = frac*vortex_particles(i)%data%alpha(2)
                            m_part(k)%data%alpha(3) = frac*vortex_particles(i)%data%alpha(3)
                            m_part(k)%work = vortex_particles(i)%work
                            m_part(k)%data%u_rk(1:3) = 0.
                            m_part(k)%data%af_rk(1:3) = 0.
                        end do
                    end do
                end do
            end do
            m_np = k

            deallocate(vortex_particles)

            call sort_remesh(m_part, m_np, m_nppm)

            np = m_np
            allocate(vortex_particles(np))
            vortex_particles(1:np) = m_part(1:m_np)

            deallocate(mesh_offset, m_part)

            call kick_out_particles()
            call reset_labels()

            total_vort = 0.
            do i = 1,np
                total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
            end do
            call MPI_ALLREDUCE(total_vort,total_vort_full_post,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

            if (my_rank == 0) then
                write(*,*) 'Vorticity before remeshing (x,y,z,norm2):', total_vort_full_pre, sqrt(dot_product(total_vort_full_pre,total_vort_full_pre))
                write(*,*) '   Vorticity after merging (x,y,z,norm2):', total_vort_full_post, sqrt(dot_product(total_vort_full_post,total_vort_full_post))
            end if

        end if

    end subroutine remeshing

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Interpolation function for remeshing
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ip_kernel(dist,c)

        implicit none

        real*8, parameter :: pi=0.3141592654D+01
        real*8, intent(in) :: dist, c
        real*8 :: ip_kernel

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
    subroutine sort_remesh(m_part, m_np, m_nppm)

        use physvars
        use treevars, only : iplace, nlev
        implicit none
        include 'mpif.h'

        type (t_particle),intent(inout) :: m_part(*)
        integer, intent(inout) :: m_np
        integer, intent(in) :: m_nppm

        integer :: i,j,ierr,k, loc_hconst, coll
        type (t_particle) :: bound_parts_loc(2), bound_parts(0:2*n_cpu-1), ship_parts(m_nppm), get_parts(m_nppm)
        integer :: prev, next, dcount_glo, dcount, nbits, handle(4), sort_step, Nscan, status(MPI_STATUS_SIZE)
        integer :: irnkl2(m_nppm), indxl(m_nppm), irnkl(m_nppm), kick_part(m_nppm,2)
        real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local, s, imba, local_work(m_nppm)
        real*8 :: xboxsize, yboxsize, zboxsize, boxsize, xmax, xmin, ymax, ymin, zmax, zmin
        integer*8 :: tmp_prev, tmp_next, ix(m_np), iy(m_np), iz(m_np), sorted_keys(m_nppm), local_keys(m_nppm)
        integer :: fposts(n_cpu+1),gposts(n_cpu+1),islen(n_cpu),irlen(n_cpu)
        integer :: load_balance, npnew, npold

        interface
            subroutine slsort_keys(nin,nmax,keys,workload,balance_weight,max_imbalance,nout,indxl,irnkl,scounts,rcounts,sdispls,rdispls,keys2,irnkl2,size,rank)
                integer,intent(in) :: nin,nmax,balance_weight,size,rank
                real*8,intent(in) :: max_imbalance
                integer,intent(out) :: nout,indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*),irnkl2(*)
                integer*8,intent(out) :: keys2(*)
                integer*8,intent(inout) :: keys(*)
                real*8,intent(inout) :: workload(*)
            end subroutine slsort_keys
        end interface

        !! Define new particle counts

        kick_part = 0
        indxl = 0
        irnkl = 0
        sort_step = 0
        dcount_glo = -1
        dcount = 0

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

        xmin_local = minval(m_part(1:m_np)%x(1))
        xmax_local = maxval(m_part(1:m_np)%x(1))
        ymin_local = minval(m_part(1:m_np)%x(2))
        ymax_local = maxval(m_part(1:m_np)%x(2))
        zmin_local = minval(m_part(1:m_np)%x(3))
        zmax_local = maxval(m_part(1:m_np)%x(3))

        ! Find global limits
        call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymin_local, ymin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymax_local, ymax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )

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

        s=boxsize/2**nlev       ! refinement length

        ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

        ix(1:m_np) = ( m_part(1:m_np)%x(1) - xmin )/s           ! partial keys
        iy(1:m_np) = ( m_part(1:m_np)%x(2) - ymin )/s           !
        iz(1:m_np) = ( m_part(1:m_np)%x(3) - zmin )/s

        ! construct keys by interleaving coord bits and add placeholder bit
        ! - note use of 64-bit constants to ensure correct arithmetic
        nbits = nlev+1
        do j = 1,m_np
            m_part(j)%key = iplace
            do i=0,nbits-1
                m_part(j)%key = m_part(j)%key &
                + 8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) )
            end do
        end do

        do while (dcount_glo .ne. 0) ! BEGIN: Loop over soring algorithm, as long as we find doublets
                                      !        (should stop after two steps!)
            sort_step = sort_step +1
            if (sort_step .ge. 3) then
                write(*,*) 'something is wrong here: still finding doublets',my_rank,dcount_glo
                call MPI_ABORT(MPI_COMM_WORLD,ierr)
                stop
            end if

            ! If we have doublets...
            if (sort_step .gt. 1) then
                ! Check for local doublets (kick_part(.,1) == 1), accumulate w and eliminate by in-place copy
                i = 0
                k = 0
                do while (i .lt. m_np)
                    i = i+1
                    if (kick_part(i,1) == 1) then
                        j = kick_part(i,2)
                        if (j .gt. m_np) then
                            write(*,*) 'something is wrong here: accumulating to bad address',my_rank,j,m_np
                            call MPI_ABORT(MPI_COMM_WORLD,ierr)
                            stop
                        end if
                        m_part(j)%data%alpha(1) = m_part(j)%data%alpha(1) + m_part(i)%data%alpha(1)
                        m_part(j)%data%alpha(2) = m_part(j)%data%alpha(2) + m_part(i)%data%alpha(2)
                        m_part(j)%data%alpha(3) = m_part(j)%data%alpha(3) + m_part(i)%data%alpha(3)
                        m_part(j)%work = m_part(j)%work + m_part(i)%work
                     else
                        k = k+1
                        m_part(k) = m_part(i)
                    end if
                end do

                if (k == 0) then
                    write(*,*) 'something is wrong here: no particles left after compression',my_rank
                    call MPI_ABORT(MPI_COMM_WORLD,ierr)
                    stop
                end if
                if (k == 1) then
                    write(*,*) 'Warning: Only one particle left on this process, candidate for m_np=0',my_rank
                end if

                !TODO: Only necessary if boundary doublets detected previously (flagging?)
                ! Gather boundary particles from every process
                bound_parts_loc(1) = m_part(1)
                bound_parts_loc(2) = m_part(k)
                call MPI_ALLGATHER(bound_parts_loc, 2, MPI_TYPE_PARTICLE, &
                bound_parts, 2, MPI_TYPE_PARTICLE, MPI_COMM_WORLD, ierr)

                ! Eliminate right boundary, if doublet with at least on neighbor
                if (bound_parts(2*my_rank+1)%key .eq. bound_parts(2*next)%key) then
                    !write(*,*) my_rank, 'is looking towards', next, 'for elimination'
                    k = k-1
                end if
                m_np = k

                ! Accumulate boundary particles, looking left for doublets
                i = prev
                do while (bound_parts(2*my_rank)%key .eq. bound_parts(2*i+1)%key)
                    !write(*,*) my_rank, 'is looking towards', i, 'for accumulation'
                    m_part(1)%data%alpha(1) = m_part(1)%data%alpha(1) + bound_parts(2*i+1)%data%alpha(1)
                    m_part(1)%data%alpha(2) = m_part(1)%data%alpha(2) + bound_parts(2*i+1)%data%alpha(2)
                    m_part(1)%data%alpha(3) = m_part(1)%data%alpha(3) + bound_parts(2*i+1)%data%alpha(3)
                    m_part(1)%work = m_part(1)%work + bound_parts(2*i+1)%work
                    i = i-1
                    if (i.lt.0) exit
                end do

            end if

            kick_part = 0
            dcount = 0

            ! Use Parallel Sort by Parallel Search (SLSORT by M. Hofmann, Chemnitz)
            npold = m_np
            npnew = m_np
            load_balance = 1
            imba = 0.05
            local_keys(1:npold) = m_part(1:npold)%key
            local_work(1:npold) = m_part(1:npold)%work
            call slsort_keys(npold,m_nppm,local_keys,local_work,load_balance,imba,&
                             npnew,indxl,irnkl,islen,irlen,fposts,gposts,sorted_keys,irnkl2,n_cpu,my_rank)

            ! Check if sort finished and find inner doublets
            do i=2,npnew
                if (sorted_keys(i) .lt. sorted_keys(i-1)) then
                    write (*,'(a,i5,2z20)') 'Sorting locally failed i,sorted_keys(i),sorted_keys(i-1)=',i,sorted_keys(i),sorted_keys(i-1)
                    call MPI_ABORT(MPI_COMM_WORLD,ierr)
                    stop
                endif
                if (sorted_keys(i) == sorted_keys(i-1)) then
                    kick_part(i-1,1) = 1
                    kick_part(i-1,2) = i
                    dcount = dcount + 1
                endif
            end do

            ! swap end items
            call MPI_ISEND(sorted_keys(1), 1, MPI_INTEGER8, prev, 1, MPI_COMM_WORLD, handle(1), ierr)
            call MPI_REQUEST_FREE(handle(1),ierr)
            call MPI_RECV(tmp_next, 1, MPI_INTEGER8, next, 1, MPI_COMM_WORLD, status, ierr)
            call MPI_ISEND(sorted_keys(npnew), 1, MPI_INTEGER8, next, 1, MPI_COMM_WORLD, handle(1), ierr)
            call MPI_REQUEST_FREE(handle(1),ierr)
            call MPI_RECV(tmp_prev, 1, MPI_INTEGER8, prev, 1, MPI_COMM_WORLD, status, ierr)

             ! Check if sort globally finished
            if (my_rank .ne. n_cpu-1) then
                if (tmp_next .lt. sorted_keys(npnew)) then
                    write (*,*) 'sorted_keys(npnew), sorted_keys(1) from',my_rank+1, '=',sorted_keys(npnew),tmp_next,npnew
                    call MPI_ABORT(MPI_COMM_WORLD,ierr)
                    stop
                endif
            endif

            ! Check for boundary identical keys
            if (tmp_next == sorted_keys(npnew)) then
                dcount = dcount + 1
                write(*,*) my_rank,' detected cross-boundary doublet with',next
            endif
            if (tmp_prev == sorted_keys(1)) then
                dcount = dcount + 1
                write(*,*) my_rank,' detected cross-boundary doublet with',prev
            endif
            call MPI_ALLREDUCE(dcount,dcount_glo,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)

            ! Set up particle structure - keys and source_pe are dummies
            m_np = npnew
            ship_parts(1:npold) = m_part(indxl(1:npold))

            ! perform permute
            call MPI_ALLTOALLV(ship_parts, islen, fposts, mpi_type_particle, &
                               get_parts, irlen, gposts, mpi_type_particle, MPI_COMM_WORLD,ierr)
            m_part(irnkl(1:m_np)) = get_parts(1:m_np)
            m_part(1:m_np)%key = sorted_keys(1:m_np)
            m_part(1:m_np)%pid = my_rank

        end do ! END: Loop over sorting algorithm, as long as we find doublets (should stop after two steps!)

    end subroutine sort_remesh


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
