!!!!!!!!!!!!!!!!!!!!
!! zufallszahlen module
!!
!! Enthaelt Methoden zur Erzeugung von verschieden verteilten Zufallszahlen
!!!!!!!!!!!!!!!!!!!!

MODULE zufall
    use variables
    implicit none

    CONTAINS

!==================================================================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> initialize random number generator
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE init_rng()
        implicit none

        integer :: rsize,i
        integer, allocatable :: rseed(:)
        real*8 ::  ran

        rng=1

        if (rng==0) then
            call random_seed(size = rsize)
            allocate(rseed(rsize))
            rseed = my_rank + [(i*144,i=1,rsize)]
            call random_seed(put = rseed)
            deallocate(rseed)
        else if (rng==1) then
            ran=par_rand(my_rank+rngseed)
        end if
    end subroutine

!======================================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> fills list with gaussian random numbers
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE random_gauss_list(list,mu,sigma)
        implicit none
  
        real*8, intent(inout) :: list(:)
        real*8  :: v(2), pi, r, p, mu, sigma
        integer :: n, i
    
        pi = 2.0_8*acos(0.0_8)
        n  = size(list)

        DO i=1, n, 2
            v(1)=rnd_num()
            v(2)=rnd_num()
            !call random_number(v)
            r = sqrt(-2.0_8 * log(v(1)))
            p = 2.0_8*pi*v(2)
            list(i)                = r * sin(p)*sigma+mu
            if((i+1)<=n) list(i+1) = r * cos(p)*sigma+mu
        END DO
  
  END SUBROUTINE

!======================================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Samples random number x distributed according to
    !> f(v)=m/T * v * exp(-m/(2T)*(v-v0)**2), vtherm=sqrt(T/m)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE random_drifting_gaussian_flux(x,F_table,v_table)
        use helper
        implicit none

        real*8, intent(inout) :: x
        real*8, intent(in)  :: F_table(:),v_table(:)
        real*8  ::y,v_u,v_l,F_u,F_l
        integer :: u,l

        u=0
        l=0
        F_u=0._8
        F_l=0._8

        y=rnd_num()
        call find_yinx(y,F_table,l,u,F_l,F_u)

        v_u = v_table(u)
        v_l = v_table(l)
        x=v_l+(v_u-v_l)*(y-F_l)/(F_u-F_l)


    END SUBROUTINE
!======================================================================================

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Fills "list" with random numbers distributed according to
    !> f(v)=m/T * v * exp(-m/(2T)*(v-v0)**2), vtherm=sqrt(T/m)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE random_drifting_gaussian_flux_list(list,F_table,v_table)
        use helper
        implicit none

        real*8, intent(inout) :: list(:)
        real*8, intent(in)  :: F_table(:),v_table(:)
        real*8  ::y,v_u,v_l,F_u,F_l
        integer :: n,i,u,l

        n  = size(list)

        DO i=1, n
            y=rnd_num()
            call find_yinx(y,F_table,l,u,F_l,F_u)
            v_u = v_table(u)
            v_l = v_table(l)
            list(i)=v_l+(v_u-v_l)*(y-F_l)/(F_u-F_l)
        END DO

    END SUBROUTINE
!======================================================================================


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Fills "list" with random numbers distributed according to
    !> f(v)=m/T * v * exp(-m/(2T)*v**2), vtherm=sqrt(T/m)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       SUBROUTINE random_gaussian_flux_list(list,vtherm)
        implicit none
  
        real*8, intent(inout) :: list(:)
        real*8  :: vtherm,y
        integer :: n,i
    
        n  = size(list)

        DO i=1, n
            y=rnd_num()
            !call random_number(y)
            list(i)=vtherm*sqrt(-2.*log(1-y))          
        END DO
  
    END SUBROUTINE
!======================================================================================


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Generates random numbers(x) distributed according to
    !> f(v)=m/T * v * exp(-m/(2T)*v**2), vtherm=sqrt(T/m)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE random_gaussian_flux(x,vtherm)
        implicit none
  
        real*8, intent(inout) :: x
        real*8  :: vtherm,y

        y=rnd_num()
        !call random_number(y)
        x=vtherm*sqrt(-2.*log(1-y))          
          
    END SUBROUTINE


!======================================================================================
    subroutine get_uniformly_distributed_number(r, opt_min, opt_max)
        implicit none
        real(KIND=8), intent(out) :: r
        real(KIND=8), intent(in), optional :: opt_min, opt_max

        real(KIND=8) :: rmin, rmax

        rmin = 0.0_8
        rmax = 1.0_8

        if (present(opt_min)) rmin = opt_min
        if (present(opt_max)) rmax = opt_max

        r = rnd_num()
        r = r*(rmax-rmin) + rmin

    end subroutine get_uniformly_distributed_number


!======================================================================================
    subroutine get_uniformly_distributed_numbers(r, opt_min, opt_max)
        implicit none
        real(KIND=8), intent(inout) :: r(:)
        real(KIND=8), intent(in), optional :: opt_min, opt_max

        real(KIND=8) :: rmin, rmax
        integer :: ir

        rmin = 0.0_8
        rmax = 1.0_8

        if (present(opt_min)) rmin = opt_min
        if (present(opt_max)) rmax = opt_max

        do ir=1, size(r)
            r(ir) = rnd_num()
            r(ir) = r(ir)*(rmax-rmin) + rmin
        end do

    end subroutine get_uniformly_distributed_numbers


!======================================================================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> portable random number generator, see numerical recipes
    !> check for the random numbers:
    !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
    !> the parameter iseed is optional
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function par_rand(iseed)
        implicit none
        real :: par_rand
        integer, intent(in), optional :: iseed

        integer, parameter :: IM1  = 2147483563
        integer, parameter :: IM2  = 2147483399
        real,    parameter :: AM   = 1.0/IM1
        integer, parameter :: IMM1 = IM1-1
        integer, parameter :: IA1  = 40014
        integer, parameter :: IA2  = 40692
        integer, parameter :: IQ1  = 53668
        integer, parameter :: IQ2  = 52774
        integer, parameter :: IR1  = 12211
        integer, parameter :: IR2  = 3791
        integer, parameter :: NTAB = 32
        integer, parameter :: NDIV = 1+IMM1/NTAB
        real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
        real,    parameter :: RNMX = 1.0 - eps_

        integer :: j, k
        integer, volatile, save :: idum  = -1
        integer, volatile, save :: idum2 =  123456789
        integer, volatile, save :: iy    =  0
        integer, volatile, save :: iv(NTAB)


        if (idum <=0 .or. present(iseed)) then
            if (present(iseed)) then
                idum = iseed
            else
                if (-idum < 1) then
                    idum = 1
                else
                    idum = -idum
                endif
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
        par_rand = AM*iy
        if (par_rand > RNMX) par_rand = RNMX

    end function par_rand

    function rnd_num()
        implicit none

        real :: rnd_num

        if (rng==0) then
            call random_number(rnd_num)
        else if (rng==1) then
            rnd_num=par_rand()
        end if

    end function

END MODULE



!==================================================================================
!==================================================================================
!==================================================================================

MODULE poisson_disc_sampling
    use variables
    implicit none

    integer, parameter :: D = 3 !Dimension
    integer :: MaxCandidates
    real(KIND=8) :: Radius
    real(KIND=8) :: CellSize

    integer, allocatable :: Grid(:,:,:)
    integer :: GridShape(3)

    real(KIND=8) :: gDomainMin(D), gDomainMax(D)
    real(KIND=8) :: gDomainSize(D)

    real(KIND=8), allocatable :: ActiveSamples(:,:), AllSamples(:,:)
    integer :: iAllSamples, iActiveSamples, nSamples, nActiveSamples

    CONTAINS


!==================================================================================
    real(KIND=8) function distance(a, b)
        implicit none

        real(KIND=8), intent(in) :: a(:), b(:)
        integer :: j

        distance = 0.
        do j = 1, D
            distance = (a(j) - b(j))**2 + distance
        end do

        distance = sqrt(distance)
    end function distance

!==================================================================================
    subroutine init()
        use zufall
        implicit none

        real(KIND=8), allocatable :: ran(:)

        !This is an educated guess. The number of sampled positions will be close to, but smaller nSamples
        Radius =  0.855 * (product(gDomainSize)/nSamples)**(1/dble(D))

        call init_grid()

        allocate(ActiveSamples(D, nSamples))
        ActiveSamples = 0.0_8
        allocate(AllSamples(D, nSamples))
        AllSamples = 0.0_8

        nActiveSamples = 1
        iActiveSamples = 1
        iAllSamples = 1

        !sample first position
        allocate(ran(D))
        call get_uniformly_distributed_number(ran(1), gDomainMin(1), gDomainMax(1))
        call get_uniformly_distributed_number(ran(2), gDomainMin(2), gDomainMax(2))
        call get_uniformly_distributed_number(ran(3), gDomainMin(3), gDomainMax(3))
        ActiveSamples(:, iActiveSamples) = ran
        AllSamples(:, iAllSamples) = ActiveSamples(:, iActiveSamples)
        call add_sample_to_grid()

    end subroutine init


!==================================================================================
    subroutine init_grid()
        implicit none

        integer :: nx, ny, nz

        CellSize = Radius/sqrt(dble(D))

        nx = ceiling(gDomainSize(1) / CellSize)
        ny = ceiling(gDomainSize(2) / CellSize)
        nz = ceiling(gDomainSize(3) / CellSize)

        allocate(Grid(-1:nx+2,-1:ny+2,-1:nz+2))
        Grid = -1
        GridShape = (/nx,ny,nz/)

    end subroutine init_grid


!==================================================================================
    subroutine add_sample_to_grid()
        implicit none

        integer :: GridPosition(D)

        GridPosition(:) = int((AllSamples(:,iAllSamples) - gDomainMin(:)) / CellSize) + 1

        if (Grid(GridPosition(1),GridPosition(2),GridPosition(3)) > -1) then
            write(*,*) "Grid Position besetzt. Fehler im Algorithmus"
            stop
        else
            Grid(GridPosition(1),GridPosition(2),GridPosition(3)) = iAllSamples
        end if

        !Add Sample to ghost cells too, if they are close to the boundary
        Grid(-1:0, :, :) = Grid(GridShape(1)-1:GridShape(1), :, :)
        Grid(GridShape(1)+1:GridShape(1)+2, :, :) = Grid(1:2, :, :)
        Grid(:, -1:0, :) = Grid(:, GridShape(2)-1:GridShape(2), :)
        Grid(:, GridShape(2)+1:GridShape(2)+2, :) = Grid(:, 1:2, :)
        Grid(:, :, -1:0) = Grid(:, :, GridShape(3)-1:GridShape(3))
        Grid(:, :, GridShape(3)+1:GridShape(3)+2) = Grid(:, :, 1:2)

    end subroutine add_sample_to_grid


!==================================================================================
    subroutine sample_new_position()
        use zufall
        implicit none

        real(KIND=8) :: ActiveSample(D), Candidate(D)
        logical :: CandidateAccepted
        integer :: iCandidate
        real(KIND=8) :: dist

        CandidateAccepted = .False.
        ActiveSample = ActiveSamples(:,iActiveSamples)

        do iCandidate = 1, MaxCandidates
            do
                call get_uniformly_distributed_numbers(Candidate, -2*Radius, 2*Radius)
                Candidate = Candidate + ActiveSample
                dist = distance(Candidate, ActiveSample)
                if ( (Radius <= dist) .and. (2*Radius >= dist) ) exit
            end do
            call test_candidate(CandidateAccepted, Candidate)
            if (CandidateAccepted) then
                iAllSamples = iAllSamples + 1
                nActiveSamples = nActiveSamples + 1
                AllSamples(:, iAllSamples) = Candidate
                ActiveSamples(:, nActiveSamples) = Candidate
                call add_sample_to_grid()
                exit
            end if
        end do

        if (CandidateAccepted .eqv. .False.) then
            ActiveSamples(:, iActiveSamples) = ActiveSamples(:, nActiveSamples)
            nActiveSamples = nActiveSamples - 1
        end if

    end subroutine sample_new_position


!==================================================================================
    subroutine test_candidate(CandidateAccepted, Candidate)
        implicit none

        logical, intent(out) :: CandidateAccepted
        real(KIND=8), intent(in) :: Candidate(D)

        real(KIND=8) :: NeighborSamplePosition(D), rMirror(3)
        integer :: GridPosition(D), NeighborGridPosition(D)
        integer :: ix,iy,iz,j
        real(KIND=8) :: dist

        GridPosition(:) = int((Candidate(:) - gDomainMin(:)) / CellSize) + 1
        rMirror = 0.0_8

        if ( (any(Candidate-gDomainMin < 0.)) .or. (any(Candidate-gDomainMax > 0.)) ) then
            CandidateAccepted = .False.
            return
        end if

        !Grid cell is occupied
        if (Grid(GridPosition(1),GridPosition(2),GridPosition(3)) > -1) then
            CandidateAccepted = .False.
            return
        end if

        !have a look at cells in the direct vicinity
        do ix = -2,2
            do iy = -2,2
                do iz = -2,2
                    NeighborGridPosition = GridPosition + (/ix,iy,iz/)

                    !this is the cell of the new sample
                    if ( (ix==0) .and. (iy==0) .and. (iz==0) ) cycle

                    !Neighbor cell is empty
                    if (Grid(NeighborGridPosition(1),NeighborGridPosition(2),&
                             NeighborGridPosition(3)) < 0) then
                        cycle
                    end if

                    !Neighbors in Ghost cells
                    rMirror = 0.0_8
                    if ( (any(NeighborGridPosition < 1)) .or. (any(NeighborGridPosition - GridShape > 0)) ) then
                        do j=1,D
                            if (NeighborGridPosition(j) < 1) rMirror(j) = -gDomainSize(j)
                            if (NeighborGridPosition(j) > GridShape(j)) rMirror(j) = gDomainSize(j)
                        end do
                    end if

                    NeighborSamplePosition = AllSamples(:, Grid(NeighborGridPosition(1), &
                                                                NeighborGridPosition(2), &
                                                                NeighborGridPosition(3)))
                    NeighborSamplePosition = NeighborSamplePosition + rMirror

                    dist = distance(Candidate, NeighborSamplePosition)
                    if (dist < Radius) then
                        CandidateAccepted = .False.
                        return
                    end if
                end do
            end do
        end do

        CandidateAccepted = .True.
        return

    end subroutine test_candidate


!==================================================================================
    subroutine poisson_sampler(DomainMin, DomainMax, positions_x, positions_y, positions_z)
        use zufall
        implicit none

        real(KIND=8), intent(in) :: DomainMin(:), DomainMax(:)
        real(KIND=8), intent(out) :: positions_x(:),positions_y(:),positions_z(:)

        real(KIND=8) :: ran
        integer :: j=1

        MaxCandidates = 100

        nSamples = size(positions_x)
        gDomainMin(:)= DomainMin(:)
        gDomainMax(:) = DomainMax(:)
        gDomainSize(:) = gDomainMax(:) - gDomainMin(:)

        call init()

        do
            do while ((nActiveSamples > 0) .and. (iAllSamples < nSamples))
                call get_uniformly_distributed_number(ran, 1._8, dble(nActiveSamples))
                iActiveSamples = int(ran)
                call sample_new_position()
            end do
            if (iAllSamples == nSamples) then !this means that the algorithm did not finish
                                              !in a sense that the maximum number of samples was reached
                                              !before all Active Samples were processed (i.e. the domain was filled)
                Radius = Radius*1.001         !increased radius will decrease the number of particles when restarting the algorithm
                write(*,'(a)', advance='yes') "Poisson disc sampler reached maximum number of samples before filling the domain. "
                write(*,'(a,i8)') "Active Samples left: ",nActiveSamples
                write(*,'(a,es12.4)')        "Restarting with larger radius R = R*1.001 = ", Radius
            else
                exit
            end if
            ActiveSamples(:,1:) = 0.0_8
            AllSamples(:,1:) = 0.0_8
            Grid = -1

            nActiveSamples = 1
            iActiveSamples = 1
            iAllSamples = 1

            call add_sample_to_grid()
        end do

        positions_x(1:iAllSamples) = AllSamples(1,1:iAllSamples)
        positions_y(1:iAllSamples) = AllSamples(2,1:iAllSamples)
        positions_z(1:iAllSamples) = AllSamples(3,1:iAllSamples)
        do j = iAllSamples + 1, nSamples
             call get_uniformly_distributed_number(positions_x(j), gDomainMin(1), gDomainMax(1))
             call get_uniformly_distributed_number(positions_y(j), gDomainMin(2), gDomainMax(2))
             call get_uniformly_distributed_number(positions_z(j), gDomainMin(3), gDomainMax(3))
        end do

        if (allocated(Grid)) deallocate(Grid)
        if (allocated(ActiveSamples)) deallocate(ActiveSamples)
        if (allocated(AllSamples)) deallocate(AllSamples)

    end subroutine poisson_sampler
END module poisson_disc_sampling
