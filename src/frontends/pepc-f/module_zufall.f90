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

        integer :: rsize,i,rngseed
        integer, allocatable :: rseed(:)
        real*8 ::  ran

        rng=1
        rngseed=0

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
