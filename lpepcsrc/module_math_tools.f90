!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Several math tools, primarily matrix manipulations,
!> Legendre polynomials and related stuff
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_math_tools
      use module_debug
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public LegendreP
      public mult_by_fac
      public div_by_fac
      public factorial
      public inverse3

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
        !> Computes the associated Legendre polynomial \f$P_l_m (x)\f$.
        !> Here m and l are integers satisfying  \f$0 \leq m \leq l\f$,
        !> while x lies in the range \f$-1 \leq x \leq 1\f$.
        !>
        !> Code fragment for \f$P_l^m(x)\f$ taken from
        !>
        !> Numerical Recipes in Fortran 77: The Art of Scientific Computing
        !>              (ISBN 0-521-43064-X)
        !> pg. 246ff
        !>
        !> and modified to give \f$P_l_m (x)\f$:
        !> \f$ P_l_m (x) = (-1)^m P_l^m (x) \f$, see
        !>
        !> Abramowitz and Stegun: Handbook of Mathematical Functions
        !> Section 8. Legendre Functions (pg. 332)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 function LegendreP(l,m,x)
          implicit none
          integer, intent(in) :: l, m
          real*8 ::x

          integer :: i,ll
          real*8 :: fact,pll,pmm,pmmp1,somx2

          if ( (m < 0) .or. (m > l) .or. (abs(x) > 1) ) then
            write(*,*) 'Invalid arguments for LegendreP'
            stop
          endif

          pmm = 1.0     ! Compute P_m^m

          if (m > 0) then
            somx2 = sqrt((1.-x)*(1.+x))
            fact  = 1.0

            do i = 1,m
               pmm  = -pmm * fact * somx2
               fact = fact+2.
            enddo
          endif

          if (l == m) then
            LegendreP = pmm
          else
            pmmp1 = x*(2*m+1)*pmm  ! Compute P_m+1^m

            if (l == m+1) then
              LegendreP = pmmp1
            else                  ! Compute P_l^m , l > m + 1
              do ll = m+2,l
                pll   = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                pmm   = pmmp1
                pmmp1 = pll
              enddo

              LegendreP = pll
            endif

            LegendreP = (-1)**m * LegendreP
          endif

          ! DEBUG (for comparison with Mathematica worksheet)
          !write(*, '("LegendreP[", I2.2, ", ", I2.2, ", ", F10.5, "],  ", D20.10)') l, m, x, LegendreP

        end function LegendreP


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Determinant of a real*8 2x2 matrix
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 function det2(mat)
          implicit none
          real*8, intent(in) :: mat(1:2,1:2)

          det2 = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

        end function det2


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Determinant of a real*8 2x2 matrix that is saved as a vector (columns first)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 function det2f(mat)
          implicit none
          real*8, intent(in) :: mat(1:4)

          det2f = mat(1)*mat(4)-mat(2)*mat(3)

        end function det2f

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Determinant of a real*8 3x3 matrix
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 function det3(mat)
          implicit none
          real*8, intent(in) :: mat(1:3,1:3)

          det3 = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + mat(1,3)*mat(2,1)*mat(3,2) &
               - mat(1,3)*mat(2,2)*mat(3,1) - mat(1,2)*mat(2,1)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2)

        end function det3


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real*8 function cofact(mat, i, j)
          implicit none
          real*8, intent(in) :: mat(1:3,1:3)
          integer, intent(in) :: i, j

          real*8 :: cf(1:2,1:2)

          cf(1:i-1,1:j-1) = mat(1  :i-1,1  :j-1)
          cf(i:   ,1:j-1) = mat(i+1:   ,1  :j-1)
          cf(1:i-1,j:   ) = mat(1  :i-1,j+1:   )
          cf(i:   ,j:   ) = mat(i+1:   ,j+1:   )

          cofact = det2(cf)

          write(*,*) 'cofact(A, i, j) = ', i, j, cofact

        end function cofact


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>  Simple 3x3-matrix inversion using Cramers rule
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function inverse3(m)
          implicit none
          real*8, intent(in) :: m(1:3,1:3)

          real*8 :: inverse3(1:3,1:3)
          real*8 :: test(1:3,1:3)
          real*8 :: det
          integer :: i,j
          logical :: failed

          det = det3(m)
          inverse3 = 0

          if (det == 0) then
            write (*,*) 'You are trying to invert a singular matrix - this is really evil'
            stop
          end if

          inverse3 = reshape([ & ! first column
                               det2f([m(2,2),m(3,2),m(2,3),m(3,3)]), &
                               det2f([m(2,3),m(3,3),m(2,1),m(3,1)]), &
                               det2f([m(2,1),m(3,1),m(2,2),m(3,2)]), &
                               ! 2nd column
                               det2f([m(1,3),m(3,3),m(1,2),m(3,2)]), &
                               det2f([m(1,1),m(3,1),m(1,3),m(3,3)]), &
                               det2f([m(1,2),m(3,2),m(1,1),m(3,1)]), &
                               ! 3rd column
                               det2f([m(1,2),m(2,2),m(1,3),m(2,3)]), &
                               det2f([m(1,3),m(2,3),m(1,1),m(2,1)]), &
                               det2f([m(1,1),m(2,1),m(1,2),m(2,2)]) ], [3,3])

          inverse3 = inverse3 / det

          ! just for testing purposes
          test = matmul(inverse3,m)
          failed = .false.

          do i=1,3
            do j=1,3
              if (i==j) then
                failed = failed .or. (abs(test(i,j)-1.D0) > 1.D-15)
              else
                failed = failed .or. (abs(test(i,j))      > 1.D-15)
              end if
            end do
          end do

          if (failed) then
            write (*,*) 'inverse3 failed for the following input matrix:', m, 'output: ', inverse3, 'test: ', test
            stop
          end if

        end function inverse3




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the factorial of the argument
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        recursive integer*8 function factorial(n) result(fact)
            implicit none
            integer, intent(in) :: n

            if (n<0) then
              write(*,*) "tried to calculate factorial of negative argument - you are evil"
              stop
            end if

            select case (n)
              case ( 0)
                fact =             1
              case ( 1)
                fact =             1
              case ( 2)
                fact =             2
              case ( 3)
                fact =             6
              case ( 4)
                fact =            24
              case ( 5)
                fact =           120
              case ( 6)
                fact =           720
              case ( 7)
                fact =          5040
              case ( 8)
                fact =         40320
              case ( 9)
                fact =        362880
              case (10)
                fact =       3628800
              case (11)
                fact =      39916800
              case (12)
                fact =     479001600
              case default
                fact = n * factorial(n-1)
            end select
        end function factorial


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Divides the argument x by n!
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function div_by_fac(x, n)
            implicit none
            complex*16, intent(in) :: x
            integer, intent(in) :: n

            integer :: i

            if (n <= 12) then
              div_by_fac = x / factorial(n)
            else
              div_by_fac = x / factorial(12)

              do i=13,n
                div_by_fac = div_by_fac / i
              end do
            end if
        end function div_by_fac


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Multiplies the argument x by n!
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function mult_by_fac(x, n)
            implicit none
            complex*16, intent(in) :: x
            integer, intent(in) :: n

            integer :: i

            if (n <= 12) then
              mult_by_fac = x * factorial(n)
            else
              mult_by_fac = x * factorial(12)

              do i=13,n
                mult_by_fac = mult_by_fac * i
              end do
            end if
        end function mult_by_fac


end module module_math_tools





