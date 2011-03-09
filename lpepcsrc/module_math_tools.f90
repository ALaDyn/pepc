!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Several math tools, primarily matrix manipulations,
!> Legendre polynomials and related stuff. Additionally
!> functions to calculate the biggest power in given interval
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
      public sort_abs
      public bpi
      public bpi_bits

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Sorts the given values with a heap sort approach
        !> in order of ther absolute value
        !> compare (Numerical Recipes f90, p1171)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine sort_abs(arr)
            implicit none
            real*8, intent(inout) :: arr(:)
            integer :: i,n

            n = size(arr)

            do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
               call sift_down(i,n)
            end do

            do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
               ! during retirement/promotion (heap selection) phase.
               call swap( arr(1),arr(i) )      ! Clear space at end of array and retire top of heap into it
               call sift_down( 1,i-1)
            end do

          contains
            subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain
              integer, intent(in) :: l,r     ! the heap structure
              integer :: j,jold    ! index
              real*8 :: a

              a = arr(l)

              jold = l
              j = l + l
              do                   ! do while j <= r
                 if (j > r) exit
                 if (j < r) then
                   if (abs(arr(j)) < abs(arr(j+1))) j = j+1
                 endif
                 if (abs(a) >= abs(arr(j))) exit       ! Found a`s level, so terminate sift-down
                 arr(jold) = arr(j)
                 jold = j                    ! Demote a and continue
                 j = j+j
              end do
              arr(jold) = a                  ! Put a into its slot

            end subroutine sift_down

            subroutine swap(p,q)
                real*8 :: p,q, dum
                dum = p
                p=q
                q = dum
            end subroutine swap

          end subroutine sort_abs


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the biggest power to \a base in a given interval with the limits \a a and \a b.
        !>
        !> @param[in] a First Limit of interval.
        !> @param[in] b Second Limit of interval.
        !> @param[in] base Base for which the power shoul be found.
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer*8 function bpi(a, b, base)

          implicit none
          include 'mpif.h'

          integer :: ierr

          integer*8,intent(inout) :: a, b
          integer*8,intent(in) :: base
          integer*8 :: k
          integer*8 :: i 
          integer*8 :: powr
          integer*8 :: res
          
          powr = floor(log(REAL(b))/log(REAL(base)))
          do i = powr, 0, -1
             
             k = b/base**i
             res = k * base**i
             if (a.lt.res) then
                bpi=res ! return
                return
             end if
          end do

          bpi=-1
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          
        end function bpi
        

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the biggest power to \a base in a given interval with the limits \a a and \a b. In this 
        !> routine it is done by bit operations. 
        !>
        !> @param[in] a First Limit of interval.
        !> @param[in] b Second Limit of interval.
        !> @param[in] base Base for which the power should be found.
        !> @param[in] levels The number of levels in the tree.\see{treevars::nlev}
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer*8 function bpi_bits(a, b, base, levels)

          implicit none
          include 'mpif.h'

          integer :: ierr

          integer*8,intent(inout) :: a, b
          integer*8,intent(in) :: base
          integer,intent(in) :: levels
          integer*8 :: i
          integer*8 :: bn, pos

          do i=1,levels
             pos=3*(levels-i)
             if(ibits(a,pos,3).ne.ibits(b,pos,3))then
                bn=base**(levels-i)
                bpi_bits=b/bn*bn  ! return, Note: is an integer division
                return
             end if
          end do
          
          bpi_bits=-1
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

        end function bpi_bits
        
 end module module_math_tools
      




