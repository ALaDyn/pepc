!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_acf
      implicit none


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      type acf
        private
          integer :: Ntau
          integer :: tau
          real*8, allocatable :: Kt(:)
          integer, allocatable :: ctr(:)
          real*8,  allocatable :: oldvals(:,:)

        contains
          procedure :: initialize => acf_initialize
          procedure :: finalize => acf_finalize
          procedure :: addval => acf_addval
          procedure :: to_file => acf_to_file

      end type acf


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

      subroutine acf_to_file(acf_, filename)
        implicit none
        class(acf) :: acf_
        character(*) :: filename

        open(47,file=filename)
        write(47,'(g18.8)') acf_%Kt / acf_%ctr
        close(47)

      end subroutine

      subroutine acf_initialize(acf_, Nt_)
        implicit none
        class(acf) :: acf_
        integer, intent(in) :: Nt_

        acf_%Ntau = Nt_
        acf_%tau = 0

        allocate(acf_%Kt(0:acf_%Ntau))
        acf_%Kt = 0.
        allocate(acf_%ctr(0:acf_%Ntau))
        acf_%ctr = 0
        allocate(acf_%oldvals(1:3,1:acf_%Ntau))
      end subroutine


      subroutine acf_finalize(acf_)
        implicit none
        class(acf) :: acf_

        call acf_%to_file("momentum_Kt.dat")

        if (allocated(acf_%Kt))      deallocate(acf_%Kt)
        if (allocated(acf_%ctr))     deallocate(acf_%ctr)
        if (allocated(acf_%oldvals)) deallocate(acf_%oldvals)
      end subroutine



      subroutine acf_addval(acf_, val)
        implicit none
        class(acf) :: acf_
        real*8, intent(in) :: val(3)
        integer :: s

        acf_%tau = acf_%tau + 1

        acf_%oldvals(1:3,acf_%tau) = val

        do s = 0,acf_%tau-1
          acf_%Kt(s) = acf_%Kt(s) + dot_product(val, acf_%oldvals(1:3,acf_%tau - s))
          acf_%ctr(s) = acf_%ctr(s) + 1
        end do

      end subroutine


end module module_acf
