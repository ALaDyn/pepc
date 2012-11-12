module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3



module data
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: RAWDATA_PX      = 1
  integer, parameter :: RAWDATA_PY      = 2
  integer, parameter :: RAWDATA_PZ      = 3
  integer, parameter :: RAWDATA_P       = 4
  
  integer, parameter :: NUM_RAWDATA = RAWDATA_P


  integer :: Nc, Nt, Nt_max
  
  integer :: NR, NTheta, NPhi
  real*8 :: maxR, rion, relectron

  ! indices: spatial component (c), time/frequency index (t/w)
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: observable
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: observableACF

  real*8, allocatable, dimension(:) :: abscissa ! will hold time/frequency values
  
  contains
  
    subroutine count_num_lines(filename, numlines)
      implicit none
      integer, intent(out) :: numlines
      character(*), intent(in) :: filename
      character(4096) :: dummy
      integer :: ios
      
      numlines = 0
      
      open(87,file=trim(filename),status='old',position='rewind',action='read',iostat=ios)

      if (ios .ne. 0) return

      do
        numlines = numlines + 1
	
        read(87,'(a)',iostat=ios) dummy
	
	if (ios .ne. 0) exit
	
      end do
      
      close(87)
    end subroutine
      
  
    subroutine allocate_data()
      implicit none
      
      allocate(observableACF(Nc+1, 2*Nt))
      
    end subroutine allocate_data
    
    
    subroutine free_data()
      implicit none
      
      if (allocated(observable))    deallocate(observable)
      if (allocated(observableACF)) deallocate(observableACF)
      if (allocated(abscissa))      deallocate(abscissa)
      
    end subroutine free_data
    
    
    logical function load_data(filename, fields, starttime_fs)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: fields(:)
      real*8, intent(in) :: starttime_fs
      integer :: iR, iTheta, iPhi, idata
      integer :: itime, Nentries
      real*8 :: time_fs
      integer :: ios
      integer :: datindex, a, c
      real*8 :: tmp(4)
      
      load_data = .false.
      
      call count_num_lines(filename, Nt_max)

      open(87,file=trim(filename),status='old',position='rewind',action='read',iostat=ios)
      
      if (ios .ne. 0) return

      Nc = size(fields)
      
      allocate(observable(Nc+1, 2*Nt_max))
      allocate(abscissa(Nt_max))

      observable = 0. ! this automatically ensures padding, i.e. the second half of the array is zero
      
      datindex = 0
      
      call progress(0, Nt_max-1)
      
      do

        read(87,'(i10,5g25.12,i12)',iostat=ios) itime, time_fs, tmp, Nentries
	
	if (ios .ne. 0) exit

	if (time_fs >= starttime_fs) then
          call progress(itime, Nt_max-1, .true.)
	  datindex = datindex + 1
	  
	  abscissa(datindex) = time_fs
          observable(1:Nc, datindex) = tmp(fields)
       	else
          call progress(itime, Nt_max, .false.)
	endif

      end do
      

      close(87)
      
      write(*,*) ! progress bar linebreak
      
      Nt = datindex ! store number of actually available timesteps
      write(*, '("# Nt       = ",   i15)') Nt

      load_data = (Nt > 2)

      if (.not. load_data) then
        write(*,'(a,/,a)') 'Did not find a sufficient number of timesteps.', 'Reduce starttime_fs or check your data'
      endif
      
    end function
  
    subroutine output_acf(filename)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      integer :: t
      character*24 :: formatstring
      
      write(formatstring,'("(g25.12,",I5.5,"(x,g25.12))")') Nc
      
      call progress(0,Nt)

      open(24,file=trim(filename),status='unknown',position='rewind',action='write')

      do t=1,Nt
        call progress(t,Nt)
        write(24, formatstring) abscissa(t), real(observableACF(1:Nc,t))
      end do  
      
      close(24)
      write(*,*) ! progress bar linebreak
      
    end subroutine


end module

module computations
  use data
  implicit none
  
  contains
  
  subroutine observable_fourier_forward()
    use progress_bar
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: c
    
    write(*,'("[STATUS] ", "observable_fourier_forward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)

    call progress(0, Nc)

    do c=1,Nc
        call progress(c, Nc)
      
        tmparray(1:2*Nt) = observable(c,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observable(c,1:2*Nt) = tmparray(1:2*Nt)
    end do
    
    write(*,*) ! progress bar line break
    
    call fftw_destroy_plan(fftw_plan)

  end subroutine


  subroutine abscissa_convert_to_frequencies()
    implicit none
    real*8 :: delta_t, wmax, delta_w
    integer :: i
    real*8, parameter :: pi = 3.14159265358979323846_8
    
    delta_t = abscissa(2)-abscissa(1)
    wmax = 2*pi*1/delta_t
    delta_w = wmax/Nt
    
    abscissa(1:Nt)        = [((i-1)*delta_w,i=1,Nt)]
    abscissa(Nt+1:Nt_max) = 0.
    
  end subroutine

  
  subroutine observable_spatial_average()
    use progress_bar
    implicit none
    integer :: a, t
    
    write(*,'("[STATUS] ", "observable_spatial_average")')

    call progress(0,2*Nt)

    do t=1,2*Nt
      call progress(t,2*Nt)
      observable(Nc+1,t) = 1./Nc * sum(observable(1:Nc, t))
    end do
    
    Nc = Nc + 1
    
    write(*,*) ! progress bar line break

  end subroutine


  subroutine observable_auto_correlate()
    use progress_bar
    implicit none
    integer :: c
    
    write(*,'("[STATUS] ", "observable_auto_correlate")')

    call progress(0,Nc)

    do c=1,Nc
      call progress(c,Nc)

      observableACF(c, 1:2*Nt) = observable(c,1:2*Nt) * conjg(observable(c,1:2*Nt))
    end do
    
    write(*,*) ! progress bar line break

  end subroutine
  
  
  subroutine observableACF_fourier_backward()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: c
    
    write(*,'("[STATUS] ", "observableACF_fourier_backward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_FORWARD, FFTW_ESTIMATE)
    
    do c=1,Nc
        tmparray(1:2*Nt) = observableACF(c,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observableACF(c,1:2*Nt) = tmparray(1:2*Nt) / (2._8*Nt)
    end do
    
    call fftw_destroy_plan(fftw_plan)
   
 end subroutine   


  subroutine observableACF_laplace_transform()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: c
    
    write(*,'("[STATUS] ", "observableACF_laplace_transform")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    do c=1,Nc
        tmparray(1:2*Nt) = observableACF(c,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observableACF(c,1:2*Nt) = tmparray(1:2*Nt)
    end do
    
    call fftw_destroy_plan(fftw_plan)
    
  end subroutine   
  
  
  subroutine analysis_workflow(filename_in, fields, starttime_fs, filename_acf)
    implicit none
    character(*), intent(in) :: filename_in
    character(*), intent(in) :: filename_acf
    integer, intent(in) :: fields(:)
    real*8, intent(in) :: starttime_fs
    
    if (load_data(filename_in, fields, starttime_fs)) then
      ! call observable_undo_volume_scaling() ! see (18) in [PRE84,036406], not necessary since we already output the observable without volume scaling, compare write_spatially_resolved_data () in module_diagnostics.f90
      call allocate_data()
      call observable_fourier_forward()
      call abscissa_convert_to_frequencies()
      call observable_spatial_average()
      call observable_auto_correlate()
      ! die folgenden beiden Schritte sind zueinander invers, sollten also obsolet sein - machen auch nach Test keinen Unterschied am Ergebnis
      !  call observableACF_fourier_backward()   ! jetzt steht D_c(t) in observableACF(c,t)
      !  call observableACF_laplace_transform()  ! jetzt steht D_c(w) in observableACF(c,w)
      call output_acf(filename_acf)
    endif

    call free_data()

  end subroutine
  
end module


program momentum_electron_acf
  use data
  use computations
  implicit none
  integer :: argc, i
  character*256, filename_in, filename_out, dirname
  
  argc = command_argument_count()

  if (argc < 1) then
    write(*,*) 'Call with directories to be processed as arguments. Exiting.'
    stop
  endif
  
  do i=1,argc
    call get_command_argument(i, dirname)
    filename_in = trim(dirname) // '/momentum_electrons.dat'
    filename_out = trim(dirname) // '/momentum_electrons_acf.dat'
    
    write(*,'(2/"Reading data from file ",a,/,"Writing data to file   ", a)') trim(filename_in), trim(filename_out)

    call analysis_workflow(trim(filename_in), [RAWDATA_PX, RAWDATA_PY, RAWDATA_PZ], 100.0_8, trim(filename_out))
  end do

  
  
 

end program
