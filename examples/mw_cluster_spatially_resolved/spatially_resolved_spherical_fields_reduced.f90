module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3


module module_data
  implicit none

  integer :: Na, Nt, Nt_max
  
  integer, parameter :: NUMCOMPONENTS = 1
  integer, parameter :: MYCOMPONENT = 4
  
  integer :: NR, NTheta, NPhi
  real*8 :: maxR, rion, relectron
  
  ! indices: spatial component (iR,iT,iP), time/frequency index (t/w)
  real*8, allocatable, dimension(:,:,:,:,:) :: observable
  real*8, allocatable, dimension(:) :: abscissa ! will hold time/frequency values
  real*8, allocatable, dimension(:,:,:,:) :: grid
  
  contains
  
    subroutine free_data()
      implicit none
      
      if (allocated(observable))    deallocate(observable)
      if (allocated(abscissa))      deallocate(abscissa)
      if (allocated(grid))          deallocate(grid)
      
    end subroutine free_data
    
    
    subroutine load_grid(filename)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      
      integer :: idx, iR, iTheta, iPhi
      real*8 :: r
      
      write(*,'("[STATUS] ", "load_grid")')

      open(87,file=trim(filename),status='old',position='rewind',action='read', FORM='unformatted')
      
      read(87) Nt_max, maxR, NR, NTheta, NPhi

      Na = (NR+1) * (NTheta+1) * (NPhi+1)

      write(*, '("# maxR     = ", g15.5)') maxR
      write(*, '("# NR       = ",   i15)') NR
      write(*, '("# NTheta   = ",   i15)') NTheta
      write(*, '("# NPhi     = ",   i15)') NPhi
      
      allocate(grid(0:NR, 0:NTheta, 0:NPhi, 1:3))
      
      idx = 0
      do iR = 0,NR
        do iTheta = 0,NTheta
          do iPhi = 0,NPhi
	    idx = idx + 1
            call progress(idx, Na)
	    read (87) grid(iR, iTheta, iPhi, 1:3)
          end do
        end do
      end do
 
      write(*,*) ! progress bar line break
      
      close(87)

    end subroutine
    
    
    logical function load_data(filename, starttime_fs)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      real*8, intent(in) :: starttime_fs
      integer :: itime
      real*8 :: time_fs
      integer :: ios, iR, iTheta, iPhi
      integer :: datindex, lineindex
      real*8, allocatable :: dummy(:,:,:,:)
      
      load_data = .false.

      write(*,'("[STATUS] ", "load_data")')

      open(87,file=trim(filename),status='old',position='rewind',action='read', FORM='unformatted',iostat=ios)
      
      if (ios .ne. 0) return
      
      read(87) Nt_max, maxR, rion, relectron, NR, NTheta, NPhi

      Na = (NR+1) * (NTheta+1) * (NPhi+1)

      write(*, '("# Nt_max    = ",   i15)') Nt_max
      write(*, '("# maxR      = ", g15.5)') maxR
      write(*, '("# rion      = ", g15.5)') rion
      write(*, '("# relectron = ", g15.5)') relectron
      write(*, '("# NR        = ",   i15)') NR
      write(*, '("# NTheta    = ",   i15)') NTheta
      write(*, '("# NPhi      = ",   i15)') NPhi
      write(*, '("# Na        = ",   i15)') Na
      
      allocate(observable(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi, Nt_max))
      allocate(dummy(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi))
      allocate(abscissa(1:Nt_max))
      
      datindex  = 1
      lineindex = 0
      
      do
        read(87, iostat=ios) itime, time_fs
	
	lineindex = lineindex + 1
	
	if (ios .ne. 0) exit
	
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
              read (87) dummy(1:4, iR, iTheta, iPhi)
            end do
          end do
        end do
	
	if (time_fs >= starttime_fs) then
	  datindex = datindex + 1
            observable(NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi, datindex) = dummy(MYCOMPONENT, 0:NR, 0:NTheta, 0:NPhi)
	    abscissa(datindex) = time_fs
	endif
	
	call progress(lineindex, Nt_max, time_fs >= starttime_fs)

      end do
      
      write(*,*) ! progress bar line break

      close(87)
      
      Nt = datindex ! store number of actually available timesteps
      write(*, '("# Nt       = ",   i15)') Nt

      load_data = (Nt > 2)

      if (.not. load_data) then
        write(*,'(a,/,a)') 'Did not find a sufficient number of timesteps.', 'Reduce starttime_fs or check your data'
      endif
      
    end function load_data
    
    subroutine observable_spectrum_to_file(filename_in)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename_in
      character(1024) :: filename

      integer :: i, iR, iTheta, iPhi, c
      
      character(8) :: endings(NUMCOMPONENTS) = ['_phi.dat']
      
      do c=1,NUMCOMPONENTS
        filename = trim(filename_in)//trim(endings(c))
 
	write(*,'("[STATUS] ", "observable_spectrum_to_file: ",a)') filename

	open(24,file=trim(filename),status='unknown',position='rewind',action='write')

	do i=1,Nt
          call progress(i, Nt)

          write(24,'(g15.5)', advance='no') abscissa(i)

          do iR = 0,NR
            do iTheta = 0,NTheta
              do iPhi = 0,NPhi
		write(24,'(x,g15.5)', advance='no') observable(c, iR, iTheta, iPhi, i) ! TODO: for now we only dump the potential
              end do	
	    end do
          end do

	  write(24,*) ! line break in output file
	end do

	close(24)

	write(*,*) ! progress bar line break
      end do
      
    end subroutine
    
  
end module


module computations
  use module_data
  use progress_bar
  implicit none
  
  contains
  
  subroutine observable_fourier_forward()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(Nt) :: tmparray
    integer :: c, iR, iTheta, iPhi
    integer :: nstep, nstep_max
    
    write(*,'("[STATUS] ", "observable_fourier_forward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)

    nstep = 0
    nstep_max = (NR+1) * (NTheta+1) * (NPhi+1)
    
    do iR = 0,NR
      do iTheta = 0,NTheta
        do iPhi = 0,NPhi
	  nstep = nstep + 1
	  call progress(nstep, nstep_max)
	  
          do c=1,NUMCOMPONENTS
            tmparray(1:Nt) = observable(c,iR, iTheta, iPhi, 1:Nt)
     
            ! perform FFT
            call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
            observable(c,iR, iTheta, iPhi, 1:Nt) = abs(tmparray(1:Nt))
	  end do
	end do
      end do
    end do
    
    call fftw_destroy_plan(fftw_plan)
    
    write(*,*) ! progress bar line break

  end subroutine


  subroutine abscissa_convert_to_frequencies()
    implicit none
    real*8 :: delta_t, wmax, delta_w
    integer :: i
    real*8, parameter :: pi = 3.14159265358979323846_8
    
    write(*,'("[STATUS] ", "abscissa_convert_to_frequencies")')

    delta_t = abscissa(2)-abscissa(1)
    wmax = 2._8*pi*1/delta_t
    delta_w = wmax/Nt
    
    abscissa(1:Nt)        = [((i-1)*delta_w,i=1,Nt)]
    abscissa(Nt+1:Nt_max) = 0._8
    
  end subroutine

  
    
  subroutine analysis_workflow(dirname_in, starttime_fs, write_vtk)
    use spherical_fourier
    implicit none
    character(*), intent(in) :: dirname_in
    real*8, intent(in) :: starttime_fs
    logical, intent(in) :: write_vtk
    
    if (load_data(trim(dirname_in)//'/field_spherical.dat', starttime_fs)) then
      call load_grid(trim(dirname_in)//'/field_spherical.dat_grid.dat')
      
      call observable_fourier_forward()
      call abscissa_convert_to_frequencies()

      call observable_spectrum_to_file(trim(dirname_in)//'/field_spherical_spectrum')
      
    endif

    call free_data()

  end subroutine
  
end module


program spatially_resolved_spherical_fields
  use module_data
  use computations
  implicit none
  include 'mpif.h'
  integer :: argc, i
  character*256 :: dirname_in
  integer rank, size, ierror
  
  argc = command_argument_count()
  
  if (argc < 1) then
    write(*,*) 'Call with directories to be processed as arguments. Exiting.'
    stop
  endif

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  
  do i=1,argc
    if (modulo(argc, size) == rank) then
      write(*,*) 'Rank ', rank, ' processing arg', i
      call get_command_argument(i, dirname_in)
    
      write(*,'(2/"Reading data from directory ",a)') trim(dirname_in)

      call analysis_workflow(trim(dirname_in), 100.0_8, .false.)
    endif
  end do

  call MPI_FINALIZE(ierror)  
   

end program
