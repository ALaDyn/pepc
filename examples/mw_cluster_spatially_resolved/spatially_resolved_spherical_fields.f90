module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3

module progress_bar
  implicit none
  
  logical, private, parameter :: showprogress            = .true.
  integer, private, parameter :: barlen                  = 100
  character, private, parameter :: progress_char         = '+'
  character, private, parameter :: progress_char_special = '*'

  contains
  
  subroutine progress(j,k,special)  
    implicit none  
    integer, intent(in):: j,k  
    logical, intent(in), optional :: special
    character(len=18) :: barfmt = "(a1,I3,a,x,???a,$)"
    
    integer :: value
    
    character :: bar(0:barlen+1)
    character :: pchar
    
    if (.not. showprogress) return
    
    pchar = progress_char
    value = (barlen*j)/k
    
    if (present(special)) then
      if (special) then
        pchar = progress_char_special
      endif
    endif
    
    bar(0)        = "|"
    bar(barlen+1) = "|"
    
    bar(1:value)        = pchar
    bar(value+1:barlen) = " "

    write(barfmt(12:14), '(I3.3)') barlen+2
    
    ! print the progress bar.  
    write(6, barfmt) char(13), (j*100)/k, "%", bar
    
  end subroutine progress  

end module



module module_data
  implicit none

  integer :: Na, Nt, Nt_max
  
  integer :: NR, NTheta, NPhi
  real*8 :: maxR
  
  ! indices: spatial component (iR,iT,iP), time/frequency index (t/w)
  real*8, allocatable, dimension(:,:,:,:,:) :: observable
  real*8, allocatable, dimension(:) :: abscissa ! will hold time/frequency values
  
  contains
  
    subroutine free_data()
      implicit none
      
      if (allocated(observable))    deallocate(observable)
      if (allocated(abscissa))      deallocate(abscissa)
      
    end subroutine free_data
    
    
    logical function load_data(filename, starttime_fs)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename
      real*8, intent(in) :: starttime_fs
      integer :: itime
      real*8 :: time_fs
      integer :: ios, iR, iTheta, iPhi
      integer :: datindex, lineindex
      real*8 :: dummy(1:4)

      open(87,file=trim(filename),status='old',position='rewind',action='read', FORM='unformatted')
      
      read(87) Nt_max, maxR, NR, NTheta, NPhi

      Na = (NR+1) * (NTheta+1) * (NPhi+1)

      write(*, '("# Nt_max   = ",   i15)') Nt_max
      write(*, '("# maxR     = ", g15.5)') maxR
      write(*, '("# NR       = ",   i15)') NR
      write(*, '("# NTheta   = ",   i15)') NTheta
      write(*, '("# NPhi     = ",   i15)') NPhi
      write(*, '("# Na       = ",   i15)') Na
      
      allocate(observable(1:4, 0:NR, 0:NTheta, 0:NPhi, Nt_max))
      
      datindex  = 0
      lineindex = 0
      
      do
        read(87, iostat=ios) itime, time_fs
	
	if (ios .ne. 0) exit
	
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
              read (87) dummy(1:4)
            end do
          end do
        end do
	
	if (time_fs >= starttime_fs) then
	  datindex = datindex + 1
            observable(1:4, iR, iTheta, iPhi, datindex) = dummy(1:4)
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
    
    
    subroutine observable_write_to_vtk(fileprefix)
      use progress_bar
      implicit none
      character(*), intent(in) :: fileprefix
      
      character*256 :: filename
      integer :: i
      integer :: iR, iTheta, iPhi
      real*8 :: theta, phi, r
      real*8, parameter :: pi = 3.14159265358979323846_8
      integer, parameter :: VTK_HEXAHEDRON = 12
      real*8, dimension(0:NR,0:Ntheta,0:NPhi,1:3) :: points
      
      write(*,'("[STATUS] ", "observable_write_to_vtk: ",a)') fileprefix
      
      ! prepare array with grid coordinates
      do iR = 0,NR
        do iTheta = 0,NTheta
          do iPhi = 0,NPhi
            theta =      pi /  Ntheta * iTheta
	    phi   = 2._8*pi /  Nphi   * iPhi
	    r     =    maxR /  Nr     * ir
	      
	    points(iR, iTheta, iPhi, 1:3) = r * [ cos(phi)*sin(theta) , sin(phi)*sin(theta), cos(theta) ]
          end do	
	end do
      end do
      

      do i = 1,Nt
        call progress(i, Nt)

        write(filename,'(a,"_",I7.7,".vtk")') fileprefix, i
	
        open(24,file=trim(filename),status='unknown',position='rewind',action='write')
	
	! write VTK header
	write(24,'("# vtk DataFile Version 2.0")')
	write(24,'("spherical fields")')
	write(24,'("ASCII")')
	! Grid definition
	write(24,'("DATASET UNSTRUCTURED_GRID")')
	! grid point coordinates
	write(24,'("POINTS ", I0, " double")') (NR+1)*(NTheta+1)*(NPhi+1)
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(3(x,g25.12))') points(iR, iTheta, iPhi, 1:3)
            end do	
	  end do
        end do
	  write(24,*)
	! grid cell types
	write(24,'("CELL_TYPES ", I0)') NR*NTheta*NPhi
        do iR = 0,NR-1
          do iTheta = 0,NTheta-1
            do iPhi = 0,NPhi-1
	      write(24,'(x,I0)', advance='no') VTK_HEXAHEDRON
            end do	
	  end do
        end do
	  write(24,*)
	! grid cell definitions
	write(24,'("CELLS ", I0,x,I0)') NR*NTheta*NPhi, (8+1)*NR*NTheta*NPhi
        do iR = 0,NR-1
          do iTheta = 0,NTheta-1
            do iPhi = 0,NPhi-1
	      write(24,'(9(x,I0))') 8, &
	                            idx(iR  , mod(iPhi+1, NPhi), iTheta+1) , &
	                            idx(iR  ,     iPhi         , iTheta+1) , &
	                            idx(iR+1,     iPhi         , iTheta+1) , &
	                            idx(iR+1, mod(iPhi+1, NPhi), iTheta+1) , &
	                            idx(iR  , mod(iPhi+1, NPhi), iTheta  ) , &
	                            idx(iR  ,     iPhi         , iTheta  ) , &
	                            idx(iR+1,     iPhi         , iTheta  ) , &
	                            idx(iR+1, mod(iPhi+1, NPhi), iTheta  )
            end do	
	  end do
        end do
	  write(24,*)
	! point data
	write(24,'("POINT_DATA ", I0)') (NR+1)*(NTheta+1)*(NPhi+1)
	! electrical field
	write(24,'("SCALARS e_field double 3")')
	write(24,'("LOOKUP_TABLE default")')

        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(3(x,g25.12))') observable(1:3, iR, iTheta, iPhi, i)
            end do	
	  end do
        end do
	
	write(24,*)

	! electrical potential
	write(24,'("SCALARS potential double 1")')
	write(24,'("LOOKUP_TABLE default")')

        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      write(24,'(g25.12)') observable(4, iR, iTheta, iPhi, i)
            end do	
	  end do
        end do
	
        close(24)
      end do

      write(*,*) ! progress bar line break
      
    contains
    
      integer function idx(iR, iPhi, iTheta)
        implicit none
	integer, intent(in) :: iR, iTheta, iPhi
	idx = iR*(NTheta+1)*(NPhi+1) + iTheta * (NPhi+1) + iPhi
      end function
      
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
	  
          do c=1,4
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
    
    delta_t = abscissa(2)-abscissa(1)
    wmax = 2*pi*1/delta_t
    delta_w = wmax/Nt
    
    abscissa(1:Nt)        = [((i-1)*delta_w,i=1,Nt)]
    abscissa(Nt+1:Nt_max) = 0.
    
  end subroutine

  
    
  subroutine analysis_workflow(filename_in, starttime_fs, filename_ev)
    implicit none
    character(*), intent(in) :: filename_in
    character(*), intent(in) :: filename_ev ! filename to output the frequency-dependent eigenvalues to
    real*8, intent(in) :: starttime_fs
    
    if (load_data(filename_in, starttime_fs)) then
      call observable_write_to_vtk('fields/spherical_fields')
      call observable_fourier_forward()
      call observable_write_to_vtk('fields/spherical_fields_fft')
    endif

    call free_data()

  end subroutine
  
end module


program spatially_resolved_ccf
  use module_data
  use computations
  implicit none
  integer :: argc, i
  character*256, filename_in, filename_out
  
  argc = command_argument_count()
  
  call system('mkdir -p fields')

  if (argc < 1) then
    write(*,*) 'Call with file to be processed as first argument. Exiting.'
    stop
  endif

  do i=1,argc
    call get_command_argument(i, filename_in)
    filename_out = trim(filename_in(1:len(trim(filename_in))-4))//'_momentum_eigenvalues.dat'
    
    write(*,'(2/"Reading data from file ",a,/,"Writing data to file   ", a)') trim(filename_in), trim(filename_out)

    call analysis_workflow(trim(filename_in), 100.0_8, trim(filename_out))
  end do

  
  
 

end program
