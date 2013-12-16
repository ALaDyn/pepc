module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3


module module_data
  implicit none

  integer :: Na, Nt, Nt_max
  
  integer, parameter :: NUMCOMPONENTS = 4
  
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
      
      datindex  = 0
      lineindex = 0
      
      do
        read(87, iostat=ios) itime, time_fs
	
	lineindex = lineindex + 1
	
	if (ios .ne. 0) exit
	
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
              read (87) dummy(1:NUMCOMPONENTS, iR, iTheta, iPhi)
            end do
          end do
        end do
	
	if (time_fs >= starttime_fs) then
	  datindex = datindex + 1
            observable(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi, datindex) = dummy(1:NUMCOMPONENTS, 0:NR, 0:NTheta, 0:NPhi)
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
      
      character(8) :: endings(4) = ['_ex.dat', '_ey.dat', '_ez.dat', '_phi.dat']
      
      do c=1,4
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
      
      write(*,'("[STATUS] ", "observable_write_to_vtk: ",a)') fileprefix

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
	      write(24,'(3(x,g25.12))') grid(iR, iTheta, iPhi, 1:3)
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



module spherical_fourier
  implicit none
  
    real*8, allocatable, dimension(:,:,:) :: Rtilda
    real*8, allocatable, dimension(:,:)   :: MM
    real*8, allocatable, dimension(:,:,:) :: P
    complex*16, allocatable, dimension(:,:)   :: E
  
  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Table access function for giving arbitrary l and m
    !>
    !> Calculates the flat index, where an entry for l and m
    !> has to be stored in a one-dimensional array
    !>
    !> @param[in] l
    !> @param[in] m
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function tblinv(l,m)
      use module_data
      implicit none
      integer, intent(in) :: l, m

      if ((l<0) .or. (m<0) .or. (m>l) .or. (l>NTheta/2)) then
        write(*,'("tblinv(l,m) - invalid arguments. l=", I0, " m=", I0)') l, m
      endif

      tblinv = l*(l+1)/2 + m

    end function tblinv
	
	  
    subroutine spherical_fourier_decomposition(freqindex, Scnlm_w)
      use module_data
      implicit none
      integer, intent(in) :: freqindex
      complex*16, intent(out) :: Scnlm_w(1:NUMCOMPONENTS, 1:NR/2, 0:NTheta/2*(NTheta/2+1)/2+NTheta/2)
      integer :: n,l,m
      integer :: iR, iTheta, iPhi, ic
      logical :: omitpoint
      
      complex*16 :: tmp
      
      Scnlm_w(:, :, :) = 0._8
      
        do n = 1,NR/2
          do l = 0,NTheta/2
	    do m=0,NTheta/2
	  
	      if (m <= l) then
                do iR = 0,NR    
                  do iTheta = 0,NTheta 
                    do iPhi = 0,NPhi-1    ! we omit the last entry since is identical to the first one and we want to avoid double counting
		    
		      omitpoint = .false.
		      ! additionally, the following compbinations are degenrate in the spherical coordinate system
		      ! and we only use the repsective first combination
		      ! origin:         iR=0,     iTheta=0..NTheta, iPhi=0..NPhi
		      omitpoint = omitpoint .or. ((iR==0)          .and. (iPhi>0))
		      ! north pole(s):  iR=fixed, iTheta=0,         iPhi=0..NPhi
		      omitpoint = omitpoint .or. ((iTheta==0)      .and. (iPhi>0))
		      ! south pole(s):  iR=fixed, iTheta=NTheta,    iPhi=0..NPhi
		      omitpoint = omitpoint .or. ((iTheta==NTheta) .and. (iPhi>0))
		      
		      if (.not. omitpoint) then
		    
		        tmp = Rtilda(n, l, iR) * MM(l, m) * P(l, m, iTheta) * E(m, iPhi)
		    
                        do ic = 1,NUMCOMPONENTS
		          Scnlm_w(ic, n, tblinv(l, m)) = Scnlm_w(ic, n, tblinv(l, m)) + observable(ic, iR, iTheta, iPhi, freqindex) * tmp
                        end do
	              endif

		    end do
		  end do
	        end do
	      
	      endif

	    end do
	  end do
        end do
      
    end subroutine
    

    subroutine spherical_fourier_decomposition_init()
      use module_data
      use math
      implicit none
      integer :: n, l, m, iR, iTheta, iPhi
      real*8 :: r, r2, rv(3), costheta, Phi
      
      write(*,'("[STATUS] ", "spherical_fourier_decomposition_init")')
      
      allocate(Rtilda(1:NR/2, 0:NTheta/2, 0:NR))
      allocate(MM(0:Ntheta/2, 0:Ntheta/2))
      allocate(P(0:NTheta/2,  0:Ntheta/2, 0:NTheta))
      allocate(E(0:Ntheta/2,  0:NPhi))
      
      ! initialize Rtilda
      do iR=0,NR
        r = spherical_r( grid(NR, NTheta, NPhi, 1:3), r2 )

        do n=1,NR/2
	  do l=0,NTheta/2
	    Rtilda(n, l, iR) = Rfunc(n,l,maxR,r) * r2
	  end do
	end do
      end do
      
      ! initialize M
      do l=0,NTheta/2
        do m=0,NTheta/2
    	    MM(l, m) = MMFunc(l,m)
	end do
      end do
      
      ! initialize P
      do iTheta=0,NTheta
        costheta = spherical_costheta( grid(NR, iTheta, NPhi, 1:3) )
	do l=0,NTheta/2
	  do m=0,NTheta/2
            P(l, m, iTheta) = Pfunc(l,m,costheta) * sqrt(1-costheta*costheta)
	  end do
	end do
      end do

      ! initialize E
      do iPhi=0,NPhi
        phi = spherical_phi( grid(NR, NTheta, iPhi, 1:3) )
        do m=0,NTheta/2
	  E(m, iPhi) = Efunc(m,Phi)
	end do
      end do        
      
    end subroutine
    
    
    subroutine spherical_fourier_decomposition_free()
      implicit none
      
      if (allocated(Rtilda)) deallocate(Rtilda)
      if (allocated(MM))     deallocate(MM)
      if (allocated(P))      deallocate(P)
      if (allocated(E))      deallocate(E)
    
    end subroutine
    
    
    subroutine spherical_fourier_decomposition_for_all(component, filename)
      use progress_bar
      use module_data
      implicit none
      complex*16 :: Scnlm_w(1:NUMCOMPONENTS, 1:NR/2, 0:NTheta/2*(NTheta/2+1)/2+NTheta/2)
      integer, intent(in) :: component
      character(*), intent(in) :: filename
      character*11 :: formatstring = '(????g15.5)'
      
      integer :: i
      
      call spherical_fourier_decomposition_init()

      write(*,'("[STATUS] ", "spherical_fourier_decomposition_for_all, component=",I0," - Wang basis")') component
      
      write(formatstring(2:5),'(I4.4)') 2*(NR/2*(NTheta/2+1)*(NTheta/2+1)) + 1 ! factor 2 since we store real and imaginary part, additional +1 for frequency in first column

      open(24,file=trim(filename),status='unknown',position='rewind',action='write')

      do i=1,Nt/2 ! second half of spectrum is boring anyway
        call progress(i, Nt/2)
        call spherical_fourier_decomposition(i, Scnlm_w)
        write(24,formatstring) abscissa(i), Scnlm_w(component, :, :)
      end do
      
      close(24)

      call spherical_fourier_decomposition_free()
      
      write(*,*) ! progress bar line break

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
      
      if (write_vtk) then
        call system('mkdir -p '//trim(dirname_in)//'/field_vtk')
        call observable_write_to_vtk(trim(dirname_in)//'/field_vtk/spherical_fields')
      endif
      
      call observable_fourier_forward()
      call abscissa_convert_to_frequencies()

      if (write_vtk) then
         call observable_write_to_vtk(trim(dirname_in)//'/field_vtk/spherical_fields_fft')
      endif
      
      call observable_spectrum_to_file(trim(dirname_in)//'/field_spherical_spectrum')
      
      call spherical_fourier_decomposition_for_all(4, trim(dirname_in)//'/field_spherical_fourier_coeffs_phi.dat')

      call spherical_fourier_decomposition_for_all(1, trim(dirname_in)//'/field_spherical_fourier_coeffs_ex.dat')

      call spherical_fourier_decomposition_for_all(2, trim(dirname_in)//'/field_spherical_fourier_coeffs_ey.dat')

      call spherical_fourier_decomposition_for_all(3, trim(dirname_in)//'/field_spherical_fourier_coeffs_ez.dat')

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
