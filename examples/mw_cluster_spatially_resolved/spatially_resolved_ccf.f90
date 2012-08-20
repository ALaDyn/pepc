module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3



module data
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: RAWDATA_PX = 1
  integer, parameter :: RAWDATA_PY = 2
  integer, parameter :: RAWDATA_PZ = 3
  integer, parameter :: RAWDATA_P  = 4
  integer, parameter :: RAWDATA_N  = 5
  
  integer, parameter :: NUM_RAWDATA = RAWDATA_N


  integer :: Nc, Na, Nt, Nt_max
  
  integer :: NR, NTheta, NPhi
  real*8 :: maxR

  ! indices: spatial component (c), location index (a), time/frequency index (t/w)
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:) :: observable
  ! indices: location index (a), location index (a'), time/frequency index (t/w)
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:) :: observableCCF
  real*8, allocatable, dimension(:) :: abscissa ! will hold time/frequency values
  
  contains
  
    subroutine allocate_data()
      implicit none
      
      allocate(observableCCF(Na, Na, 2*Nt))
      
    end subroutine allocate_data
    
    
    subroutine free_data()
      implicit none
      
      if (allocated(observable))    deallocate(observable)
      if (allocated(observableCCF)) deallocate(observableCCF)
      if (allocated(abscissa))      deallocate(abscissa)
      
    end subroutine free_data
    
    
    logical function load_data(filename, fields, starttime_fs)
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: fields(:)
      real*8, intent(in) :: starttime_fs
      integer :: ncolumns
      integer :: itime
      real*8 :: time_fs
      real*8, allocatable, dimension(:) :: rawdata
      character*26 :: formatstring
      integer :: ios
      integer :: datindex, a, c

      open(87,file=trim(filename),status='old',position='rewind',action='read')
      
      read(87, '("# Nt     = ",   i15)') Nt_max
      read(87, '("# maxR   = ", g15.5)') maxR
      read(87, '("# NR     = ",   i15)') NR
      read(87, '("# NTheta = ",   i15)') NTheta
      read(87, '("# NPhi   = ",   i15)') NPhi
      
      Nc = size(fields)
      Na = NR * NTheta * NPhi
      ncolumns = NUM_RAWDATA * Na ! one entry of each single observable per position

      write(*, '("# Nt_max   = ",   i15)') Nt_max
      write(*, '("# maxR     = ", g15.5)') maxR
      write(*, '("# NR       = ",   i15)') NR
      write(*, '("# NTheta   = ",   i15)') NTheta
      write(*, '("# NPhi     = ",   i15)') NPhi
      write(*, '("# Nc       = ",   i15)') Nc
      write(*, '("# Na       = ",   i15)') Na
      write(*, '("# ncolumns = ",   i15)') ncolumns
      
      ! skip 4 lines in input file
      read(87,*)
      read(87,*)
      read(87,*)
      read(87,*)
      
      write(formatstring, '("(i10,g25.12,",I5.5,"(g25.12))")') ncolumns
      allocate(rawdata(ncolumns))
      allocate(observable(Nc, Na, 2*Nt_max))
      allocate(abscissa(Nt_max))
      observable = 0. ! this automatically ensures padding, i.e. the second half of the array is zero
      
      datindex = 0
      
      do
        read(87,formatstring,iostat=ios) itime, time_fs, rawdata
	
	if (ios .ne. 0) exit
      
        !write(*,*) itime, time_fs
	
	if (time_fs >= starttime_fs) then
	  datindex = datindex + 1
	  
	  abscissa(datindex) = time_fs
	  
	  do a=1,Na
	    do c=1,Nc
              observable(c, a, datindex) = rawdata((a-1)*NUM_RAWDATA+fields(c))
	    end do
	  end do
	endif
      end do
      
      deallocate(rawdata)

      close(87)
      
      Nt = datindex ! store number of actually available timesteps
      write(*, '("# Nt       = ",   i15)') Nt

      load_data = (Nt > 2)

      if (.not. load_data) then
        write(*,'(a,/,a)') 'Did not find a sufficient number of timesteps.', 'Reduce starttime_fs or check your data'
      endif
      
    end function
  
end module

module computations
  use data
  implicit none
  
  logical, parameter :: showprogress = .false.
  
  contains
  
  subroutine observable_fourier_forward()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: a,c
    
    write(*,'("[STATUS] ", "observable_fourier_forward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)

    do c=1,Nc
      do a=1,Na
        tmparray(1:2*Nt) = observable(c,a,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observable(c,a,1:2*Nt) = tmparray(1:2*Nt)
      end do
    end do
    
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
    implicit none
    integer :: a, t
    
    write(*,'("[STATUS] ", "observable_spatial_average")')

    do a=1,Na
      do t=1,2*Nt
        observable(1,a,t) = 1./3. * sum(observable(1:3, a, t))
      end do
    end do
    
  end subroutine


  subroutine observable_cross_correlate()
    implicit none
    integer :: a, aprime
    
    write(*,'("[STATUS] ", "observable_cross_correlate")')

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a, aprime)
    do a=1,Na
      if (showprogress) write(*, '(15x,i0,"/",i0)') a, Na
      do aprime = 1,Na
        observableCCF(a, aprime, 1:2*Nt) = observable(1,a,1:2*Nt) * conjg(observable(1,aprime,1:2*Nt))
      end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine
  
  
  subroutine observableCCF_fourier_backward()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: a, aprime
    
    write(*,'("[STATUS] ", "observableCCF_fourier_backward")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_FORWARD, FFTW_ESTIMATE)
    
    do a=1,Na
      if (showprogress) write(*, '(15x,i0,"/",i0)') a, Na
      do aprime=1,Na
        tmparray(1:2*Nt) = observableCCF(a,aprime,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observableCCF(a,aprime,1:2*Nt) = tmparray(1:2*Nt) / (2._8*Nt)
      end do
    end do
    
    call fftw_destroy_plan(fftw_plan)
   
 end subroutine   


  subroutine observableCCF_laplace_transform()
    use, intrinsic :: iso_c_binding
    use FFTW3
    implicit none
    
    type(C_PTR) :: fftw_plan
    complex(C_DOUBLE_COMPLEX), dimension(2*Nt) :: tmparray
    integer :: a, aprime
    
    write(*,'("[STATUS] ", "observableCCF_laplace_transform")')

    fftw_plan = fftw_plan_dft_1d(size(tmparray), tmparray, tmparray, FFTW_BACKWARD, FFTW_ESTIMATE)
    
    do a=1,Na
      if (showprogress) write(*, '(15x,i0,"/",i0)') a, Na
      do aprime=1,Na
        tmparray(1:2*Nt) = observableCCF(a,aprime,1:2*Nt)
      
        ! perform FFT
        call fftw_execute_dft(fftw_plan, tmparray, tmparray)
      
        observableCCF(a,aprime,1:2*Nt) = tmparray(1:2*Nt)
      end do
    end do
    
    call fftw_destroy_plan(fftw_plan)
    
  end subroutine   
    
    
  subroutine identify_eigenvalues_for_all_frequencies(num_ev,filename_ev)
    implicit none
    
    integer :: t,a,aprime
    double precision, dimension(Na*(Na+1)/2) :: tmparray
    double precision, dimension(Na) :: D
    double precision, dimension(Na-1) :: E, TAU
    integer :: INFO
    double precision :: VL = 0., VU=0., ABSTOL=1.0
    integer, intent(in) :: num_ev ! suche nach groessten 10 Eigenwerten
    integer :: num_ev_found
    double precision, dimension(Na) :: eigenvalues
    double precision, dimension(Na, num_ev) :: eigenvectors
    integer, dimension(2*num_ev) :: ISUPPZ
    double precision, allocatable :: WORK(:)
    integer, allocatable :: IWORK(:)
    double precision :: TESTWORK(1)
    integer :: TESTIWORK(1)
    logical :: TRYAC = .TRUE.
    integer :: LWORK, LIWORK
    character*24 :: formatstring
    character(*), intent(in) :: filename_ev

    write(*,'("[STATUS] ", "identify_eigenvalues_for_all_frequencies")')
    write(formatstring,'("(g25.12,",I5.5,"(x,g25.12))")') num_ev
    open(24,file=trim(filename_ev),status='unknown',position='rewind',action='write')
    
    do t=1,Nt
   
      do a=1,Na
        do aprime=a,Na                         ! we are only interested inthe real part here
	  tmparray(a + (aprime-1)*aprime/2) = (real(observableCCF(a,aprime,t)) + real(observableCCF(aprime,a,t)))/2._8 ! average over both (symmetric) contributions
	end do
      end do
	  
      call DSPTRD( 'U', Na, tmparray, D, E, TAU, INFO )
	  
      if (INFO .ne. 0) write (*,'("Error while calling ZHPTRD:", I0)') INFO
  
      ! query size of work arrays
      LWORK  = -1
      LIWORK = -1
      call DSTEMR( 'V', 'I', Na, D, E, VL, VU, Na-num_ev+1, Na, num_ev_found, eigenvalues, eigenvectors, Na, num_ev, ISUPPZ, TRYAC, TESTWORK, LWORK, TESTIWORK, LIWORK, INFO )
      LWORK  = int( TESTWORK(1))
      LIWORK =     TESTIWORK(1)
      
      allocate( WORK( LWORK))
      allocate(IWORK(LIWORK))

      call DSTEMR( 'N', 'I', Na, D, E, VL, VU, Na-num_ev+1, Na, num_ev_found, eigenvalues, eigenvectors, Na, num_ev, ISUPPZ, TRYAC, WORK, LWORK, IWORK, LIWORK, INFO )
	  ! call with 'V' as first argument for also getting the eigenvectors
      deallocate(WORK)
      deallocate(IWORK)

      if (INFO .ne. 0) write (*,'("Error while calling DSTEGR:", I0)') INFO

      if (showprogress) write(*, '(15x,i0,"/",i0," - Found ",i0,"/",i0," eigenvalues.")') t, Nt, num_ev_found, num_ev

      ! output them in reverse order, i.e. largest first
      write(24,formatstring) abscissa(t), eigenvalues(num_ev:1:-1)
    end do
    
    close(24)
    
  end subroutine
  
  
  subroutine analysis_workflow(filename_in, fields, starttime_fs, num_ev, filename_ev)
    implicit none
    character(*), intent(in) :: filename_in
    character(*), intent(in) :: filename_ev ! filename to output the frequency-dependent eigenvalues to
    integer, intent(in) :: num_ev ! number of eigenvalues to be searched for
    integer, intent(in) :: fields(:)
    real*8, intent(in) :: starttime_fs
    
    if (load_data(filename_in, fields, starttime_fs)) then
      ! call observable_undo_volume_scaling() ! see (18) in [PRE84,036406], not necessary since we already output the observable without volume scaling, compare write_spatially_resolved_data () in module_diagnostics.f90
      call allocate_data()
      call observable_fourier_forward()
      call abscissa_convert_to_frequencies()
      call observable_spatial_average()
      call observable_cross_correlate()
      ! die folgenden beiden Schritte sind zueinander invers, sollten also obsolet sein - machen auch nach Test keinen Unterschied am Ergebnis
      !  call observableCCF_fourier_backward()   ! jetzt steht D_aa'(t) in observableCCF(a,a',t)
      !  call observableCCF_laplace_transform()  ! jetzt steht D_aa'(w) in observableCCF(a,a',w)
      !call observableCCF_use_symmetries() ! TODO: see (15-17) in  [PRE84,036406]
      call identify_eigenvalues_for_all_frequencies(num_ev, filename_ev)
      ! call output eigenvalues()
    endif

    call free_data()

  end subroutine
  
end module


program spatially_resolved_ccf
  use data
  use computations
  implicit none
  integer :: argc, i
  character*256, filename_in, filename_out
  
  argc = command_argument_count()

  if (argc < 1) then
    write(*,*) 'Call with file to be processed as first argument. Exiting.'
    stop
  endif

  do i=1,argc
    call get_command_argument(i, filename_in)
    filename_out = trim(filename_in(1:len(trim(filename_in))-4))//'_momentum_eigenvalues.dat'
    
    write(*,'(2/"Reading data from file ",a,/,"Writing data to file   ", a)') trim(filename_in), trim(filename_out)

    call analysis_workflow(trim(filename_in), [RAWDATA_PX, RAWDATA_PY, RAWDATA_PZ], 100.0_8, 10, trim(filename_out))
  end do

  
  
 

end program
