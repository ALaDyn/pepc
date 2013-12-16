module FFTW3
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
end module FFTW3



module data
  use, intrinsic :: iso_c_binding
  implicit none
  
  integer :: Nt, Nev, Nfiles
  
  real*8, allocatable, dimension(:,:,:) :: rawdata_all
  logical, allocatable, dimension(:)  :: valid
  integer, allocatable, dimension(:)  :: nvalues
  real*8, allocatable, dimension(:,:) :: theta_hat, stddev
  
  contains
  
    subroutine allocate_data(Nfiles_, Nt_, Nev_)
      implicit none
      integer, intent(in) :: Nfiles_, Nt_, Nev_
      
      Nt     = Nt_
      Nev    = Nev_
      Nfiles = Nfiles_
      
      allocate(rawdata_all(Nfiles, Nt, Nev + 1)) ! first column is left for frequency argument
      rawdata_all = 0.
      
      allocate(valid(NFiles))
      allocate(nvalues(NFiles))
      
      allocate(theta_hat(Nt, Nev + 1))
      allocate(stddev(Nt, Nev + 1))
      
    end subroutine
    
    
    subroutine free_data()
      implicit none
      
      deallocate(rawdata_all)
      deallocate(valid)
      deallocate(nvalues)
      deallocate(theta_hat)
      deallocate(stddev)
      
    end subroutine
    
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
    
  
    subroutine load_data(index, filename)
      implicit none
      character(*), intent(in) :: filename
      integer, intent(in) :: index
      character*24 :: formatstring
      integer :: ios
      integer :: datindex

      datindex = 0
      
      open(87,file=trim(filename),status='old',position='rewind',action='read',iostat=ios)
            
      if (ios .eq. 0) then

        write(formatstring,'("(g25.12,",I5.5,"(x,g25.12))")') Nev
            
        do
          datindex = datindex + 1
	  
	  if (datindex > Nt) then
	    datindex = 0
	    exit
	  endif
	
          read(87,formatstring,iostat=ios) rawdata_all(index, datindex, :)
	
  	  if (ios .ne. 0) exit
	
        end do
      
        close(87)
      endif

      valid(index)   = (datindex == Nt) ! TODO: validity is only checked through comparing the number of lines in data tables
      nvalues(index) = datindex
      
    end subroutine
    
    subroutine write_data(filename_avg, filename_stddev)
      implicit none
      character(*), intent(in) :: filename_avg, filename_stddev
      character*24 :: formatstring
      integer :: i
      
      write(*,'("[STATUS] ", "write_data")')
      
      write(formatstring,'("(g25.12,",I5.5,"(x,g25.12))")') Nev
      
      open(24,file=trim(filename_avg),status='unknown',position='rewind',action='write')
      do i=1,Nt-1
        write(24, formatstring) theta_hat(i,:)
      end do
      close(24)
      
      open(24,file=trim(filename_stddev),status='unknown',position='rewind',action='write')
      do i=1,Nt-1
        write(24, formatstring) stddev(i,:)
      end do
      close(24)

   end subroutine
  
end module


module computations
  use data
  implicit none
  
  logical, parameter :: showprogress = .false.
  
  contains
  
  subroutine jackknife()
    use progress_bar
    implicit none
    integer :: nvalid, i, j
    real*8, dimension(Nt, Nev + 1) :: theta_dot
    real*8, dimension(Nfiles, Nt, Nev + 1) :: theta_i
    
    write(*,'("[STATUS] ", "jackknife")')

    nvalid = 0
    
    theta_hat   = 0.
    theta_dot   = 0.
    theta_i     = 0.
    stddev      = 0.
    
    call progress(0, Nfiles)
    
    do i=1,Nfiles
      call progress(i, Nfiles)

      if (valid(i)) then
         nvalid = nvalid + 1
      
        ! classical average
        theta_hat(:,:) = theta_hat(:,:) + rawdata_all(i, :, :)
	
	! every jackknife-dataset is made form (nvalid-1) original datasets
	do j=1,Nfiles
	  if ((j .ne. i) .and. (valid(j))) then
            theta_i(nvalid,:,:) = theta_i(nvalid,:,:) + rawdata_all(j, :, :)
	  endif
	end do
      endif
    end do

    write(*,*) ! progress bar linebreak
    
    theta_i   = theta_i   / (1._8*(nvalid-1))

    call progress(0, nvalid)
    do i=1,nvalid
      call progress(i, nvalid)
      ! jackknife-average
      theta_dot(:,:) = theta_dot(:,:) + theta_i(i,:,:)
    end do
    write(*,*) ! progress bar linebreak
    
    theta_hat = theta_hat / (1._8*nvalid)
    theta_dot = theta_dot / (1._8*nvalid)
    
    ! correct bias and store bias-corrected average in theta_hat (now theta_tilda in our notation)
    theta_hat = theta_hat - ((nvalid - 1) * (theta_dot - theta_hat))
    
    ! compute variance/stddev
    call progress(0, nvalid)
    do i=1,nvalid
      call progress(i, nvalid)
      stddev = stddev + ((theta_i(i,:,:) - theta_dot(:,:))*(theta_i(i,:,:) - theta_dot(:,:)))
    end do
    write(*,*) ! progress bar linebreak
    
    stddev = sqrt(stddev)
        
  end subroutine
   
  
end module


program spatially_resolved_ccf_jackknife
  use data
  use computations
  implicit none
  integer :: argc, i
  character*256, filename_in, dirname
  integer :: numlines
  
  argc = command_argument_count()

  if (argc < 1) then
    write(*,*) 'Call with file to be processed as first argument. Exiting.'
    stop
  endif
  
  ! we use the number of lines in the last file as array size indicator
  call get_command_argument(argc, dirname)
  filename_in = trim(dirname) // '/momentum_electrons_acf.dat'
  call count_num_lines(trim(filename_in), numlines)
  
  call allocate_data(argc, numlines, 4) ! TODO: determine the parameters automatically

  do i=1,argc
    call get_command_argument(i, dirname)
    filename_in = trim(dirname) // '/momentum_electrons_acf.dat'
    
    write(*,'("Reading data from file ",a)') trim(filename_in)

    call load_data(i, trim(filename_in))
  end do
  
  call jackknife()
  call write_data('./momentum_electrons_acf_jackknife.dat', './momentum_electrons_acf_jackknife_stddev.dat')
  write(*,*) argc, valid
  call free_data()
 

end program
