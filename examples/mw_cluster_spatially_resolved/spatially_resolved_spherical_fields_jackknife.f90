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
      
      open(87,file=trim(filename),status='old',position='rewind',action='read')

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
      character*13 :: formatstring
      integer :: ios
      integer :: datindex

      open(87,file=trim(filename),status='old',position='rewind',action='read', iostat=ios)
      
      if (ios .eq. 0) then
            
        write(formatstring,'("(",I5.5,"g15.5)")') Nev + 1
            
        datindex = 0
      
        do
          datindex = datindex + 1
	
          read(87,formatstring,iostat=ios) rawdata_all(index, datindex, :)
	
	  if (ios .ne. 0) exit
	
        end do
      
        close(87)

        valid(index)   = (datindex == Nt) ! TODO: validity is only checked through comparing the number of lines in data tables
        nvalues(index) = datindex
      else
        valid(index)   = .false.
	nvalues(index) = 0
      endif
      
    end subroutine
    
    subroutine write_data(filename_avg, filename_stddev)
      use progress_bar
      implicit none
      character(*), intent(in) :: filename_avg, filename_stddev
      character*24 :: formatstring
      integer :: i
      
      write(*,'("[STATUS] ", "write_data")')
      
      write(formatstring,'("(g25.12,",I5.5,"(x,g25.12))")') Nev
      
      call progress(0,2*(Nt-1))
      
      open(24,file=trim(filename_avg),status='unknown',position='rewind',action='write')
      do i=1,Nt-1
        call progress(i,2*(Nt-1))
        write(24, formatstring) theta_hat(i,:)
      end do
      close(24)
      
      open(24,file=trim(filename_stddev),status='unknown',position='rewind',action='write')
      do i=1,Nt-1
        call progress(i+Nt-1,2*(Nt-1),.true.)
        write(24, formatstring) stddev(i,:)
      end do
      close(24)

      write(*,*) ! progress bar line break
   end subroutine
  
end module


module computations
  use data
  use progresS_bar
  implicit none
  
  logical, parameter :: showprogress = .false.
  
  contains
  
  subroutine jackknife()
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
    
    call progress(0,Nfiles)

    do i=1,Nfiles
      call progress(i,Nfiles,valid(i))
      
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
    
    theta_i   = theta_i   / (1._8*(nvalid-1))

    do i=1,nvalid
      ! jackknife-average
      theta_dot(:,:) = theta_dot(:,:) + theta_i(i,:,:)
    end do
    
    theta_hat = theta_hat / (1._8*nvalid)
    theta_dot = theta_dot / (1._8*nvalid)
    
    ! correct bias and store bias-corrected average in theta_hat (now theta_tilda in our notation)
    theta_hat = theta_hat - ((nvalid - 1) * (theta_dot - theta_hat))
    
    ! compute variance/stddev
    do i=1,nvalid
      stddev = stddev + ((theta_i(i,:,:) - theta_dot(:,:))*(theta_i(i,:,:) - theta_dot(:,:)))
    end do
    
    stddev = sqrt(stddev)
    
    write(*,*) ! progress bar line break
        
  end subroutine
   
  
end module


program spatially_resolved_spherical_fields_jackknife
  use data
  use computations
  use progress_bar
  implicit none
  integer :: argc, i, f
  character*256, filename_in, dirname_in
  integer :: numlines
  
  character*256 :: filelist(8)
  
  filelist = ['field_spherical_fourier_coeffs_wang_phi',    &
              'field_spherical_fourier_coeffs_raitza_phi',  &
	      'field_spherical_fourier_coeffs_wang_ex',     &
              'field_spherical_fourier_coeffs_raitza_ex',   &
	      'field_spherical_fourier_coeffs_wang_ey',     &
              'field_spherical_fourier_coeffs_raitza_ey',   &
	      'field_spherical_fourier_coeffs_wang_ez',     &
              'field_spherical_fourier_coeffs_raitza_ez']
  
  argc = command_argument_count()

  if (argc < 1) then
    write(*,*) 'Call with directories to be processed as arguments. Exiting.'
    stop
  endif
  
  do f=1,size(filelist)
    write(*,'(/,/,"Processing file",x,a)') trim(filelist(f))
  
    ! we use the number of lines in the first file as array size indicator
    call get_command_argument(argc, dirname_in)
  
    filename_in = trim(dirname_in)//'/'//trim(filelist(f))//'.dat'
  
    call count_num_lines(trim(filename_in), numlines)
  
    call allocate_data(argc, numlines, 81) ! TODO: determine the parameters automatically
    
    call progress(0,argc)

    do i=1,argc
  
      call get_command_argument(i, dirname_in)
      filename_in = trim(dirname_in)//'/'//trim(filelist(f))//'.dat'
    
      call load_data(i, trim(filename_in))
      
      call progress(i, argc, valid(i))
    end do
  
    write(*,*) ! progress bar line break
  
    call jackknife()
    call write_data('./'//trim(filelist(f))//'_jackknife.dat', './'//trim(filelist(f))//'_jackknife_stddev.dat')

    call free_data()
    
  end do
 

end program
