!>
!> provides functions for producing histograms 
!> of 3D vectorial data with variable number of components
!>
module module_histogram
  implicit none
  save
  private
  
  public write_histogram
  
  public t_threedhist  
  
  type t_threedhist
    real*8, private, allocatable :: raw(:,:)
    integer, private :: counter
    integer, private :: vals
    logical, public :: initialized = .false.
  contains
    procedure :: init     => threedhist_init
    procedure :: add      => threedhist_add
    procedure :: dump     => threedhist_dump
    procedure :: finalize => threedhist_finalize
  end type
  
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine threedhist_init(self, nvals)
    use module_debug
    implicit none
    class(t_threedhist) :: self
    integer, intent(in) :: nvals
    
    self%vals = nvals
    self%counter = 0
    
    allocate(self%raw(1:self%vals,1:3))
    
    self%initialized = .true.
  end subroutine
  
  subroutine threedhist_add(self, vec)
    use module_debug
    implicit none
    class(t_threedhist) :: self
    real*8, intent(in) :: vec(1:3)
    
    if (.not. self%initialized) then
      DEBUG_ERROR(*,'Histogram is not initialized.')
    endif
    
    self%counter = self%counter + 1
    
    if (self%counter <= self%vals) then
      self%raw(self%counter, 1:3) = vec(1:3)
    else
      DEBUG_ERROR(*,'Added too many entries to histogram')
    endif
  end subroutine
  
  subroutine threedhist_dump(self, sigmavals, sigmabins, filename, my_rank, comm, state)
    use module_debug
    implicit none
    class(t_threedhist) :: self
    character(*), intent(in) :: filename, state
    integer, intent(in) :: comm, my_rank
    real*8, intent(in) :: sigmavals, sigmabins
    
    if (.not. self%initialized) then
      DEBUG_ERROR(*,'Histogram is not initialized.')
    endif

    call write_histogram(self%raw(1:self%counter,1:3), sigmavals, sigmabins, 3, self%counter, [.True., .True., .True., .False.], filename, my_rank, comm, state)
  end subroutine
 
  subroutine threedhist_finalize(self)
    implicit none
    class(t_threedhist) :: self

    if (self%initialized) then
      self%initialized = .false.
      deallocate(self%raw)
    endif
    
  end subroutine
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  subroutine get_limits(rawdata, ncols, nvals, maxvals, minvals, mean, nvals_tot, comm)
    implicit none
    include 'mpif.h'
    real*8, intent(in)   :: rawdata(1:nvals,1:ncols)
    real*8, intent(out)  :: maxvals(1:ncols), minvals(1:ncols), mean(1:ncols)
    integer, intent(out) :: nvals_tot
    integer, intent(in)  :: ncols, nvals
    integer, intent(in)  :: comm
    integer :: ierr, i
    
    ! find local minima and maxima
    do i=1,ncols
      maxvals(i) = maxval(rawdata(:,i))
      minvals(i) = minval(rawdata(:,i))
      mean(i)    = sum(rawdata(:,i))
    end do
    
    nvals_tot = nvals
    
    ! find global minima and maxima
    call MPI_ALLREDUCE(MPI_IN_PLACE, maxvals, ncols, MPI_REAL8,   MPI_MAX, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, minvals, ncols, MPI_REAL8,   MPI_MIN, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, mean,    ncols, MPI_REAL8,   MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, nvals_tot,   1, MPI_INTEGER, MPI_SUM, comm, ierr)
    
    mean(1:ncols) = mean(1:ncols) / nvals_tot
    
  end subroutine
  

  !> determines optimal bin width according to
  !>       David W. Scott: On optimal and data-based histogram. 
  !>                       Biometrika. 3, Nr. 66, 1979, p. 605-610
  !>                       doi:10.1093/biomet/66.3.605
  subroutine get_binwidth(sigma, ncols, nvals_tot, binwidth)
    implicit none
    real*8, intent(in)  :: sigma
    real*8, intent(out) :: binwidth
    integer, intent(in) :: ncols, nvals_tot
    
    binwidth = (3.49/((1.*nvals_tot)**(1./3.))) * sigma   
        
  end subroutine
  
  
  subroutine write_histogram(rawdata, sigmavals, sigmabins, ncols, nvals, correctmean, filename, my_rank, comm, state)
    implicit none
    include 'mpif.h'
    real*8, intent(in)  :: rawdata(1:nvals,1:ncols)
    real*8, intent(in)  :: sigmavals, sigmabins
    integer, intent(in) :: ncols, nvals
    logical, intent(in) :: correctmean(1:ncols)
    integer, intent(in) :: my_rank, comm
    character(*), intent(in) :: filename, state
    integer, parameter :: filept = 97
    
    integer :: ierr, c, i, idx
    real*8 :: binwidth
    real*8 :: maxvals(1:ncols), minvals(1:ncols), mean(1:ncols)
    integer :: numbins(1:ncols), nvals_tot, nbins
    integer, allocatable :: counter(:,:)
    
    call get_limits(rawdata, ncols, nvals,     maxvals, minvals, mean, nvals_tot, comm)
    call get_binwidth(sigmabins, ncols, nvals_tot, binwidth)
    
    where (correctmean)
      minvals = minvals - mean
      maxvals = maxvals - mean
    end where
    
    ! round minvals downwards towards next multiple of binwidth - with constant binwidth (i.e. constant sigmabins) bins will be reproducible now
    minvals = floor(minvals/binwidth)*binwidth
    
    numbins(1:ncols) = ceiling( (maxvals(1:ncols) - minvals(1:ncols)) / binwidth )
    nbins            = maxval(numbins)
    
    ! data field for storing the results
    allocate(counter(1:ncols,nbins))
    counter(:,:) = 0
    
    ! local sums
    do c=1,ncols
      do i=1,nvals
        if (correctmean(c)) then
          idx             = floor( ( (rawdata(i, c) - mean(c)) - minvals(c)) / binwidth ) + 1
        else
          idx             = floor( ( (rawdata(i, c)          ) - minvals(c)) / binwidth ) + 1
        endif
          
        counter(c, idx) = counter(c,idx) + 1
      end do
    end do


    ! global sum    
    call MPI_ALLREDUCE(MPI_IN_PLACE, counter, ncols*nbins, MPI_INTEGER, MPI_SUM, comm, ierr)

    ! write data to file
    if (my_rank == 0) then
    
      open(filept, FILE=trim(filename), STATUS='UNKNOWN', POSITION = 'REWIND')
      
      write(filept,*) '# ', state
      write(filept,*) '# minvals', minvals(1:ncols)
      write(filept,*) '# maxvals', maxvals(1:ncols)
      write(filept,*) '# mean   ', mean(1:ncols)
      write(filept,*) '# numbins', numbins(1:ncols)
      write(filept,*) '# sigmavals', sigmavals
      write(filept,*) '# nvals_tot', nvals_tot
      write(filept,*) '# nbins', nbins
      write(filept,*) '# ncols', ncols
      write(filept,*) '# binwidth', binwidth
      write(filept,*) '# sigmabins', sigmabins
      write(filept,*) '#'
        
      do i=1,nbins
        write(filept,'(i0)', advance='no') i
        do c=1,ncols
          write(filept,'(2(1x,1pe20.10))',advance='no') minvals(c)+(i-0.5)*binwidth, counter(c,i) / (1.*nvals_tot*binwidth)
        end do
        write(filept,'(/)', advance='no')
      end do
      
      close(filept)

      ! write json file with parameters
      open(filept, FILE=trim(filename)//'.params', STATUS='UNKNOWN', POSITION = 'REWIND')
      write(filept,*) state
      write(filept,*) minvals(1:ncols)
      write(filept,*) maxvals(1:ncols)
      write(filept,*) mean(1:ncols)
      write(filept,*) numbins(1:ncols)
      write(filept,*) sigmavals
      write(filept,*) nvals_tot
      write(filept,*) nbins
      write(filept,*) ncols
      write(filept,*) binwidth
      write(filept,*) sigmabins
      close(filept)
    endif


    deallocate(counter)  
    
  end subroutine


end module module_histogram
