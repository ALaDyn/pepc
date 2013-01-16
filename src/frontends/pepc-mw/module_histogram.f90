!>
!> provides functions for producing histograms 
!> of 3D vectorial data with variable number of components
!>
module module_histogram
  implicit none
  save
  private
  
  public write_histogram
  public fourdhist_init
  public fourdhist_add
  public fourdhist_finalize
  
  real*8, allocatable :: fourdhist(:,:)
  integer :: fourdhist_counter
  integer :: fourdhist_nvals
  
  contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fourdhist_init(nvals)
    use module_debug
    implicit none
    integer, intent(in) :: nvals
    
    fourdhist_nvals = nvals
    
    allocate(fourdhist(1:fourdhist_nvals,1:4))
    fourdhist_counter = 0
  end subroutine
  
  subroutine fourdhist_add(vec, mag2)
    use module_debug
    implicit none
    real*8, intent(in) :: vec(1:3), mag2
    
    fourdhist_counter = fourdhist_counter + 1
    
    if (fourdhist_counter <= fourdhist_nvals) then
      fourdhist(fourdhist_counter,1:3) = vec(1:3)
      fourdhist(fourdhist_counter,  4) = sqrt(mag2)
    else
      DEBUG_ERROR(*,'Added too many entries to histogram')
    endif
  end subroutine
  
  subroutine fourdhist_finalize(sigma, filename, my_rank, comm, state)
    implicit none
    character(*), intent(in) :: filename, state
    integer, intent(in) :: comm, my_rank
    real*8, intent(in) :: sigma
    
    call write_histogram(fourdhist(1:fourdhist_counter,1:4), sigma*[1., 1., 1., 1.], 4, fourdhist_counter, [.True., .True., .True., .False.], filename, my_rank, comm, state)
    deallocate(fourdhist)    
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
    real*8, intent(in)  :: sigma(1:ncols)
    real*8, intent(out) :: binwidth(1:ncols)
    integer, intent(in) :: ncols, nvals_tot
    
    binwidth(1:ncols) = (3.49/((1.*nvals_tot)**(1./3.))) * sigma(1:ncols)   
        
  end subroutine
  
  
  subroutine write_histogram(rawdata, sigma, ncols, nvals, correctmean, filename, my_rank, comm, state)
    implicit none
    include 'mpif.h'
    real*8, intent(in)  :: rawdata(1:nvals,1:ncols)
    real*8, intent(in)  :: sigma(1:ncols)
    integer, intent(in) :: ncols, nvals
    logical, intent(in) :: correctmean(1:ncols)
    integer, intent(in) :: my_rank, comm
    character(*), intent(in) :: filename, state
    integer, parameter :: filept = 97
    
    integer :: ierr, c, i, idx
    real*8 :: binwidth(1:ncols)
    real*8 :: maxvals(1:ncols), minvals(1:ncols), mean(1:ncols)
    integer :: numbins(1:ncols), nvals_tot, nbins
    integer, allocatable :: counter(:,:)
    
    call get_limits(rawdata, ncols, nvals,     maxvals, minvals, mean, nvals_tot, comm)
    call get_binwidth(sigma, ncols, nvals_tot, binwidth)
    
    where (correctmean)
      minvals = minvals - mean
      maxvals = maxvals - mean
    end where
    
    numbins(1:ncols) = ceiling( (maxvals(1:ncols) - minvals(1:ncols)) / binwidth(1:ncols) )
    nbins            = maxval(numbins)
    
    ! data field for storing the results
    allocate(counter(1:ncols,nbins))
    counter(:,:) = 0
    
    ! local sums
    do c=1,ncols
      do i=1,nvals
        if (correctmean(c)) then
          idx             = floor( ( (rawdata(i, c) - mean(c)) - minvals(c)) / binwidth(c) ) + 1
        else
          idx             = floor( ( (rawdata(i, c)          ) - minvals(c)) / binwidth(c) ) + 1
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
      write(filept,*) '# sigma',   sigma(1:ncols)
      write(filept,*) '# nvals_tot', nvals_tot
      write(filept,*) '# nbins', nbins
      write(filept,*) '# ncols', ncols
      write(filept,*) '# binwidth', binwidth
      write(filept,*) '#'
        
      do i=1,nbins
        write(filept,'(i0)', advance='no') i
        do c=1,ncols
          write(filept,'(2(1x,1pe20.10))',advance='no') minvals(c)+(i-0.5)*binwidth(c), counter(c,i) / (1.*nvals_tot*binwidth(c))
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
      write(filept,*) sigma(1:ncols)
      write(filept,*) nvals_tot
      write(filept,*) nbins
      write(filept,*) ncols
      write(filept,*) binwidth
      close(filept)
    endif


    deallocate(counter)  
    
  end subroutine


end module module_histogram
