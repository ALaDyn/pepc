module field_helper
   use module_pepc_kinds
   use module_globals !,only: offset,extent
   implicit none

  character(*), parameter, private :: field_dir = "fields"

  type field_grid_nml_t
    integer, dimension(3) :: n
    real(kind_physics), dimension(3) :: offset, extent
  end type field_grid_nml_t



contains


   subroutine pepc_setup(pepc_pars)
      use encap
      use module_globals, only: FRONTEND_NAME,my_rank,n_ranks,communicator
      use module_pepc
      use module_pepc_types, only: t_particle
      implicit none

      type(pepc_pars_t), intent(out) :: pepc_pars
      character(255)                 :: para_file_name
      logical                        :: para_file_available

!      call pepc_initialize(FRONTEND_NAME, pepc_pars%pepc_comm%mpi_rank, &
!        pepc_pars%pepc_comm%mpi_size, .true., comm = pepc_pars%pepc_comm%mpi_comm)

!      call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)
!      call read_in_params(pepc_nml, para_file_available, para_file_name)

      ! Pass MPI stuff to parameters
      pepc_pars%pepc_comm%mpi_rank      = my_rank
      pepc_pars%pepc_comm%mpi_size      = n_ranks
      pepc_pars%pepc_comm%mpi_comm      = communicator


   end subroutine pepc_setup


  subroutine setup_field_grid(field_grid, pepc_comm)
    !use helper, only: my_rank,n_ranks,root
    use module_utils, only: create_directory
    use encap
    implicit none

    type(field_grid_t), intent(out) :: field_grid
    type(pepc_comm_t), intent(in)   :: pepc_comm
    integer(kind_default)           :: mpi_size, mpi_rank
    type(field_grid_nml_t)          :: field_grid_nml
    integer(kind_particle)          :: ipl, ipg, ix, iy,iz
    character(len = 255)            :: file_name

    call create_directory(trim(folder))
    call create_directory(trim(folder)//trim("fields"))
    call read_in_field_grid_params(field_grid_nml, file_name)

    mpi_rank = int(pepc_comm%mpi_rank, kind = kind(mpi_rank))
    mpi_size = int(pepc_comm%mpi_size, kind = kind(mpi_size))

    field_grid%n = field_grid_nml%n
    field_grid%ntot = product(field_grid%n)
    field_grid%nl = field_grid%ntot / mpi_size
    if (mpi_rank < mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) &
      field_grid%nl = field_grid%nl + 1
    field_grid%offset = field_grid_nml%offset
    field_grid%extent = field_grid_nml%extent
    field_grid%dx = field_grid%extent / field_grid%n


    allocate(field_grid%p(field_grid%nl), &
      field_grid%ne(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
      field_grid%ni(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
      field_grid%nb(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
      field_grid%Ia(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
      field_grid%Ip(field_grid%n(1), field_grid%n(2), field_grid%n(3))   &
!      field_grid%qe(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
!      field_grid%qi(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
!      field_grid%qb(field_grid%n(1), field_grid%n(2), field_grid%n(3)) , &
!      field_grid%vex(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vey(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vez(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vix(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%viy(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%viz(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vbx(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vby(field_grid%n(1), field_grid%n(2), field_grid%n(3)), &
!      field_grid%vbz(field_grid%n(1), field_grid%n(2), field_grid%n(3))  &
                                             
      )


    do ipl = 1, field_grid%nl

      ipg = ipl + min(int(mpi_rank, kind=kind_particle), mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) * &
        (field_grid%ntot / mpi_size + 1) + &
        max(0_kind_particle, mpi_rank - mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) * &
        (field_grid%ntot / mpi_size)

      ix = mod(ipg - 1, field_grid%n(1)) + 1
      iy = mod( (ipg - 1) / field_grid%n(1) , field_grid%n(2) )  + 1
      iz = (ipg - 1) / ( field_grid%n(1) * field_grid%n(2) )  + 1

      field_grid%p(ipl)%x(1) = (ix - 0.5) * field_grid%dx(1) + field_grid%offset(1)
      field_grid%p(ipl)%x(2) = (iy - 0.5) * field_grid%dx(2) + field_grid%offset(2)
      field_grid%p(ipl)%x(3) = (iz - 0.5) * field_grid%dx(3) + field_grid%offset(3)
      field_grid%p(ipl)%label = ipg
      

    end do


  end subroutine setup_field_grid

  
   subroutine read_in_field_grid_params(field_grid_namelist, file_name)
      use encap
      !use helper, only: root
      use module_pepc, only: pepc_read_parameters_from_first_argument
      implicit none

      namelist /field_grid_nml/  offset,extent,n
      character(len = 255), intent(in)    :: file_name
      type(field_grid_nml_t), intent(out) :: field_grid_namelist
      logical                             :: file_available



      !namelist /field_grid_nml/ n, offset, extent

      integer, parameter :: para_file_id = 10

      call pepc_read_parameters_from_first_argument(file_available, file_name)
      if (file_available) then
         open(para_file_id,file=file_name)!,action='read')
         rewind(para_file_id)
         read(para_file_id, NML=field_grid_nml)
         close(para_file_id)

      else
         if(root) write(*,*) " == no param file, using default parameter "

      end if

      field_grid_namelist%n      = n
      field_grid_namelist%offset = offset
      field_grid_namelist%extent = extent

   end subroutine read_in_field_grid_params


  subroutine compute_field(pepc_pars, field_grid, p)
  !subroutine compute_field(field_grid, p)
    !include 'mpif.h'
    use module_pepc     , only: pepc_particleresults_clear, pepc_traverse_tree
    use mpi
    use module_pepc
    use module_tree     , only : tree_allocated
    use module_interaction_specific, only : force_law
    use module_shortcut , only: half,pi,one,prec,zero
    use module_globals  , only: atilde,btilde,etilde,jtilde,phitilde
    use encap
    !use physics_helper
    implicit none

    type(pepc_pars_t)              , intent(in)    :: pepc_pars
    type(field_grid_t)             , intent(inout) :: field_grid
    type(t_particle)  , allocatable, intent(inout) :: p(:)

    integer(kind_particle) :: ipl
    integer(kind_default)  :: mpi_err
    integer, dimension(3)  :: ic
    real(kind_physics)     :: da, rda
    
    call pepc_particleresults_clear(field_grid%p)

    if (.not. tree_allocated(global_tree)) then
      call pepc_grow_and_traverse_for_others(p, field_grid%p, 0)
    else
      call pepc_traverse_tree(field_grid%p)
    endif
    
    
    do ipl = 1,3
        field_grid%p(:)%results%E(ipl)    = field_grid%p(:)%results%E(ipl)    * etilde
!        field_grid%p(:)%results%dxE(ipl)  = field_grid%p(:)%results%dxE(ipl)  * etilde
!        field_grid%p(:)%results%dyE(ipl)  = field_grid%p(:)%results%dyE(ipl)  * etilde
        field_grid%p(:)%results%B(ipl)    = field_grid%p(:)%results%B(ipl)    * btilde
        field_grid%p(:)%results%A(ipl)    = field_grid%p(:)%results%A(ipl)    * atilde
        field_grid%p(:)%results%dxA(ipl)  = field_grid%p(:)%results%dxA(ipl)  * atilde
        field_grid%p(:)%results%dyA(ipl)  = field_grid%p(:)%results%dyA(ipl)  * atilde
        field_grid%p(:)%results%dzA(ipl)  = field_grid%p(:)%results%dzA(ipl)  * atilde
!        field_grid%p(:)%results%dxxA(ipl) = field_grid%p(:)%results%dxxA(ipl) * atilde
!        field_grid%p(:)%results%dxyA(ipl) = field_grid%p(:)%results%dxyA(ipl) * atilde
!        field_grid%p(:)%results%dyyA(ipl) = field_grid%p(:)%results%dyyA(ipl) * atilde
        field_grid%p(:)%results%J(ipl)    = field_grid%p(:)%results%J(ipl)    * jtilde
        field_grid%p(:)%results%Jirr(ipl) = field_grid%p(:)%results%Jirr(ipl) * jtilde   
    enddo
    
    field_grid%p(:)%results%pot    = field_grid%p(:)%results%pot  * phitilde
!    field_grid%p(:)%results%rho    = field_grid%p(:)%results%rho  * rhotilde
    
    field_grid%ne = 0
    field_grid%ni = 0
    field_grid%nb = 0
!    field_grid%qe = 0
!    field_grid%qi = 0
!    field_grid%qb = 0
!    field_grid%vex = zero
!    field_grid%vey = zero
!    field_grid%vez = zero
!    field_grid%vix = zero
!    field_grid%viy = zero
!    field_grid%viz = zero
!    field_grid%vbx = zero
!    field_grid%vby = zero
!    field_grid%vbz = zero
!    field_grid%Jix = zero
!    field_grid%Jiy = zero
!    field_grid%Jiz = zero
!    field_grid%Jbx = zero
!    field_grid%Jby = zero
!    field_grid%Jbz = zero
!    field_grid%Bex = zero
!    field_grid%Bey = zero
!    field_grid%Bez = zero
!    field_grid%Bix = zero
!    field_grid%Biy = zero
!    field_grid%Biz = zero
!    field_grid%Bbx = zero
!    field_grid%Bby = zero
!    field_grid%Bbz = zero
    
    da = field_grid%dx(1)
    if ( field_grid%dx(2) .gt. zero ) da = da*field_grid%dx(2) 
    if ( field_grid%dx(3) .gt. zero ) da = da*field_grid%dx(3) 
    rda = one/da
    
    ! make a spatial histogram of local particle numbers and velocities per species
    do ipl = 1, size(p)
      ! do not count particles outside the grid
      
             
      if (any(p(ipl)%x(1:3) < field_grid%offset(1:3)) .or. &
        any((p(ipl)%x(1:3) - field_grid%offset(1:3)) >= field_grid%extent(1:3))) cycle
        
      ic(force_law:)      = 1  
!      write(*,*) "x-offset: ", p(ipl)%x(1:2) , field_grid%offset(1:2)
       
      ic(1:force_law) = ceiling( (p(ipl)%x(1:force_law) - field_grid%offset(1:force_law) ) / field_grid%dx(1:force_law) ) 
      
      
      if (p(ipl)%label .eq. 1 ) then ! electrons
        field_grid%ne(ic(1), ic(2), ic(3))  = field_grid%ne(ic(1), ic(2), ic(3))  +  1
!        field_grid%qe(ic(1), ic(2), ic(3))  = field_grid%qe(ic(1), ic(2), ic(3))  +  p(ipl)%data%q
!        field_grid%vex(ic(1), ic(2), ic(3)) = field_grid%vex(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(1)/p(ipl)%data%g
!        field_grid%vey(ic(1), ic(2), ic(3)) = field_grid%vey(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(2)/p(ipl)%data%g
!        field_grid%vez(ic(1), ic(2), ic(3)) = field_grid%vez(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(3)/p(ipl)%data%g
        
      else if (p(ipl)%label .eq. 2 ) then ! beam
   
        field_grid%nb(ic(1), ic(2), ic(3))  = field_grid%nb(ic(1), ic(2), ic(3))  +  1
!        field_grid%qb(ic(1), ic(2), ic(3))  = field_grid%qb(ic(1), ic(2), ic(3))  +  p(ipl)%data%q
!        field_grid%vbx(ic(1), ic(2), ic(3)) = field_grid%vbx(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(1)/p(ipl)%data%g
!        field_grid%vby(ic(1), ic(2), ic(3)) = field_grid%vby(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(2)/p(ipl)%data%g
!        field_grid%vbz(ic(1), ic(2), ic(3)) = field_grid%vbz(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(3)/p(ipl)%data%g
         
      else if (p(ipl)%label .eq. 3 ) then ! ions
        field_grid%ni(ic(1), ic(2), ic(3))  = field_grid%ni(ic(1), ic(2), ic(3))  +  1
!        field_grid%qi(ic(1), ic(2), ic(3))  = field_grid%qi(ic(1), ic(2), ic(3))  +  p(ipl)%data%q
!        field_grid%vix(ic(1), ic(2), ic(3)) = field_grid%vix(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(1)/p(ipl)%data%g
!        field_grid%viy(ic(1), ic(2), ic(3)) = field_grid%viy(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(2)/p(ipl)%data%g
!        field_grid%viz(ic(1), ic(2), ic(3)) = field_grid%viz(ic(1), ic(2), ic(3)) +  p(ipl)%data%v(3)/p(ipl)%data%g
        
      end if
    end do

    ! accumulate the histograms
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ne, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ni, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%nb, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
      
!      call mpi_allreduce(MPI_IN_PLACE, field_grid%qe, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%qi, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%qb, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
!
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vex, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vey, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vez, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vix, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%viy, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%viz, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vbx, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vby, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
!    call mpi_allreduce(MPI_IN_PLACE, field_grid%vbz, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
!      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
    call mpi_allreduce(MPI_IN_PLACE, field_grid%Ia, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)  
    call mpi_allreduce(MPI_IN_PLACE, field_grid%Ip, int(field_grid%ntot, kind = kind_default), MPI_KIND_PHYSICS, &
      MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)        

    ! normalize to fluid quantities
    ! total velocity -> mean velocity
!    where (field_grid%ne .gt. 0)
!      field_grid%vex = field_grid%vex / field_grid%ne
!      field_grid%vey = field_grid%vey / field_grid%ne
!      field_grid%vez = field_grid%vez / field_grid%ne
!
!    elsewhere
!      field_grid%vex = zero
!      field_grid%vey = zero
!      field_grid%vez = zero
!
!    end where
!
!    where (field_grid%ni .gt. 0)
!      field_grid%vix = field_grid%vix / field_grid%ni
!      field_grid%viy = field_grid%viy / field_grid%ni
!      field_grid%viz = field_grid%viz / field_grid%ni
!
!    elsewhere
!      field_grid%vix = zero
!      field_grid%viy = zero
!      field_grid%viz = zero
!
!    end where
!    
!    where (field_grid%nb .gt. 0)
!      field_grid%vbx = field_grid%vbx / field_grid%nb
!      field_grid%vby = field_grid%vby / field_grid%nb
!      field_grid%vby = field_grid%vby / field_grid%nb
!    elsewhere
!      field_grid%vbx = zero
!      field_grid%vby = zero
!      field_grid%vbz = zero
!    end where
    
    
    where (field_grid%nb .gt. zero)
          
        field_grid%Ia(:,:,:) = field_grid%Ia(:,:,:)/field_grid%nb(:,:,:)
        
    elsewhere
        
          field_grid%Ia(:,:,:) = zero
          
    end where
    
    
    field_grid%Ip(:,:,:) = field_grid%Ip * rda
    
    ! particle number -> particle density
    field_grid%ne(:,:,:) = field_grid%ne * rda
    field_grid%ni(:,:,:) = field_grid%ni * rda
    field_grid%nb(:,:,:) = field_grid%nb * rda
    
!    field_grid%qe(:,:,:) = field_grid%qe * rda
!    field_grid%qi(:,:,:) = field_grid%qi * rda
!    field_grid%qb(:,:,:) = field_grid%qb * rda
    
    
  end subroutine compute_field

  
  
!  subroutine write_field_on_grid(pepc_comm, time_pars, step, field_grid)
  subroutine write_field_on_grid(pepc_comm, step, field_grid)    
    use mpi
    use encap
    implicit none

    type(pepc_comm_t), intent(in)  :: pepc_comm
!    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in)            :: step
!    type(physics_pars_t), intent(in) :: physics_pars
    type(field_grid_t), intent(in) :: field_grid

    integer :: fh, mpi_err
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    integer(kind_particle) :: my_count, fflatshape(1)
    integer(kind = MPI_OFFSET_KIND) :: my_offset
    real(kind_physics), dimension(:), allocatable :: fflat
    real(kind = 8), dimension(:), allocatable :: real8_buf

    allocate(fflat(field_grid%ntot))
    fflatshape(1) = field_grid%ntot
    allocate(real8_buf(field_grid%nl))

    my_count = field_grid%nl
    my_offset = min(int(my_rank, kind=kind_particle), &
      mod(field_grid%ntot, int(n_ranks, kind=kind_particle))) &
      * (field_grid%ntot / n_ranks + 1) + &
      max(0_kind_particle, my_rank - mod(field_grid%ntot, int(n_ranks, kind=kind_particle))) * &
      (field_grid%ntot / n_ranks)
!
    call write_quantity_on_grid("potential", field_grid%p(:)%results%pot)
!    call write_quantity_on_grid("ex", field_grid%p(:)%results%e(1))
!    call write_quantity_on_grid("ey", field_grid%p(:)%results%e(2))
!    call write_quantity_on_grid("ez", field_grid%p(:)%results%e(3))
    call write_quantity_on_grid("Jx", field_grid%p(:)%results%J(1))
    call write_quantity_on_grid("Jy", field_grid%p(:)%results%J(2))
    call write_quantity_on_grid("Jz", field_grid%p(:)%results%J(3))
    call write_quantity_on_grid("Bx", field_grid%p(:)%results%B(1))
    call write_quantity_on_grid("By", field_grid%p(:)%results%B(2))
    call write_quantity_on_grid("Bz", field_grid%p(:)%results%B(3))
      
    call flatten_and_write_quantity_on_grid("ne", field_grid%ne)
    call flatten_and_write_quantity_on_grid("ni", field_grid%ni)
    call flatten_and_write_quantity_on_grid("nb", field_grid%nb)
    call flatten_and_write_quantity_on_grid("Ia", field_grid%Ia)
    call flatten_and_write_quantity_on_grid("Ip", field_grid%Ip)
!    call flatten_and_write_quantity_on_grid("qe", field_grid%qe)
!    call flatten_and_write_quantity_on_grid("qi", field_grid%qi)
!    call flatten_and_write_quantity_on_grid("qb", field_grid%qb)
!    call flatten_and_write_quantity_on_grid("vex", field_grid%vex)
!    call flatten_and_write_quantity_on_grid("vey", field_grid%vey)
!    call flatten_and_write_quantity_on_grid("vez", field_grid%vez)
!    call flatten_and_write_quantity_on_grid("vix", field_grid%vix)
!    call flatten_and_write_quantity_on_grid("viy", field_grid%viy)
!    call flatten_and_write_quantity_on_grid("viz", field_grid%viz)
!    call flatten_and_write_quantity_on_grid("vbx", field_grid%vbx)
!    call flatten_and_write_quantity_on_grid("vby", field_grid%vby)
!    call flatten_and_write_quantity_on_grid("vbz", field_grid%vbz)
    deallocate(fflat, real8_buf)

    contains

    subroutine flatten_and_write_quantity_on_grid(yname, y)
      implicit none

      character(*), intent(in) :: yname
      real(kind_physics), dimension(:,:,:), intent(in) :: y

      fflat(:) = reshape(y, fflatshape)
      call write_quantity_on_grid(yname, fflat(1 + my_offset : 1 + my_offset + field_grid%nl - 1))
    end subroutine

    subroutine write_quantity_on_grid(yname, y)
      use mpi
      use module_debug
      use module_globals, only:folder
      implicit none

      character(*), intent(in) :: yname
      real(kind_physics), dimension(:), intent(in) :: y

      integer(kind = MPI_OFFSET_KIND), parameter :: field_header_size = 512
      integer(kind = 4), parameter :: field_header_magic = 1

      character(1024) :: filename
      integer :: nx_, ny_, nz_ 
      integer(kind = MPI_OFFSET_KIND) :: mpi_disp

      real8_buf(:) = real(y(:), kind = 8)

      write(filename, '(a,"/",a,"_",i6.6,".bin")') trim(folder)//trim(field_dir), yname, step

      call mpi_file_open(MPI_COMM_WORLD, filename, &
        ior(MPI_MODE_RDWR, MPI_MODE_CREATE), MPI_INFO_NULL, fh, mpi_err)
        
      if (mpi_err .ne. MPI_SUCCESS) then
        DEBUG_ERROR(*, 'write_field_on_grid(): I/O error for ', filename)
      end if

      nx_ = int(field_grid%n(1))
      ny_ = int(field_grid%n(2))
      nz_ = int(field_grid%n(3))

      call mpi_file_set_view(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, &
        'native', MPI_INFO_NULL, mpi_err)

      ! write file header on rank 0
      ! int32   magic
      ! int32   nx
      ! int32   ny
      ! float64 offset_x
      ! float64 offset_y
      ! float64 extent_x
      ! float64 extent_y
      ! int32   step
      ! float64 t
      ! float64 B0
      ! float64 vte
      ! float64 vti
      ! float64 qe
      ! float64 qi
      ! float64 me
      ! float64 mi
      ! float64 min
      ! float64 max
        
      if (pepc_comm%mpi_rank == 0) then
        call header_write_integer4(field_header_magic)
        call header_write_integer4(nx_)
        call header_write_integer4(ny_)
        call header_write_integer4(nz_)
        call header_write_real(field_grid%offset(1))
        call header_write_real(field_grid%offset(2))
        call header_write_real(field_grid%offset(3))
        call header_write_real(field_grid%extent(1))
        call header_write_real(field_grid%extent(2))
        call header_write_real(field_grid%extent(3))
        call header_write_integer4(step)
!        call header_write_real(step * time_pars%dt)
!        call header_write_real(physics_pars%vte)
!        call header_write_real(physics_pars%vti)
!        call header_write_real(physics_pars%qe)
!        call header_write_real(physics_pars%qi)
!        call header_write_real(physics_pars%me)
!        call header_write_real(physics_pars%mi)
        
        call header_write_real(minval(y))
        call header_write_real(maxval(y))

        call mpi_file_get_position(fh, mpi_disp, mpi_err)
        if (mpi_disp > field_header_size) then
          DEBUG_ERROR(*, "write_field_on_grid(): field_header_size is too small: ", field_header_size, " < ", mpi_disp)
        end if
      end if

      call mpi_file_set_view(fh, field_header_size, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, mpi_err)
      call mpi_file_write_at_all(fh, my_offset, real8_buf, int(field_grid%nl, kind = kind_default), MPI_REAL8, mpi_stat, mpi_err)

      call mpi_file_sync(fh, mpi_err)
      call mpi_file_close(fh, mpi_err)
      
    end subroutine

    subroutine header_write_real(x)
      implicit none

      real(kind_physics), intent(in) :: x

      real(kind = 8) :: real8_buf

      real8_buf = real(x, kind = 8)
      call mpi_file_write(fh, real8_buf, 1, MPI_REAL8, mpi_stat, mpi_err)
    end subroutine

    subroutine header_write_integer4(i)
      implicit none

      integer(kind = 4), intent(in) :: i

      call mpi_file_write(fh, i, 1, MPI_INTEGER4, mpi_stat, mpi_err)
    end subroutine

  end subroutine
  
  
    subroutine write_field_on_grid3D(field_grid, step,realtime)
    use module_interaction_specific, only : force_law
    use module_vtk_helpers
    use module_pepc_types
    use mpi
    use module_shortcut            , only: three, half
    use encap
    implicit none

    type(field_grid_t)             , intent(in)    :: field_grid
    integer                        , intent(in)    :: step
    real(kind_particle)            , intent(in)    :: realtime
    integer(kind_particle)                         :: ip,rc=15
    integer                       , dimension(2,3) :: dims
    real(kind_particle)           , allocatable    :: x(:),y(:),z(:)
    
    
    if (allocated(x) ) deallocate(x)
    allocate(x(field_grid%n(1)+1), stat=rc)
    if (allocated(y) ) deallocate(y)
    allocate(y(field_grid%n(2)+1), stat=rc)
    if (allocated(z) ) deallocate(z)
    allocate(z(field_grid%n(3)+1), stat=rc)
    
    do ip = 1,field_grid%n(1)+1
        x(ip)  = (ip - three*half) * field_grid%dx(1) + field_grid%offset(1)
    enddo
    
    do ip = 1,field_grid%n(2)+1
        y(ip)  = (ip - three*half) * field_grid%dx(2) + field_grid%offset(2)
    enddo
    
    do ip = 1,field_grid%n(3)+1
        z(ip)  = (ip - three*half) * field_grid%dx(3) + field_grid%offset(3)
    enddo
    
    dims(:, :)           = 0
    dims(2, 1:force_law) = field_grid%n(1:force_law) 
    
!    call vtk_write_scalar_on_grid("Electron", step, realtime, 0, dims, dims, &
!        x, y, z, field_grid%ne,  'n_el' , MPI_COMM_SELF )

    call vtk_write_densities_on_grid("Electron", step, realtime, 0, dims, dims, &
        x, y, z, field_grid%ne,  'n_el' ,field_grid%ni,  'n_io' , MPI_COMM_SELF )        
    
        
        
        contains


        subroutine vtk_write_scalar_on_grid(filename, step, tsim, vtk_step, globaldims, mydims, xcoords, ycoords, zcoords, &
                          dens1, name1, mpi_comm, coord_scale)
          use module_vtk
          implicit none
          include 'mpif.h'
          character(*), intent(in) :: filename, name1!, name2
          integer, intent(in) :: step
          integer, intent(in) :: vtk_step
          real*8, intent(in) :: tsim
          integer, dimension(2,3), intent(in) :: globaldims, mydims
          real*8, intent(in) :: xcoords(:), ycoords(:), zcoords(:)
          real*8, intent(in) :: dens1(:,:,:)
          integer, intent(in) :: mpi_comm
          real*8, intent(in), optional :: coord_scale

          integer :: mpi_rank, mpi_size, ierr

          type(vtkfile_rectilinear_grid) :: vtk

          call MPI_Comm_rank(mpi_comm, mpi_rank, ierr)
          call MPI_Comm_size(mpi_comm, mpi_size, ierr)

          call vtk%create_parallel(trim(filename), step, mpi_rank, mpi_size, tsim, vtk_step)
            call vtk%set_communicator(mpi_comm)
            call vtk%write_headers(globaldims, mydims)
              call vtk%startcoordinates()
                call vtk%write_data_array("x_coordinate", xcoords, coord_scale)
                call vtk%write_data_array("y_coordinate", ycoords, coord_scale)
                call vtk%write_data_array("z_coordinate", zcoords, coord_scale)
              call vtk%finishcoordinates()
              call vtk%startpointdata()
                ! no point data here
              call vtk%finishpointdata()
              call vtk%startcelldata()
                call vtk%write_data_array(name1, dens1)
      !          call vtk%write_data_array(name2, dens2)
              call vtk%finishcelldata()
            call vtk%write_final()
          call vtk%close()
        end subroutine
        
    end subroutine write_field_on_grid3D
    
  

end module field_helper
