!  ================================
!
!         SETUP
!
!   $Revision$
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine pepc_setup()
  use physvars
  use tree_utils
  use module_fmm_framework
  use tree_walk_pthreads
  implicit none

  !  Default input set
  db_level        =   0
  np_mult         = -45
  weighted        =   1

  ! physics stuff
  force_const = 1.
  mac         = 0
  theta       = 0.6
  eps           = 0.01

  if (n_cpu.eq.1) then
     nppm=int(1.5*npart_total + 1000)  ! allow for additional ghost particles for field plots
!  else if (np_mult<0) then
!     nppm = abs(np_mult)*max(npart_total/n_cpu,1000) ! allow 50% fluctuation
  else
     nppm = int(1.5*max(npart_total/n_cpu,1000)) ! allow 50% fluctuation
  end if


  allocate ( pelabel(nppm), work(nppm) )

end subroutine pepc_setup




