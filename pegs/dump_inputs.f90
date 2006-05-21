!  ================================
!
!         DUMP_INPUTS
!
!   $Revision:  $
!
!   Echo inputs to run protocol
!     
!
!  ================================


subroutine dump_inputs
  use treevars
  use physvars
  use utils

  implicit none
  include 'mpif.h'
  integer :: ifile, machinebits
  character*7 :: ensembles(1:5)=(/ &
       'const U','Te, Ti ','glob Te','loc  Te','Ti only' /)


  if (my_rank==0) then
     do ifile = 6,15,9

        write (ifile,*) ' Scheme: ',ensembles(scheme)
        write (ifile,'(a20,1pe12.3)') ' mass of disc particles: ',qe
        write (ifile,'(a20,1pe12.3)') ' Max timestep: ',0.45*sqrt(3.)*eps**2/qe
        write (ifile,'(a20,1pe12.3)') ' Neighbour search radius: ',r_neighbour


        write (ifile,*) ' no.gas/dust: ', ndisc
        write (ifile,*) ' no.stars: ', nstar

        write (ifile,*) ' Particles per PE: ', npp


        write (ifile,'(/a/a)') ' Switches:','--------'
        write (ifile,'(a20,l3)') ' Gravitational forces: ',coulomb

        write (ifile,'(a20,l3)') ' load balance: ',load_balance
        write (ifile,'(a20,l3)') ' walk balance: ',walk_balance
        write (ifile,'(a20,l3)') ' restart: ',restart
        write (ifile,'(a20,i3)') ' tree_debug: ',domain_debug

        write (ifile,'(a20,l3)') ' domain debug: ',domain_debug
        write (ifile,'(a20,l3)') ' walk debug: ',walk_debug
        write (ifile,'(a20,l3)') ' dump tree: ',dump_tree
        write (ifile,'(a20,l3)') ' performance anal.: ',perf_anal
        write (ifile,'(a20,l3)') ' visit: ',vis_on
        write (ifile,'(a20,l3/)') ' steering: ',steering
        write (ifile,*) 'Max address in #-table: ',2**nbaddr-1
        machinebits = bit_size(1_8)    ! # bits in integer variable (hardware) 
        write (ifile,*) 'Machine bit-size = ',machinebits
        write (ifile,*) 'Max permitted particles / PE: ', nppm
        write (ifile,*) 'Max size of interaction lists: ', nintmax
	write (ifile,*) 'Shortlist length: ',nshortm
        write (ifile,*) 'Memory needed for lists = ',4*nintmax*nshortm*8/2**20,' Mbytes'



     end do
  endif



end subroutine dump_inputs
