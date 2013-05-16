module module_interaction_partners
  use module_pepc_types
  use variables
  use module_interaction_specific
  use module_treediags
  use module_pepc

  implicit none

  type(t_particle), allocatable :: probe_particles(:)

  contains
    subroutine init_probes(num_of_probes)
      implicit none

      integer :: num_of_probes

      integer :: i

      if (root) then
          allocate(probe_particles(num_of_probes))
          allocate(interaction_keylist(num_of_probes,1000000))
          allocate(no_interaction_partners(num_of_probes))
          allocate(interaction_vbox(num_of_probes,1000000,3))
      else
          allocate(probe_particles(1))
          allocate(interaction_keylist(1,1000000))
          allocate(no_interaction_partners(1))
          allocate(interaction_vbox(1,1000000,3))
      end if
      no_interaction_partners=0

      if (root) then
        do i=1,num_of_probes
            probe_particles(i)%label       = i
            probe_particles(i)%data%q      = 1.0_8
            probe_particles(i)%data%m      = 1.0_8

            probe_particles(i)%results%e   = 0.0_8
            probe_particles(i)%results%pot = 0.0_8
            probe_particles(i)%work        = 1.0_8
            probe_particles(i)%data%species= 12

            probe_particles(i)%data%v(1:3)     =0.0_8
            probe_particles(i)%x(1) =(i-1)*dx/(num_of_probes-1.)
            probe_particles(i)%x(2) =(i-1)*dy/(num_of_probes-1.)
            probe_particles(i)%x(3) =(i-1)*dz/(num_of_probes-1.)

            probe_particles(i)%data%B(1)=Bx
            probe_particles(i)%data%B(2)=By
            probe_particles(i)%data%B(3)=Bz
        end do
      end if
    end subroutine init_probes


    subroutine get_interaction_partners(num_of_probes)
      use module_interaction_specific, only: interaction_nodelist, &
        no_interaction_partners, interaction_vbox
      implicit none

      integer :: num_of_probes

      integer :: i,help_num_of_probes

      if (root) then
          help_num_of_probes=num_of_probes
      else
          help_num_of_probes=0
      end if
      force_law=6
      call pepc_traverse_tree(probe_particles(1:help_num_of_probes))
      force_law=3
      do i=1,help_num_of_probes
          call write_interaction_partners_to_vtk(step, i,0.0_8, -1, interaction_nodelist, no_interaction_partners, interaction_vbox)
          write(*,*)i,no_interaction_partners(i)
      end do
      no_interaction_partners=0
    end subroutine get_interaction_partners

end module module_interaction_partners
