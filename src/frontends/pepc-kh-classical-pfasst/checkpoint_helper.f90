module checkpoint_helper
  implicit none

  contains

  subroutine write_checkpoint(pepc_pars, time_pars, step, physics_pars, field_grid, p)
    use module_pepc_types
    use module_checkpoint
    use encap
    use pepc_helper
    use time_helper
    use physics_helper
    use field_helper
    implicit none

    type(pepc_pars_t), intent(in) :: pepc_pars
    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in) :: step
    type(physics_pars_t), intent(in) :: physics_pars
    type(field_grid_t), intent(in) :: field_grid
    type(t_particle), dimension(:), intent(in) :: p

    character(len = 255) :: file_name

    call write_particles_mpiio(pepc_pars%pepc_comm%comm_space, &
      step, pepc_pars%np, p, file_name)

    if (pepc_pars%pepc_comm%rank_space == 0) then
      call write_field_grid_params(field_grid, file_name)
      call write_physics_params(physics_pars, file_name)
      call write_time_params(time_pars, step, file_name)
      call write_params(pepc_pars, file_name)
    end if

  end subroutine write_checkpoint

end module checkpoint_helper
