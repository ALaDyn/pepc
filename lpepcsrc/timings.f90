module timings

  ! global timings
  real*8 :: t_domains=0., t_allocate=0., t_build=0., t_branches=0., t_fill=0., t_properties=0., t_restore=0., &
            t_walk=0., t_walkc=0., t_force=0., t_deallocate=0., t_all=0., t_local=0., t_exchange=0., t_global=0.
  real*8 :: t0_domains=0., t0_allocate=0., t0_build=0., t0_branches=0., t0_fill=0., t0_properties=0., t0_restore=0., &
            t0_walk=0., t0_walkc=0., t0_force=0., t0_deallocate=0., t0_all=0., t0_local=0., t0_exchange=0., t0_global=0.


  ! fields internal
  real*8 :: t_fields_begin=0., t_fields_tree=0., t_fields_nshort=0., t_fields_passes=0., t_fields_stats=0., t_restore_async
  
  ! tree_domains
  real*8 :: t_domains_keys=0., t_domains_sort=0., t_domains_sort_pure=0., t_domains_ship=0.,t_domains_bound=0.

  ! tree_allocate
  real*8 :: t_allocate_async=0.

  ! tree_build
  real*8 :: t_build_neigh=0., t_build_part=0., t_build_byte=0.

  ! tree_branches
  real*8 :: t_branches_find=0., t_branches_exchange=0., t_branches_integrate=0.

  ! tree_fill
  real*8 :: t_fill_local=0., t_fill_global=0.

  ! tree_props
  real*8 :: t_props_leafs=0., t_props_twigs=0., t_props_branches=0., t_props_global=0.


end module timings
