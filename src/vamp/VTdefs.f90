
SUBROUTINE VTdefs
  INCLUDE "VT.inc"
  INCLUDE 'VTcommon.h'
  INTEGER IE
  INTEGER IWALKC
  integer ibuildc
  integer initc
  integer iforcec
  integer idiagc
  integer iconc

  VTNOSCL=VT_NOSCL
  call VTCLASSDEF('F',ICLASSH,IE)
  call VTCLASSDEF('WALK',IWALKC,IE)
  call VTCLASSDEF('BUILD',IBUILDC,IE)
  call VTCLASSDEF('DIAG',IDIAGC,IE)
  call VTCLASSDEF('FOR',IFORCEC,IE)
  call VTCLASSDEF('CTRL',ICONC,IE)
  call VTCLASSDEF('INIT',INITC,IE)



  call VTFUNCDEF('treemp',ICONC,IF_treemp,IE)
  call VTFUNCDEF('bal',ICONC,IF_forces_bal,IE)
  call VTFUNCDEF('beam',ICONC,IF_beam,IE)
  call VTFUNCDEF('beam_control',ICONC,IF_beam_control,IE)

  call VTFUNCDEF('setup',INITC,IF_setup,IE)
  call VTFUNCDEF('configure',INITC,IF_configure,IE)
  call VTFUNCDEF('scramble_v',INITC,IF_scramble_v,IE)
  call VTFUNCDEF('predef',INITC,IF_predef_parts,IE)
  call VTFUNCDEF('mc_config',INITC,IF_mc_config,IE)
  call VTFUNCDEF('cold_start',INITC,IF_cold_start,IE)
  call VTFUNCDEF('randion',INITC,IF_randion,IE)
  call VTFUNCDEF('reset_ions',INITC,IF_reset_ions,IE)
  call VTFUNCDEF('openfiles',INITC,IF_openfiles,IE)

  call VTFUNCDEF('potenergy',IFORCEC,IF_potenergy,IE)
  call VTFUNCDEF('sumpot',IFORCEC,IF_sumpot,IE)
  call VTFUNCDEF('sumf',IFORCEC,IF_sum_force,IE)
  call VTFUNCDEF('pond',IFORCEC,IF_fpond,IE)
  call VTFUNCDEF('vel',IFORCEC,IF_velocities,IE)
  call VTFUNCDEF('push',IFORCEC,IF_push,IE)
  call VTFUNCDEF('unbal',IFORCEC,IF_forces,IE)
  call VTFUNCDEF('sum_lj',IFORCEC,IF_sum_lennardjones,IE)

  call VTFUNCDEF('domains',IBUILDC,IF_make_domains,IE)
  call VTFUNCDEF('build',IBUILDC,IF_tree_build,IE)
  call VTFUNCDEF('branches',IBUILDC,IF_make_branches,IE)
  call VTFUNCDEF('fill',IBUILDC,IF_tree_fill,IE)
  call VTFUNCDEF('check',IBUILDC,IF_check_table,IE)
  call VTFUNCDEF('props',IBUILDC,IF_tree_properties,IE)
  call VTFUNCDEF('prefetch',IBUILDC,IF_tree_prefetch,IE)

  call VTFUNCDEF('walk_local',IWALKC,IF_tree_walk,IE)
  call VTFUNCDEF('walk_fetch',IWALKC,IF_tree_nlswap,IE)


  call VTFUNCDEF('diags',IDIAGC,IF_diagnostics,IE)
  call VTFUNCDEF('kine',IDIAGC,IF_kinenergy,IE)
  call VTFUNCDEF('domains',IDIAGC,IF_draw_domains,IE)
  call VTFUNCDEF('tree2d',IDIAGC,IF_draw_tree2d,IE)
  call VTFUNCDEF('econs',IDIAGC,IF_energy_cons,IE)
  call VTFUNCDEF('draw2d',IDIAGC,IF_draw2d,IE)
  call VTFUNCDEF('diagtree',IDIAGC,IF_diagnose_tree,IE)
  call VTFUNCDEF('visitd',IDIAGC,IF_visit_dump,IE)
  call VTFUNCDEF('hash',IDIAGC,IF_draw2d_hash,IE)
  call VTFUNCDEF('lists',IDIAGC,IF_draw_lists,IE)
  call VTFUNCDEF('dump',IDIAGC,IF_dump,IE)
  call VTFUNCDEF('slices',IDIAGC,IF_slices,IE)
  call VTFUNCDEF('vis_f',IDIAGC,IF_vis_fields,IE)
  call VTFUNCDEF('vis_p',ICLASSH,IF_vis_parts,IE)

  call VTFUNCDEF('F_psrsperm_i8',ICLASSH,IF_psrsperm_i8,IE)
  call VTFUNCDEF('F_phase',ICLASSH,IF_phase,IE)
  call VTFUNCDEF('F_key2addr',ICLASSH,IF_key2addr,IE)
  call VTFUNCDEF('F_psrsperm_i4',ICLASSH,IF_psrsperm_i4,IE)
  call VTFUNCDEF('F_blank6',ICLASSH,IF_blank6,IE)
  call VTFUNCDEF('F_genran',ICLASSH,IF_genran,IE)
  call VTFUNCDEF('F_maxwell1',ICLASSH,IF_maxwell1,IE)
  call VTFUNCDEF('F_make_hashentry',ICLASSH,IF_make_hashentry,IE)
  call VTFUNCDEF('F_sort_i',ICLASSH,IF_sort_i,IE)
  call VTFUNCDEF('F_cput',ICLASSH,IF_cput,IE)
  call VTFUNCDEF('F_key2node',ICLASSH,IF_key2node,IE)
  call VTFUNCDEF('F_nwaymrg',ICLASSH,IF_nwaymrg,IE)
  call VTFUNCDEF('F_driver',ICLASSH,IF_driver,IE)
  call VTFUNCDEF('F_indsort_i',ICLASSH,IF_indsort_i,IE)
  call VTFUNCDEF('F_psrsperm_r8',ICLASSH,IF_psrsperm_r8,IE)

  call VTFUNCDEF('F_keytest',ICLASSH,IF_keytest,IE)
  call VTFUNCDEF('F_blankn',ICLASSH,IF_blankn,IE)
  call VTFUNCDEF('F_swap_ab',ICLASSH,IF_swap_ab,IE)
  call VTFUNCDEF('F_sift_down',ICLASSH,IF_sift_down,IE)
  call VTFUNCDEF('F_rano',ICLASSH,IF_rano,IE)
  call VTFUNCDEF('F_closefiles',ICLASSH,IF_closefiles,IE)
  call VTFUNCDEF('F_constrain',ICLASSH,IF_constrain,IE)
  call VTFUNCDEF('F_pswssort',ICLASSH,IF_pswssort,IE)
  call VTFUNCDEF('F_psrssort',ICLASSH,IF_psrssort,IE)




      END
