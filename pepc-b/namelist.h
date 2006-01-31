


  namelist /pepcdata/ nep, nip, ne, ni, &
       theta, mac, mass_ratio, eps, force_tolerance, &
       plasma_config, target_geometry, velocity_config, ispecial, &
       Te_keV, Ti_keV, T_scale, Zion, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, force_const, &
       n_layer, x_layer, y_layer, z_layer, r_layer, Zion_layer, rho_layer, mratio_layer, &
       xl, yl, zl, displace, bond_const, fnn, rho_min, lolam, &
       beam_config_in, np_beam, np_error, idim, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, ivis_domains, iprot, itrack, ncpu_merge, ngx, ngy, ngz, &
       vis_on, steering,  mc_init, restart, scheme, particle_bcs, &
       coulomb,  bfields,  bonds, lenjones, target_dup, ramp, &
       debug_level, debug_tree, ncpu_merge, balance, ifreeze, &
       constrain_proof, len_tripod, struct_step, uthresh, domain_cut, &
       debug_rank, np_mult, fetch_mult,nbuf_max, te_perturb, tpert, kpert, &
       q_factor
