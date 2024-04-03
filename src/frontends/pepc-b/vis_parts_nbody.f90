! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

! ======================
!
!   VIS_PARTS_NBODY
!
!   Send particle data to VISIT for real-time visualisation
!
!
! ======================

subroutine vis_parts_nbody(vcount)

   use physvars
   use treevars
   use module_spacefilling
   use mpi
   implicit none

   ! Buffer arrays vbuffer, vbuf_local allocated in setup_arrays

   integer, dimension(num_pe) :: nparts_pe, recv_strides, nbuf_pe  ! array of npp on each PE
   integer, dimension(num_pe) :: istrid, ships, fetches  ! arrays of ships/fetches
   integer :: icolour(nbuf_max)
   integer :: lvisit_active, nskip, nproot, ne_buf, ni_buf, npart_buf, nbufe, nbufi
   integer :: i, j, k, p, ioffset, ixd, iyd, izd, ilev, lcount, wfdatai

   real :: s, simtime, dummy, xd, yd, zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz, phipond
   real :: plasma1, plasma2, plasma3, t_display, wfdatar, u2, amp_las, box_x, box_y, box_z
   integer :: nship, ierr, cpuid
   integer :: type
   integer :: vbufcols = 22, incdf, ndom_vis, ivisdom, ndomain_vis = 5000, npart_vis
   real :: lbox, work_ave, upmax, uproton_max, uxmax, boxz_min
   logical :: vis_debug = .false.
   logical :: pick = .false.
   integer :: vcount, iglobal, ibox, max_ships, max_fetches, min_ships, min_fetches
   integer :: nproc_vis, proc_vis(num_pe)
   real :: rlev, ave_ships, ave_fetches
   integer*8 :: bkey
   convert_mu = 1.
   simtime = dt * (itime + itime_start)

   pick_particles: select case (mod(vis_select, 10))
   case (1, 3)
      nproot = npart
   case (2)
      nproot = n_layer(1)
   case default
      nproot = 0
   end select pick_particles

   npart_vis = nbuf_max - ndomain_vis  ! max # particles sent to visualization
   if (nproot .gt. npart_vis) then
      nskip = nproot / npart_vis + 2
   else
      nskip = 1
   end if

   if (me .eq. 0 .and. nproot .gt. npart_vis) then
      write (*, '(a)') "VIS_PARTS   | # particles > vis nbuf_max - reducing number shipped"
      write (*, '(a23,i8,a10,i8)') "VIS_PARTS   | nbuf_max=", nbuf_max, " nskip=", nskip
   end if

   if (beam_config .eq. 4) then
      amp_las = vosc * min(1., simtime / tpulse)
   else
      amp_las = vosc
   end if

   !    if (.not. launch) then
   t_display = vcount
   !    else if (beam_config>=3 .and. beam_config<=6) then
   !       t_display = tlaser*convert_fs
   !    else
   !      t_display = simtime
   !    endif

   if (beam_config .le. 3) then
      plasma1 = Ukine * convert_keV
      plasma2 = Ukini * convert_keV
   else
      plasma1 = sigma
      plasma2 = tpulse
   end if

   if (nskip .gt. 1 .and. uthresh .eq. -1) uthresh = vte**2

   nbufe = 0
   nbufi = 0

   if (me .eq. 0) then
      nship = 1  ! leave gap on root for info block (filled in later)
   else
      nship = 0
   end if

   ! Filter out particles for shipping
   type = 1
   upmax = 0.
   uxmax = 0.
   uproton_max = 0.
   ! TODO:  separate interesting ions and electrons
   do i = 1, npp
      u2 = 0.5 * 0.511 * mratio_layer(1) * (ux(i)**2 + uy(i)**2 + uz(i)**2) ! in MeV
      pick_part: select case (mod(vis_select, 10))
      case (1)  ! select from bulk
         pick = (npart .lt. nbuf_max .or. (npart .gt. nbuf_max .and. mod(pelabel(i), nskip) .eq. 0))
      case (2)  ! select from 2nd layer
         pick = (mod(pelabel(i), nskip) .eq. 0 .and. pelabel(i) .gt. npart - n_layer(1))
      case (3)  ! select above energy threshold
         pick = (u2 .gt. uthresh)
      case (4)  ! select protons only
         pick = (mod(pelabel(i), nskip) .eq. 0 .and. q(i) .gt. 0 .and. m(i) .eq. mass_proton)
      case default
         pick = .false.
      end select pick_part

      if (pick) then
         nship = nship + 1
         if (q(i) .lt. 0) then
            nbufe = nbufe + 1
         else
            nbufi = nbufi + 1
         end if

         upmax = max(upmax, u2)
         uxmax = max(1.d0 * uxmax, ux(i))
         ! Store attributes for visualizing
         vbuf_local(0, nship) = t_display
         vbuf_local(1, nship) = dt
         ! scaled coordinates
         vbuf_local(2, nship) = convert_mu * x(i)
         vbuf_local(3, nship) = convert_mu * y(i)
         vbuf_local(4, nship) = convert_mu * z(i)
         ! velocities
         vbuf_local(5, nship) = ux(i)
         vbuf_local(6, nship) = uy(i)
         vbuf_local(7, nship) = uz(i)
         ! E-field
         vbuf_local(8, nship) = Ex(i)
         vbuf_local(9, nship) = Ey(i)
         vbuf_local(10, nship) = Ez(i)
         vbuf_local(11, nship) = 0.
         vbuf_local(12, nship) = 0.
         vbuf_local(13, nship) = 0.
         vbuf_local(14, nship) = pelabel(i)
         vbuf_local(15, nship) = 0.
         vbuf_local(16, nship) = type ! particle type
         vbuf_local(17, nship) = q(i)
         vbuf_local(18, nship) = m(i)
         vbuf_local(19, nship) = pepid(i)
         vbuf_local(20, nship) = 0.
         vbuf_local(21, nship) = 0.
      end if
   end do

   ! Find # particles to collect from each cpu
   call MPI_ALLGATHER(nship, 1, MPI_INTEGER, nparts_pe, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(upmax, uproton_max, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)

   ! Total buffer lengths/pe
   nbuf_pe = nparts_pe * attrib_max

   recv_strides(1:num_pe) = (/0, (SUM(nbuf_pe(1:i - 1)), i=2, num_pe)/)

   npart_buf = (SUM(nparts_pe))  ! global total # particles to visualize

   call MPI_ALLREDUCE(nbufe, ne_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(nbufi, ni_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_GATHER(sum_fetches, 1, MPI_INTEGER, fetches, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_GATHER(sum_ships, 1, MPI_INTEGER, ships, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (me .eq. 0) write (*, '(a28,i8,a10,i8)') "VIS_PARTS   | selected nbuf=", npart_buf, " out of ", npart_vis
   if (me .eq. 0) write (*, '(a28,i8)') "VIS_PARTS   | total # branches = ", nbranch_sum

   if (npart_buf .gt. 0 .and. npart_buf .lt. npart_vis) then
      ! send  particle data from all PEs
      ! increase momentum threshold if close to max ship #
      if (1.0 * npart_buf / npart_vis .gt. 0.9) uthresh = uthresh * 1.1

      call MPI_GATHERV(vbuf_local, nship * attrib_max, MPI_REAL, vbuffer, nbuf_pe, recv_strides, &
                       MPI_REAL, 0, MPI_COMM_WORLD, ierr)

      ! Ship to visualization
      if (me .eq. 0) then

         ! Stick info block in front of particle data on root

         type = 8

         vbuffer(0, 1) = t_display
         !              vbuffer(1,1) = npart_buf-1
         vbuffer(1, 1) = npart
         ! scaled coordinates
         vbuffer(2, 1) = t_display
         vbuffer(3, 1) = 0
         vbuffer(4, 1) = 1
         ! velocities
         vbuffer(5, 1) = npart_buf - 1
         vbuffer(6, 1) = npart
         vbuffer(7, 1) = ops_per_sec
         vbuffer(8, 1) = uproton_max
         vbuffer(9, 1) = 0.
         vbuffer(10, 1) = 0.
         vbuffer(11, 1) = 0.
         vbuffer(12, 1) = 0.
         vbuffer(13, 1) = 0.
         vbuffer(14, 1) = 0.
         vbuffer(15, 1) = 0.
         vbuffer(16, 1) = type ! particle type
         vbuffer(17, 1) = 0.
         vbuffer(18, 1) = 0.
         vbuffer(19, 1) = 0.
         vbuffer(20, 1) = 0.
         vbuffer(21, 1) = 0.

         work_ave = SUM(work_loads) / num_pe

         ! Add branch nodes to show domains
         ndom_vis = 0

         !              if (ndomain_vis>0) then

         pick_domain: select case (vis_select / 10)

         case (1)

            ! Show all branch boxes below z-cutoff

            ! First count how many to be shipped
            boxz_min = zl
            do k = 1, nbranch_sum

               ilev = level_from_key(1.*branch_key(k))
               ixd = SUM((/(2**i * ibits(branch_key(k), 3 * i, 1), i=0, ilev - 1)/))
               iyd = SUM((/(2**i * ibits(branch_key(k), 3 * i + 1, 1), i=0, ilev - 1)/))
               izd = SUM((/(2**i * ibits(branch_key(k), 3 * i + 2, 1), i=0, ilev - 1)/))
               lbox = boxsize / 2**(ilev)          !  box length
               box_x = xmin + lbox * (ixd + .5) ! box centres
               box_y = ymin + lbox * (iyd + .5) ! box centres
               box_z = zmin + lbox * (izd + .5) ! box centres
               boxz_min = min(boxz_min, box_z)
               if (box_z .lt. domain_cut) ndom_vis = ndom_vis + 1
            end do

            ivisdom = 0
            do k = 1, nbranch_sum

               ilev = level_from_key(branch_key(k))
               ixd = SUM((/(2**i * ibits(branch_key(k), 3 * i, 1), i=0, ilev - 1)/))
               iyd = SUM((/(2**i * ibits(branch_key(k), 3 * i + 1, 1), i=0, ilev - 1)/))
               izd = SUM((/(2**i * ibits(branch_key(k), 3 * i + 2, 1), i=0, ilev - 1)/))
               lbox = boxsize / 2**(ilev)          !  box length
               box_x = xmin + lbox * (ixd + .5) ! box centres
               box_y = ymin + lbox * (iyd + .5) ! box centres
               box_z = zmin + lbox * (izd + .5) ! box centres

               if (box_z .lt. domain_cut) then   ! only show domains below z=domain_cut
                  ivisdom = ivisdom + 1
                  j = npart_buf + 1 + ivisdom
                  ! Store attributes for visualizing
                  vbuffer(0, j) = t_display
                  vbuffer(1, j) = box_x
                  vbuffer(2, j) = box_y
                  vbuffer(3, j) = box_z
                  vbuffer(4, j) = lbox
                  vbuffer(5, j) = branch_owner(k)  ! cpu id of branch node
                  cpuid = branch_owner(k) + 1
                  vbuffer(6, j) = work_loads(cpuid) / work_ave   ! total work load of branch node/cpu
                  vbuffer(7, j) = npps(cpuid) ! total # particles on cpu
                  vbuffer(8, j) = 0.
                  vbuffer(9, j) = 0.
                  vbuffer(10, j) = 0.
                  vbuffer(11, j) = 0.
                  vbuffer(12, j) = 0.
                  vbuffer(13, j) = 0.
                  vbuffer(14, j) = 0.
                  vbuffer(15, j) = 0.
                  vbuffer(16, j) = 16    ! Domain type
                  vbuffer(17, j) = ndom_vis  ! Total # branch nodes
                  vbuffer(18, j) = 0.
                  vbuffer(19, j) = 0.
                  vbuffer(20, j) = 0.
                  vbuffer(21, j) = 0.
               end if
               !        write (*,'(7f13.4)') vbuffer(1,j),vbuffer(2,j), vbuffer(4,j), vbuffer(5,j), &
               !             vbuffer(6,j), vbuffer(7,j), vbuffer(17,j)
            end do

         case (2)

            ! Show first box in each processor domain

            ! precalculate strides
            write (*, '(a28)') "VIS_PARTS   | case 2"
            ivisdom = 0

            istrid(1:num_pe) = (/0, (SUM(nbranches(1:i - 1)), i=2, num_pe)/)
            ave_ships = SUM(ships(1:num_pe)) / num_pe
            ave_fetches = SUM(fetches(1:num_pe)) / num_pe
            max_ships = maxval(ships(1:num_pe))
            max_fetches = maxval(fetches(1:num_pe))
            min_ships = minval(ships(1:num_pe))
            min_fetches = minval(fetches(1:num_pe))
            write (*, '(a28,1pe12.2,a2,1pe12.2,a2,1pe12.2,a2,1pe12.2)') "VIS_PARTS   | ranges ships/fetches: ", &
               min_ships / ave_ships, "-", max_ships / ave_ships, ",", min_fetches / ave_fetches, "-", max_fetches / ave_fetches
            ! Pick biggest box in each domain
            ndom_vis = num_pe

            do k = 1, num_pe

               iglobal = istrid(k) + 1  ! branches are sorted, so 1st one will have lowest level #
               rlev = level_from_key(branch_key(iglobal))
               ilev = rlev
               ibox = iglobal
               bkey = branch_key(iglobal)

               !                    write(*,'(a28,3i10)') "VIS_PARTS   | branch ",i,iglobal, ilev

               ixd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i, 1), i=0, ilev - 1)/))
               iyd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i + 1, 1), i=0, ilev - 1)/))
               izd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i + 2, 1), i=0, ilev - 1)/))
               lbox = boxsize / 2**(ilev)          !  box length
               box_x = xmin + lbox * (ixd + .5) ! box centres
               box_y = ymin + lbox * (iyd + .5) ! box centres
               box_z = zmin + lbox * (izd + .5) ! box centres
               ivisdom = ivisdom + 1
               j = npart_buf + 1 + ivisdom
               ! Store attributes for visualizing
               vbuffer(0, j) = t_display
               vbuffer(1, j) = box_x
               vbuffer(2, j) = box_y
               vbuffer(3, j) = box_z
               vbuffer(4, j) = lbox
               vbuffer(5, j) = branch_owner(ibox)  ! cpu id of branch node
               cpuid = branch_owner(ibox) + 1
               vbuffer(6, j) = 1.*ships(cpuid) / ave_ships   ! # ships made by this cpu
               vbuffer(7, j) = 1.*fetches(cpuid) / ave_fetches ! # fetches made by this cpu
               write (*, '(a28,i8,a25,i8,o22,1pe12.3,i10,i10,1pe12.2,i10,1pe12.2)') &
                  "VIS_PARTS   | pe: ", k, "box, level, ships, fetches ", ibox, bkey, rlev, ilev, ships(cpuid), &
                  vbuffer(6, j), fetches(cpuid), vbuffer(7, j)
               vbuffer(8, j) = 0.
               vbuffer(9, j) = 0.
               vbuffer(10, j) = 0.
               vbuffer(11, j) = 0.
               vbuffer(12, j) = 0.
               vbuffer(13, j) = 0.
               vbuffer(14, j) = 0.
               vbuffer(15, j) = 0.
               vbuffer(16, j) = 16    ! Domain type
               vbuffer(17, j) = ndom_vis  ! Total # branch nodes
               vbuffer(18, j) = 0.
               vbuffer(19, j) = 0.
               vbuffer(20, j) = 0.
               vbuffer(21, j) = 0.
            end do

         case (3)

            ! Show domain boxes of selected processors

            ! precalculate strides
            write (*, '(a28)') "VIS_PARTS   | case 3"
            ivisdom = 0

            istrid(1:num_pe) = (/0, (SUM(nbranches(1:i - 1)), i=2, num_pe)/)
            ave_ships = SUM(ships(1:num_pe)) / num_pe
            ave_fetches = SUM(fetches(1:num_pe)) / num_pe
            max_ships = maxval(ships(1:num_pe))
            max_fetches = maxval(fetches(1:num_pe))
            min_ships = minval(ships(1:num_pe))
            min_fetches = minval(fetches(1:num_pe))
            write (*, '(a28,1pe12.2,a2,1pe12.2,a2,1pe12.2,a2,1pe12.2)') &
               "VIS_PARTS   | ranges ships/fetches: ", min_ships / ave_ships, "-", max_ships / ave_ships, ",", &
               min_fetches / ave_fetches, "-", max_fetches / ave_fetches
            ivisdom = 0
            ndom_vis = 0
            nproc_vis = 5
            proc_vis(1:nproc_vis) = (/5, 297, 563, 804, 1023/)  ! selected cpus for vis

            !  Determine total box #
            do p = 1, nproc_vis
               ndom_vis = ndom_vis + nbranches(proc_vis(p))
            end do

            do p = 1, nproc_vis

               do k = 1, nbranches(proc_vis(p))
                  iglobal = istrid(proc_vis(p)) + k  ! global branch index
                  rlev = level_from_key(branch_key(iglobal))
                  ilev = rlev
                  ibox = iglobal
                  bkey = branch_key(iglobal)

                  !                    write(*,'(a28,3i10)') "VIS_PARTS   | branch ",i,iglobal, ilev

                  ixd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i, 1), i=0, ilev - 1)/))
                  iyd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i + 1, 1), i=0, ilev - 1)/))
                  izd = SUM((/(2**i * ibits(branch_key(ibox), 3 * i + 2, 1), i=0, ilev - 1)/))
                  lbox = boxsize / 2**(ilev)          !  box length
                  box_x = xmin + lbox * (ixd + .5) ! box centres
                  box_y = ymin + lbox * (iyd + .5) ! box centres
                  box_z = zmin + lbox * (izd + .5) ! box centres
                  ivisdom = ivisdom + 1
                  j = npart_buf + 1 + ivisdom
                  ! Store attributes for visualizing
                  vbuffer(0, j) = t_display
                  vbuffer(1, j) = box_x
                  vbuffer(2, j) = box_y
                  vbuffer(3, j) = box_z
                  vbuffer(4, j) = lbox
                  vbuffer(5, j) = branch_owner(ibox)  ! cpu id of branch node
                  cpuid = branch_owner(ibox) + 1
                  vbuffer(6, j) = 1.*ships(cpuid) / ave_ships   ! # ships made by this cpu
                  vbuffer(7, j) = 1.*fetches(cpuid) / ave_fetches ! # fetches made by this cpu
                  write (*, '(a28,2i8,a25,i8,o22,1pe12.3,i10,i10,1pe12.2,i10,1pe12.2)') &
                     "VIS_PARTS   | pe: ", p, k, "box, level, ships, fetches ", ibox, bkey, rlev, ilev, &
                     ships(cpuid), vbuffer(6, j), fetches(cpuid), vbuffer(7, j)
                  vbuffer(8, j) = 0.
                  vbuffer(9, j) = 0.
                  vbuffer(10, j) = 0.
                  vbuffer(11, j) = 0.
                  vbuffer(12, j) = 0.
                  vbuffer(13, j) = 0.
                  vbuffer(14, j) = 0.
                  vbuffer(15, j) = 0.
                  vbuffer(16, j) = 16    ! Domain type
                  vbuffer(17, j) = ndom_vis  ! Total # branch nodes
                  vbuffer(18, j) = 0.
                  vbuffer(19, j) = 0.
                  vbuffer(20, j) = 0.
                  vbuffer(21, j) = 0.
               end do
            end do

         case default
            ! skip domain boxes
            ndom_vis = 0
         end select pick_domain

         ! Fill out dummy values for netcdf
         do j = npart_buf + 1 + ndom_vis, nbuf_max
            vbuffer(0:attrib_max - 1, j) = 0.
         end do

         if (me .eq. 0) then
            write (*, '(2a,f12.1)') 'VIS_PARTS   |', ' Display time/count ', t_display
            write (*, '(2a,2i8)') 'VIS_PARTS   |', ' # particles shipped ', npart_buf, nship
            write (*, '(2a,i8,a2,i8)') 'VIS_PARTS   |', ' # branches shipped ', ndom_vis, '/', nbranch_sum
            write (*, '(2a,f12.1,a2,f12.1)') 'VIS_PARTS   |', ' box min ', boxz_min, 'cutoff', domain_cut
            write (*, '(2a,i8,a2,i8)') 'VIS_PARTS   |', ' Total # objects shipped :', ndom_vis + 1 + npart_buf, ' /', nbuf_max
            write (*, '(2a,f12.2)') 'VIS_PARTS   |', ' u_thresh: (MeV)     ', uthresh
            write (*, '(2a,3f12.2)') 'VIS_PARTS   |', ' ux, up_max: (MeV)     ', uxmax, upmax, uproton_max
            write (*, '(2a,i8)') 'VIS_PARTS   |', ' vis_select=  ', vis_select
            if (vis_debug) then
               write (*, '(2a)') 'VIS_PARTS   |', ' t,1,x,y,z,q,label,owner,type :'
               do j = 1, npart_buf
                  write (*, '(i6,6f12.4,2i5,3i6)') j, vbuffer(0, j), vbuffer(1, j), vbuffer(2, j), vbuffer(3, j), vbuffer(4, j), &
                     vbuffer(17, j), int(vbuffer(14, j)), int(vbuffer(19, j)), int(vbuffer(16, j))
               end do
            end if
         end if

         ! ---- Preprocess VISIT setup -----------
#ifdef VISIT_NBODY

         call flvisit_nbody2_check_connection(lvisit_active)
         ! send particles and branch boxes together
         if (lvisit_active .ne. 0) then
            call flvisit_nbody2_partstep_send(vbuffer, npart_buf + ndom_vis + 1, attrib_max)
         end if
#else
         !  --- No VISIT installed ----
#endif
         !  --- end preprocess --------

         ! netcdf needs fixed buffer size, so take max used for initialisation
#ifdef NETCDFLIB
         if (netcdf) then
            call ncnbody_put(ncid, vbuffer, nbuf_max, attrib_max, incdf)
            write (*, '(a)') "VIS_PARTS   | Writing particles to netcdf"
         end if
#else
         !  --- No NETCDF installed ----
#endif

         !
         !        write (90,*) 'local',i
         !        write (90,'((22(f12.5/)//))') vbuf_local(0:attrib_max-1,i)
         !        write (90,*) 'global',i
         !        write (90,'((22(f15.8/)//))') vbuffer(0:attrib_max-1,i)
         !        write (90,*) vbuffer(0:attrib_max-1,i)
         !              end do

      end if

   else
      ! Just send particles on root

      ! TODO: this needs removing or fixing

      if (me .eq. 0) then
         if (npart_buf .gt. nbuf_max - 1000) then
            write (*, '(a)') 'VIS_PARTS   | Too many particles to ship: npart_buf= ', npart_buf, '/', nbuf_max
            uthresh = uthresh * 2
            write (*, *) 'VIS_PARTS   | Increasing momentum threshold to: ', sqrt(abs(uthresh))
         else if (npart_buf .eq. 0.and .ne. gt.0) then
            uthresh = vte * 2
            write (*, *) 'VIS_PARTS   | Reducing momentum threshold to: ', sqrt(abs(uthresh))
         end if
         nproot = 0.8 * npart / num_pe ! fixed # parts close to npp

         !           call flvisit_spk_check_connection(lvisit_active)
         !              call flvisit_nbody2_check_connection(lvisit_active)
         !        call flvisit_spk_info_send(npp,xl,yl,zl,zl,nep,nip,np_beam,itime+itime_start)
         wfdatai = int(tlaser)
         wfdatar = tlaser
         !           call flvisit_spk_info_send(nproot,xl,yl,zl, t_display, &
         !                x_crit, amp_las, sigma, tpulse, &
         !                ne,ni,npart,wfdatai)
         !           call flvisit_nbody2_info_send(t_display,xl,yl,zl, xl, &
         !             t_display, amp_las, plasma1, plasma2, xl, &
         !             xl, xl, xl, xl, xl, &
         !             xl, xl, xl, xl, xl )
         !           call flvisit_spk_particles_send(t_display,x,y,z,ux,uy,uz,q,pepid,pelabel,npp)
         !        call flvisit_nbody2_particles_send(t_display,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart_buf)
      end if
   end if

   ! Make sure everyone else knows about new momentum threshold
   call MPI_BCAST(uthresh, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

   vcount = vcount + 1

end subroutine vis_parts_nbody

