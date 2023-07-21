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
!>
!> This is a little test to check for consistency while refactoring the interaction kernels.
!> Input/Output to routines added by hand in the hope to reflect the actual math.
!>
#ifdef TEST_INTERACTION
program test_interaction
   use module_interaction_specific
   use module_pepc_types
   use treevars
   implicit none

   type(t_tree_node_interaction_data) :: node
   type(t_particle) :: particle
   integer(kind_node) :: node_idx
   real(kind_physics) :: delta(3), dist2
   real(kind_physics) :: vbox(3)

   write (*, *) "Testing interaction kernels"

   ! Set up interaction example
   ! All numbers 'monkey random'
   particle = t_particle( &                                                                  !&
                 [-0.75777575896313332, -0.92130003029573582, -0.40461090836553604], &       !&
                 5.0000000000000000_8, &                                                     !&
                 1157908711597454372_kind_key, &                                             !&
                 5_kind_node, &                                                              !&
                 11507_kind_particle, &                                                      !&
                 t_particle_data( &                                                          !&
                    [7.0759621470156551E-005_kind_physics, -5.8397610295706905E-005_kind_physics, 0._kind_physics], &  !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics], &                   !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics], &                   !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics], &                   !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics] &                    !&
                 ), &                                                                        !&
                 t_particle_results( &                                                       !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics], &                   !&
                    [0._kind_physics, 0._kind_physics, 0._kind_physics], &                   !&
                    0._kind_physics &                                                        !&
                 ) &                                                                         !&
                 )
   node = t_tree_node_interaction_data( &
          [-1.1002075655195405, -0.72729974965783273, 0.34083351528516737], &
          9.1745326414205950E-002, 5.0830903238316002E-001, &
          -7.6376856408379043E-003, 6.473425384543E-001, &
          9.4745346463262350E-002, 2.2524523238316002E-005, &
          4.1545323245205950E-004, 1.4830903463636002E-004, &
          2.6745323463675950E-002, 3.1411123526409867E-002, &
          9.1445926252527474E-003, 8.0830363668316002E-002, &
          5.7785326414436950E-005, 7.0836482343523652E-003, &
          2.3745325254745950E-006, 6.0830903234567402E-002, &
          7.3745325254203630E-008, 2.0830245245216002E-005, &
          6.1342525414236340E-002, 7.0830903245234522E-003, &
          8.2545326414236363E-002, 3.0830903225216002E-002, &
          2.3454574424525950E-003, 9.0830235263456302E-006, &
          1.4223646413363950E-007, 2.0830903246246202E-005, &
          4.2438745414205950E-005, 7.0830924624624602E-002, &
          5.4524652423463650E-006, 3.0152303238316002E-005, &
          2.1298526400363360E-003, 1.0046536363424200E-003)
   node_idx = 45601

   delta = [2.142431806556E-002, -1.194000280637E-001, -4.045444423650E-001]
   dist2 = sum(delta * delta)
   sig2 = 0.5

   force_law = 21
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 0.60973719248139013_kind_physics)
   call test(particle%results%af, -1.5452074863397835E-004_kind_physics)
   force_law = 22
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 1.5117595239300063_kind_physics)
   call test(particle%results%af, -4.4039727129395576E-004_kind_physics)
   force_law = 61
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 5.1779418086786198_kind_physics)
   call test(particle%results%af, -2.4269170113274392E-003_kind_physics)
   force_law = 62
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 8.7918855918903134_kind_physics)
   call test(particle%results%af, -4.0444647618862435E-003_kind_physics)
   force_law = 22
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 9.6881792618759057_kind_physics)
   call test(particle%results%af, -4.1731028164021314E-003_kind_physics)
   force_law = 62
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 15.634842375207171_kind_physics)
   call test(particle%results%af, -3.5482499572673954E-003_kind_physics)

   delta = [2.142431806556E-001, -1.194000280637E-002, -4.045444423650E-002]
   dist2 = sum(delta * delta)
   sig2 = 0.25

   force_law = 21
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 16.884920877826609_kind_physics)
   call test(particle%results%af, -3.9382934831662317E-003_kind_physics)
   force_law = 22
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 19.150838935027245_kind_physics)
   call test(particle%results%af, -4.4243081482363443E-003_kind_physics)
   force_law = 61
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 39.420917266392777_kind_physics)
   call test(particle%results%af, -1.1504300363356355E-003_kind_physics)
   force_law = 62
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 53.197615642744132_kind_physics)
   call test(particle%results%af, -1.3644748524379038E-003_kind_physics)
   force_law = 22
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 57.103468790835102_kind_physics)
   call test(particle%results%af, -3.6693810348134709E-003_kind_physics)
   force_law = 62
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 86.196062710981153_kind_physics)
   call test(particle%results%af, -4.3606871964010474E-002_kind_physics)

   delta = [2.142431806556E-001, -1.194000280637E-002, -4.045444423650E-002]
   dist2 = sum(delta * delta)
   sig2 = 0.0

   force_law = 21
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 101.83882023952870_kind_physics)
   !call test(particle%results%af, NaN)
   force_law = 22
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 117.48157776807625_kind_physics)
   !call test(particle%results%af, NaN)
   force_law = 61
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 133.12433529662380_kind_physics)
   !call test(particle%results%af, NaN)
   force_law = 62
   call calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 148.76709282517135_kind_physics)
   !call test(particle%results%af, NaN)
   force_law = 22
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 168.10032148599817_kind_physics)
   !call test(particle%results%af, NaN)
   force_law = 62
   call calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
   write (*, *) sum(particle%results%u), sum(particle%results%af)
   call test(particle%results%u, 187.43355014682496_kind_physics)
   !call test(particle%results%af, NaN)

contains
   subroutine test(vec, val)
      real(kind_physics), intent(in) :: vec(3)
      real(kind_physics), intent(in) :: val

      if (abs((sum(vec) - val) / val) .gt. 0.01) then
         write (*, *) 'Results outside 1% margin!'
         write (*, *) sum(vec), val
         stop 666
      end if
   end subroutine

end program test_interaction
#endif
