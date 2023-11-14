! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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
!> Encapsulates the low-level kernels for Coulomb- and similar interactions
!>
module module_coulomb_kernels
   use module_pepc_kinds
   use module_pepc_types
   use module_shortcut
   implicit none
   save
   private

   public calc_force_darwin_2D3V
   public calc_force_darwin_2D3V_direct
   public calc_force_darwin_3D
   public calc_force_darwin_3D_direct

contains

   !>
   !> Calculates 3D Coulomb interaction of particle p with tree node inode
   !> that is shifted by the lattice vector vbox
   !> results are returned in exyz, phi
   !>
   subroutine calc_force_coulomb_3D(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi

      real(kind_physics) :: rd, dx, dy, dz, r, dx2, dy2, dz2, dx3, dy3, dz3, rd2, rd3, rd5, rd7, fd1, fd2, fd3, fd4, fd5, fd6, m2rd5, m5rd7, &
                            pre1, pre2x, pre2y, pre2z, preQ1, preQ2, preQ3

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd = one / r
      rd2 = rd * rd
      rd3 = rd * rd2
      rd5 = rd3 * rd2
      rd7 = rd5 * rd2

      dx2 = dx * dx
      dy2 = dy * dy
      dz2 = dz * dz
      dx3 = dx * dx2
      dy3 = dy * dy2
      dz3 = dz * dz2

      fd1 = three * dx2 * rd5 - rd3
      fd2 = three * dy2 * rd5 - rd3
      fd3 = three * dz2 * rd5 - rd3
      fd4 = three * dx * dy * rd5
      fd5 = three * dy * dz * rd5
      fd6 = three * dx * dz * rd5

      m2rd5 = two * rd5
      m5rd7 = five * rd7
      pre1 = m5rd7 * dx * dy * dz
      pre2x = m5rd7 * dx2 - rd5
      pre2y = m5rd7 * dy2 - rd5
      pre2z = m5rd7 * dz2 - rd5
      preQ1 = pre2x * t%quad(1)
      preQ2 = pre2y * t%quad(2)
      preQ3 = pre2z * t%quad(3)

      phi = t%charge * rd &                                                     !&  monopole term
            - (dx * t%dip(1) + dy * t%dip(2) + dz * t%dip(3)) * rd3 &           !&  dipole
            + fd1 * t%quad(1) + fd2 * t%quad(2) + fd3 * t%quad(3) &             !&  quadrupole
            + fd4 * t%xyquad + fd5 * t%yzquad + fd6 * t%zxquad

      exyz(1) = t%charge * dx * rd3 &                                           !& monopole term
                - (fd1 * t%dip(1) + fd4 * t%dip(2) + fd6 * t%dip(3)) &          !& dipole term
                + three * ( &                                                   !& quadrupole term
                          dx * ( &                                              !&
                               (pre2x - m2rd5) * t%quad(1) &                    !&
                               + preQ2 &                                        !&
                               + preQ3 &                                        !&
                               ) &                                              !&
                          + dy * pre2x * t%xyquad &                             !&
                          + dz * pre2x * t%zxquad &                             !&
                          + pre1 * t%yzquad &                                   !&
                          )                                                     !&

      exyz(2) = t%charge * dy * rd3 &                                           !&
                - (fd2 * t%dip(2) + fd4 * t%dip(1) + fd5 * t%dip(3)) &          !&
                + three * ( &                                                   !&
                          dy * ( &                                              !&
                               (pre2y - m2rd5) * t%quad(2) &                    !&
                               + preQ1 &                                        !&
                               + preQ3 &                                        !&
                               ) &                                              !&
                          + dx * pre2y * t%xyquad &                             !&
                          + dz * pre2y * t%yzquad &                             !&
                          + pre1 * t%zxquad &                                   !&
                          )                                                     !&

      exyz(3) = t%charge * dz * rd3 &                                           !&
                - (fd3 * t%dip(3) + fd5 * t%dip(2) + fd6 * t%dip(1)) &          !&
                + three * ( &                                                   !&
                          dz * ( &                                              !&
                               +(pre2z - m2rd5) * t%quad(3) &                   !&
                               + preQ2 &                                        !&
                               + preQ1 &                                        !&
                               ) &                                              !&
                          + dx * pre2z * t%zxquad &                             !&
                          + dy * pre2z * t%yzquad &                             !&
                          + pre1 * t%xyquad &                                   !&
                          )                                                     !&
   end subroutine calc_force_coulomb_3D

   !>
   !> Calculates 2D Coulomb interaction of particle p with tree node inode
   !> that is shifted by the lattice vector vbox
   !> results are returned in exy, phi
   !> Unregularized force law is:
   !>   Phi = -2q log R
   !>   Ex = -dPhi/dx = 2 q x/R^2 etc
   !>
   subroutine calc_force_coulomb_2D(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(1:2), phi

      real(kind_physics) :: dx, dy, rd2, rd4, rd6, dx2, dy2, dx3, dy3

      dx = d(1)
      dy = d(2)

      rd2 = one / d2
      rd4 = rd2 * rd2
      rd6 = rd4 * rd2

      dx2 = dx * dx
      dy2 = dy * dy
      dx3 = dx2 * dx
      dy3 = dy2 * dy

      phi = -half * t%charge * log(d2) &                                     !  monopole term
            - (dx * t%dip(1) + dy * t%dip(2)) * rd2 &                        !  dipole
            + t%quad(1) * (two * dx2 * rd4 - rd2) &                          ! quadrupole
            + t%quad(2) * (two * dy2 * rd4 - rd2) &
            + two * t%xyquad * dx * dy * rd4

      exy(1) = t%charge * dx * rd2 &                                         ! monopole
               - t%dip(1) * (two * dx2 * rd4 - rd2) &                        ! dipole
               - t%dip(2) * two * dx * dy * rd4 &
               + t%quad(1) * (eight * dx3 * rd6 - six * dx * rd4) &          ! quadrupole
               + t%quad(2) * (eight * dx * dy2 * rd6 - two * dx * rd4) &
               + t%xyquad * (eight * dx2 * dy * rd6 - two * dy * rd4)

      exy(2) = t%charge * dy * rd2 &                                         ! monopole
               - t%dip(2) * (two * dy2 * rd4 - rd2) &                        ! dipole
               - t%dip(1) * two * dx * dy * rd4 &
               + t%quad(2) * (eight * dy3 * rd6 - six * dy * rd4) &          ! quadrupole
               + t%quad(1) * (eight * dy * dx2 * rd6 - two * dy * rd4) &
               + t%xyquad * (eight * dy2 * dx * rd6 - two * dx * rd4)

   end subroutine calc_force_coulomb_2D

   subroutine calc_force_darwin_2D(t, d, d2, eps2, phi, E, A, J, Jirr, B, Ex, Ey, Ax, Ay, Axx, Axy, Ayy)
      use module_tool, only: cross_product, double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  E(1:2), A(1:3), J(1:3), Jirr(1:3), phi, Ax(1:3), Ay(1:3), B, Axx(1:3), Axy(1:3), Ayy(1:3), Ex(1:2), Ey(1:2)

      real(kind_physics) :: dx, dy, rd2, rd4, rd6, rd8, dx2, dy2, dx3, dy3, r2, logR2e, over2, over4, over6, over8, &
                            phi_mono, phi_dip(1:2), phi_quad(1:3), rho_mono, rho_dip(1:2), rho_quad(1:3), &
                            E_mono(1:3), Ex_dip(1:3), Ey_dip(1:3), Exx_quad(1:3), Eyy_quad(1:3), logTmp, &
                            A1_mono, A1_dip(1:2), A1_quad(1:3), A2_mono, A2_dip(1:2), A2_quad(1:3), &
                            dist(1:3), Jir_mono, Jir_dip(1:2), Jir_quad(1:3), over, Exy_quad(1:3), &
                            A3_mono(1:3), A3x_dip(1:3), A3y_dip(1:3), A3xx_quad(1:3), A3yy_quad(1:3), &
                            A3xy_quad(1:3), A1_oct(1:4), A3xxx_oct(1:3), A3xxy_oct(1:3), &
                            A3xyy_oct(1:3), A3yyy_oct(1:3), Btmp(1:3), Exxx_oct(1:3), Exxy_oct(1:3), &
                            Exyy_oct(1:3), Eyyy_oct(1:3), A4_mono(1:3), A4x_dip(1:3), A4y_dip(1:3), &
                            A1_sed(1:5), A3xxxx_sed(1:3), A3xxxy_sed(1:3), A3xxyy_sed(1:3), &
                            A3xyyy_sed(1:3), A3yyyy_sed(1:3), over5, over7, over9, dx4, dy4, dx5, dy5, over10, &
                            A2_oct(1:4), A2_sed(1:5)

      dx = d(1)
      dy = d(2)
      dist(1:2) = d(1:2)
      dist(3) = zero
      r2 = dx**2 + dy**2
      over = sqrt(one / r2)
      over2 = one / r2
      over4 = over2 * over2
      over5 = over**5
      over6 = over4 * over2
      over7 = over**7
      over8 = over6 * over2
      over9 = over**9
      over10 = over8 * over2
      rd2 = one / d2
      rd4 = rd2 * rd2
      rd6 = rd4 * rd2
      rd8 = rd6 * rd2

      dx2 = dx * dx
      dy2 = dy * dy
      dx3 = dx2 * dx
      dy3 = dy2 * dy
      dx4 = dx3 * dx
      dy4 = dy3 * dy
      dx5 = dx4 * dx
      dy5 = dy4 * dy

      logTmp = log(r2 / eps2 + one)
      logR2e = eps2 * over2 * logTmp

      !!! Rho/J must be multiply by eps2/pi
      rho_mono = rd4
      rho_dip(1:2) = -four * d(1:2) * rd6
      rho_quad(1) = twentyfour * dx2 * rd8 - four * rd6
      rho_quad(2) = twentyfour * dy2 * rd8 - four * rd6
      rho_quad(3) = twentyfour * dy * dx * rd8

      !!! Jirr must be multiply by 2
      Jir_mono = rd2
      Jir_dip(1:2) = two * d(1:2) * rd4
      Jir_quad(1) = two * rd4 * (four * dx2 * rd2 - one)
      Jir_quad(2) = two * rd4 * (four * dy2 * rd2 - one)
      Jir_quad(3) = eight * dx * dy * rd6

      phi_mono = -log(d2)
      phi_dip = -two * d(1:2) * rd2
      phi_quad(1) = four * dx2 * rd4 - two * rd2
      phi_quad(2) = four * dy2 * rd4 - two * rd2
      phi_quad(3) = four * dy * dx * rd4

      E_mono(1:2) = -phi_dip
      E_mono(3) = zero
      Ex_dip(1:3) = -(/phi_quad(1), phi_quad(3), zero/)
      Ey_dip(1:3) = -(/phi_quad(3), phi_quad(2), zero/)

      Exx_quad(1) = sixteen * dx3 * rd6 - twelve * dx * rd4
      Exx_quad(2) = sixteen * dx2 * dy * rd6 - four * dy * rd4  !! Tshis is E_xy(1)
      Exx_quad(3) = zero
      Eyy_quad(1) = sixteen * dx * dy2 * rd6 - four * dx * rd4  !! This is E_xy(2)
      Eyy_quad(2) = sixteen * dy3 * rd6 - twelve * dy * rd4
      Eyy_quad(3) = zero

      Exy_quad(1) = Exx_quad(2)
      Exy_quad(2) = Eyy_quad(1)
      Exy_quad(3) = zero

      Exxx_oct(1) = ninetysix * dx2 * rd6 - twelve * rd4 - ninetysix * dx2**2 * rd8
      Exxx_oct(2) = fortyeight * dx * dy * rd6 - ninetysix * dx3 * dy * rd8
      Exxx_oct(3) = zero

      Exxy_oct(1) = fortyeight * dx * dy * rd6 - ninetysix * dx3 * dy * rd8
      Exxy_oct(2) = sixteen * dx2 * rd6 - four * rd4 + sixteen * dy2 * rd6 - ninety * dx2 * dy2 * rd8
      Exxy_oct(3) = zero

      Exyy_oct(1) = sixteen * dx2 * rd6 - four * rd4 + sixteen * dy2 * rd6 - ninety * dx2 * dy2 * rd8
      Exyy_oct(2) = fortyeight * dx * dy * rd6 - ninetysix * dx * dy3 * rd8
      Exyy_oct(3) = zero

      Eyyy_oct(1) = fortyeight * dx * dy * rd6 - ninetysix * dx * dy3 * rd8
      Eyyy_oct(2) = ninetysix * dy2 * rd6 - twelve * rd4 - ninetysix * dy2**2 * rd8
      Eyyy_oct(3) = zero

      A1_mono = logR2e - one
      A1_dip(1:2) = two * eps2 * d(1:2) * rd2 * over2 - two * over2 * logR2e * d(1:2)
      A1_quad(1) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dx2 * rd2 * over2 &
                                         - two * dx2 * rd4 + four * logTmp * dx2 * over4)
      A1_quad(2) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dy2 * rd2 * over2 &
                                         - two * dy2 * rd4 + four * logTmp * dy2 * over4)
      A1_quad(3) = four * over2 * dx * dy * (two * logR2e * over2 - eps2 * rd4 - two * eps2 * rd2 * over2)

      A1_oct(1) = twelve * eps2 * dx * over2 * (four * dx2 * rd2 * over4 - two * rd2 * over2 - rd4 + &
                                                two * logTmp * over4 + two * dx2 * rd4 * over2 - four * dx2 * logTmp * over6) + sixteen * eps2 * dx3 * rd6 * over2 ! xxx

      A1_oct(2) = four * eps2 * dy * over2 * (two * logTmp * over4 - rd4 - two * rd2 * over2 + twelve * dx2 * rd2 * over4 + &
                                              six * dx2 * rd4 * over2 + four * dx2 * rd6 - twelve * dx2 * logTmp * over6)                                ! xxy

      A1_oct(3) = four * eps2 * dx * over2 * (two * logTmp * over4 - rd4 - two * rd2 * over2 + twelve * dy2 * rd2 * over4 + &
                                              six * dy2 * rd4 * over2 + four * dy2 * rd6 - twelve * dy2 * logTmp * over6)                                ! xxy

      A1_oct(4) = twelve * eps2 * dy * over2 * (four * dy2 * rd2 * over4 - two * rd2 * over2 - rd4 + &
                                                two * logTmp * over4 + two * dy2 * rd4 * over2 - four * dy2 * logTmp * over6) + sixteen * eps2 * dy3 * rd6 * over2 ! yyy

      A1_sed(1) = twentyfour * twelve * dx2 * over2 * rd2 * eps2 - twentyfour * eps2 * over4 * rd2 - twentyfour * sixteen * eps2 * dx4 * over8 * rd2 & !xxxx
                  - twelve * eps2 * over2 * rd4 + twentyfour * eps2 * over6 * logTmp - twentyfour * twelve * dx2 * eps2 * over8 * logTmp &
                  + twentyfour * sixteen * eps2 * dx4 * over5 * logTmp + twelve**2 * dx2 * over4 * eps2 * rd4 + ninetysix * dx2 * over2 * rd6 &
                  - twentyfour * eight * dx4 * eps2 * rd4 * over6 - two * eight**2 * dx4 * eps2 * rd6 * over4 - ninetysix * dx4 * eps2 * rd8 * over2

      A1_sed(2) = twelve**2 * dx * dy * eps2 * rd2 * over6 - twentyfour * sixteen * eps2 * dx3 * dy * rd2 * over8 - twelve**2 * dx * dy * eps2 * over8 * logTmp & !xxxy
                  - twentyfour * eight * dx3 * dy * eps2 * rd6 * over4 - ninetysix * dx3 * dy * eps2 * rd8 * over2 + twentyfour * sixteen * eps2 * dx3 * dy * logTmp * over10 &
                  + nine * eight * dx * dy * eps2 * rd4 * over4 + fortyeight * dx * dy * eps2 * rd6 * over2

      A1_sed(3) = fortyeight * eps2 * rd2 * over4 - eight * eps2 * rd2 * over4 - four * eps2 * rd4 * over2 + eight * eps2 * logTmp * over6 & !xxyy
                  - fortyeight * eps2 * logTmp * over6 - twentyfour * sixteen * dx2 * dy2 * eps2 * rd2 * over8 + twentyfour * eps2 * over2 * rd4 &
                  + sixteen * eps2 * over2 * rd6 - twentyfour * eight * dx2 * dy2 * eps2 * rd4 * over6 - two * eight**2 * dx2 * dy2 * eps2 * rd6 * over4 &
                  - ninetysix * dx2 * dy2 * eps2 * rd8 * over2 + twentyfour * sixteen * dx2 * dy2 * eps2 * logTmp * over10

      A1_sed(4) = twelve**2 * dx * dy * eps2 * rd2 * over6 - twentyfour * sixteen * eps2 * dx * dy3 * rd2 * over8 - twelve**2 * dx * dy * eps2 * over8 * logTmp & !xyyy
                  - twentyfour * eight * dx * dy3 * eps2 * rd6 * over4 - ninetysix * dx * dy3 * eps2 * rd8 * over2 + twentyfour * sixteen * eps2 * dx * dy3 * logTmp * over10 &
                  + nine * eight * dx * dy * eps2 * rd4 * over4 + fortyeight * dx * dy * eps2 * rd6 * over2

      A1_sed(5) = twentyfour * twelve * dy2 * over2 * rd2 * eps2 - twentyfour * eps2 * over4 * rd2 - twentyfour * sixteen * eps2 * dy4 * over8 * rd2 & !yyyy
                  - twelve * eps2 * over2 * rd4 + twentyfour * eps2 * over6 * logTmp - twentyfour * twelve * dy2 * eps2 * over8 * logTmp &
                  + twentyfour * sixteen * eps2 * dy4 * over5 * logTmp + twelve**2 * dy2 * over4 * eps2 * rd4 + ninetysix * dy2 * over2 * rd6 &
                  - twentyfour * eight * dy4 * eps2 * rd4 * over6 - two * eight**2 * dy4 * eps2 * rd6 * over4 - ninetysix * dy4 * eps2 * rd8 * over2

      A2_mono = -A1_mono + one + phi_mono + log(eps2)
      A2_dip(1:2) = -A1_dip(1:2) + phi_dip(1:2)
      A2_quad(1:3) = -A1_quad(1:3) + phi_quad(1:3)
      A2_oct(1) = -A1_oct(1) - Exx_quad(1)
      A2_oct(2) = -A1_oct(2) - Exy_quad(1)
      A2_oct(3) = -A1_oct(2) - Eyy_quad(1)
      A2_oct(4) = -A1_oct(3) - Eyy_quad(2)
      A2_sed(1) = -A1_sed(1) - Exxx_oct(1)
      A2_sed(2) = -A1_sed(2) - Exxy_oct(1)
      A2_sed(3) = -A1_sed(3) - Exyy_oct(1)
      A2_sed(4) = -A1_sed(4) - Eyyy_oct(1)
      A2_sed(5) = -A1_sed(5) - Eyyy_oct(2)

      A3_mono(1:3) = over * dist(1:3)
      A3x_dip(1) = over - dx2 * over * over2
      A3x_dip(2) = -dx * dy * over * over2
      A3x_dip(3) = zero
      A3y_dip(1) = -dx * dy * over * over2
      A3y_dip(2) = over - dy2 * over * over2
      A3y_dip(3) = zero
      A3xx_quad(1) = -three * over * over2 * dx * (one - dx2 * over2)
      A3xx_quad(2) = -over * over2 * dy * (one - three * dx2 * over2)   ! This is equal to A3_xy(1)
      A3xx_quad(3) = zero
      A3yy_quad(1) = -over * over2 * dx * (one - three * dy2 * over2)     ! This is equal to A3_xy(2)
      A3yy_quad(2) = -three * over * over2 * dy * (one - dy2 * over2)
      A3yy_quad(3) = zero
      A3xy_quad(1) = A3xx_quad(2)
      A3xy_quad(2) = A3yy_quad(1)
      A3xy_quad(3) = zero
      A3xxx_oct(1) = three * over**3 * (six * dx2 * over2 - one - five * dx2**2 * over4)
      A3xxx_oct(2) = three * over**5 * dx * dy * (three - five * dx2 * over2)
      A3xxx_oct(3) = zero

      A3xxy_oct(1) = three * over**5 * dx * dy * (three - five * dx2 * over2)
      A3xxy_oct(2) = over**3 * (two - fifteen * dx2 * dy2 * over4)
      A3xxy_oct(3) = zero

      A3xyy_oct(1) = over**3 * (two - fifteen * dx2 * dy2 * over4)
      A3xyy_oct(2) = three * over**5 * dx * dy * (three - five * dy2 * over2)
      A3xyy_oct(3) = zero

      A3yyy_oct(1) = three * over**5 * dx * dy * (three - five * dy2 * over2)
      A3yyy_oct(2) = three * over**3 * (six * dy2 * over2 - one - five * dy2**2 * over4)
      A3yyy_oct(3) = zero

      A3xxxx_sed(1) = fortyfive * dx * over5 - hundredfifty * dx3 * over7 + hundredfive * dx5 * over9
      A3xxxx_sed(2) = nine * dy * over5 - ninety * dx2 * dy * over7 + hundredfive * dx4 * dy * over9
      A3xxxx_sed(3) = zero

      A3xxxy_sed(1) = nine * dy * over5 - ninety * dx2 * dy * over7 + hundredfive * dx4 * dy * over9
      A3xxxy_sed(2) = nine * dx * over5 - fifteen * dx3 * over7 + hundredfive * dx3 * dy2 * over9 - fortyfive * dx * dy2 * over7
      A3xxxy_sed(3) = zero

      A3xxyy_sed(1) = nine * dx * over5 - fifteen * dx3 * over7 + hundredfive * dx3 * dy2 * over9 - fortyfive * dx * dy2 * over7
      A3xxyy_sed(2) = nine * dy * over5 - fifteen * dy3 * over7 + hundredfive * dx2 * dy3 * over9 - fortyfive * dx2 * dy * over7
      A3xxyy_sed(3) = zero

      A3xyyy_sed(1) = nine * dy * over5 - fifteen * dy3 * over7 + hundredfive * dx2 * dy3 * over9 - fortyfive * dx2 * dy * over7
      A3xyyy_sed(2) = nine * dx * over5 - ninety * dx * dy2 * over7 + hundredfive * dx * dy4 * over9
      A3xyyy_sed(3) = zero

      A3yyyy_sed(1) = nine * dx * over5 - ninety * dx * dy2 * over7 + hundredfive * dx * dy4 * over9
      A3yyyy_sed(2) = fortyfive * dy * over5 - hundredfifty * dy3 * over7 + hundredfive * dy5 * over9
      A3yyyy_sed(3) = zero

      A4_mono(1:3) = dist(1:3)
      A4x_dip(1) = one
      A4x_dip(2) = zero
      A4x_dip(3) = zero
      A4y_dip(1) = zero
      A4y_dip(2) = one
      A4y_dip(3) = zero

      phi = (t%charge * phi_mono) &! Monopole
            + (dot_product(t%dip(1:2), phi_dip(1:2))) &! Dipole
            + (dot_product(t%quad(1:2), phi_quad(1:2))) &! Quadrupole
            + (phi_quad(3) * t%xyquad)                                                                 ! Quadrupole
      phi = half * phi

      E(1:2) = (t%charge * E_mono(1:2)) &! Monopole
               + (t%dip(1) * Ex_dip(1:2) + t%dip(2) * Ey_dip(1:2)) &! Dipole
               + (t%quad(1) * Exx_quad(1:2) + t%quad(2) * Eyy_quad(1:2)) &! Quadrupole
               + (t%xyquad * Exy_quad(1:2))                                                               ! Quadrupole

      E(1:2) = half * E(1:2)

      Ex(1:2) = (t%charge * Ex_dip(1:2)) &! Monopole
                + (t%dip(1) * Exx_quad(1:2) + t%dip(2) * Exy_quad(1:2)) &! Dipole
                + (t%quad(1) * Exxx_oct(1:2) + t%quad(2) * Exyy_oct(1:2)) &! Quadrupole
                + (t%xyquad * Exxy_oct(1:2))                                                               ! Quadrupole

      Ex(1:2) = half * Ex(1:2)

      Ey(1:2) = (t%charge * Ey_dip(1:2)) &! Monopole
                + (t%dip(1) * Exy_quad(1:2) + t%dip(2) * Eyy_quad(1:2)) &! Dipole
                + (t%quad(1) * Exxy_oct(1:2) + t%quad(2) * Eyyy_oct(1:2)) &! Quadrupole
                + (t%xyquad * Exyy_oct(1:2))                                                               ! Quadrupole

      Ey(1:2) = half * Ey(1:2)

      J(1:3) = (t%monoj(1:3) * rho_mono) &! Monopole
               + (rho_dip(1) * t%dipjx(1:3) + rho_dip(2) * t%dipjy(1:3)) &! Dipole
               + (t%quadjx(1:3) * rho_quad(1) + t%quadjy(1:3) * rho_quad(2)) &! Quadrupole
               + (t%quadjxy(1:3) * rho_quad(3))                                                          ! Quadrupole

      J(1:3) = half * eps2 / pi * J(1:3)

      Jirr(1:3) = (t%monoj(1:3) * (three * half * Jir_mono - eps2 * rho_mono) &
                   + double_cross_product_left(E_mono(1:3), E_mono(1:3), t%monoj(1:3))) &! Monopole
                  + (t%dipjx(1:3) * (three * half * Jir_dip(1) - eps2 * rho_dip(1))) &! Dipole - x
                  + (double_cross_product_left(Ex_dip(1:3), E_mono(1:3), t%dipjx(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Ex_dip(1:3), t%dipjx(1:3))) &
                  + (t%dipjy(1:3) * (three * half * Jir_dip(2) - eps2 * rho_dip(2))) &! Dipole - y
                  + (double_cross_product_left(Ey_dip(1:3), E_mono(1:3), t%dipjy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Ey_dip(1:3), t%dipjy(1:3))) &
                  + (t%quadjx(1:3) * (three * half * Jir_quad(1) - eps2 * rho_quad(1))) &! Quadrupole - xx
                  + (double_cross_product_left(Exx_quad(1:3), E_mono(1:3), t%quadjx(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Exx_quad(1:3), t%quadjx(1:3))) &
                  + (two * double_cross_product_left(Ex_dip(1:3), Ex_dip(1:3), t%quadjx(1:3))) &
                  + (t%quadjy(1:3) * (three * half * Jir_quad(2) - eps2 * rho_quad(2))) &! Quadrupole - yy
                  + (double_cross_product_left(Eyy_quad(1:3), E_mono(1:3), t%quadjy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Eyy_quad(1:3), t%quadjy(1:3))) &
                  + (two * double_cross_product_left(Ey_dip(1:3), Ey_dip(1:3), t%quadjy(1:3))) &
                  + (t%quadjxy(1:3) * (three * half * Jir_quad(3) - eps2 * rho_quad(3))) &! Quadrupole - xy
                  + (double_cross_product_left(Exy_quad(1:3), E_mono(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Exy_quad(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(Ex_dip(1:3), Ey_dip(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(Ey_dip(1:3), Ex_dip(1:3), t%quadjxy(1:3)))

      Jirr(1:3) = half * oneoverpi * Jirr(1:3)

      Btmp(1:3) = cross_product(t%monoj(1:3), E_mono(1:3)) &! Monopole
                  + (cross_product(t%dipjx(1:3), Ex_dip(1:3))) &! Dipole x
                  + (cross_product(t%dipjy(1:3), Ey_dip(1:3))) &! Dipole y
                  + (cross_product(t%quadjx(1:3), Exx_quad(1:3))) &! Quadrupole xx
                  + (cross_product(t%quadjy(1:3), Eyy_quad(1:3))) &! Quadrupole yy
                  + (cross_product(t%quadjxy(1:3), Exy_quad(1:3)))                                          ! Quadrupole xy

      Btmp(1:3) = half * Btmp(1:3)

      A(1:3) = half * t%monoj(1:3) * (one - A1_mono - log(d2) + log(eps2)) &
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3)) &! Monopole
               + half * (t%dipjx(1:3) * (-A1_dip(1) + phi_dip(1))) &! Dipole - x
               + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
               + A1_mono * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
               + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3))) &
               + half * (t%dipjy(1:3) * (-A1_dip(2) + phi_dip(2))) &! Dipole - y
               + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
               + A1_mono * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
               + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3))) &
               + half * (t%quadjx(1:3) * (-A1_quad(1) + phi_quad(1))) &! Quadrupole - xx
               + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
               + A1_mono * (two * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
               !
               + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
               + two * A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
               + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3)) &
               + half * (t%quadjy(1:3) * (-A1_quad(2) + phi_quad(2))) &! Quadrupole - yy
               + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
               + A1_mono * (two * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
               + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
               + two * A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
               + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3)) &
               + half * (t%quadjxy(1:3) * (-A1_quad(3) + phi_quad(3))) &! Quadrupole - xy
               + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
               + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
               + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
               !
               + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
               + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
               + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
               + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
               + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))
!
      A(3) = t%monoj(3) * phi_mono &! Monopole
             + t%dipjx(3) * phi_dip(1) + t%dipjy(3) * phi_dip(2) &! Dipole
             + t%quadjx(3) * phi_quad(1) + t%quadjy(3) * phi_quad(2) &! Quadrupole
             + t%quadjxy(3) * phi_quad(3)                                                               ! Quadrupole

      A(1:3) = half * A(1:3)

      Ax(1:3) = half * (t%monoj(1:3) * (-A1_dip(1) + phi_dip(1))) &! Monopole
                + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3))) &
                + half * (t%dipjx(1:3) * (-A1_quad(1) + phi_quad(1))) &! Dipole - x
                + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                + two * A1_mono * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                + half * (t%dipjy(1:3) * (-A1_quad(3) + phi_quad(3))) &! Dipole - y
                + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                + half * (t%quadjx(1:3) * (-A1_oct(1) - Exx_quad(1))) &! Quadrupole - xx
                + A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + A1_oct(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + three * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + three * A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + three * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + three * A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + six * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + half * (t%quadjy(1:3) * (-A1_oct(3) - Exy_quad(2))) &! Quadrupole - yy
                + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + two * A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + half * (t%quadjxy(1:3) * (-A1_oct(2) - Exy_quad(1))) &! Quadrupole - xy
                + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + two * A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3)))

      Ax(3) = t%monoj(3) * phi_dip(1) &! Monopole
              + t%dipjx(3) * phi_quad(1) + t%dipjy(3) * phi_quad(3) &! Dipole
              - t%quadjx(3) * Exx_quad(1) - t%quadjy(3) * Eyy_quad(1) &! Quadrupole
              - t%quadjxy(3) * Exx_quad(2)                                                               ! Quadrupole

      Ay(1:3) = half * (t%monoj(1:3) * (-A1_dip(2) + phi_dip(2))) &! Monopole
                + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3))) &
                + half * (t%dipjx(1:3) * (-A1_quad(3) + phi_quad(3))) &! Dipole - x
                + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                + half * (t%dipjy(1:3) * (-A1_quad(2) + phi_quad(2))) &! Dipole - y
                + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                + two * A1_mono * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                + half * (t%quadjx(1:3) * (-A1_oct(2) - Exy_quad(1))) &! Quadrupole - xx
                + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + two * A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                + half * (t%quadjy(1:3) * (-A1_oct(4) - Eyy_quad(2))) &! Quadrupole - yy
                + A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + A1_oct(4) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + three * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + three * A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + three * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                + three * A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + six * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                + half * (t%quadjxy(1:3) * (-A1_oct(3) - Exy_quad(2))) &! Quadrupole - xy
                + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                + two * A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + two * A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + two * A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                + two * A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3)))

      Ay(3) = t%monoj(3) * phi_dip(2) &! Monopole
              + t%dipjx(3) * phi_quad(3) + t%dipjy(3) * phi_quad(2) &! Dipole
              - t%quadjx(3) * Exx_quad(2) - t%quadjy(3) * Eyy_quad(2) &! Quadrupole
              - t%quadjxy(3) * Exy_quad(2)                                                                 ! Quadrupole

      Ax(1:3) = half * Ax(1:3)
      Ay(1:3) = half * Ay(1:3)

!
!
      Axx(1:3) = half * (t%monoj(1:3) * A2_quad(1)) &! Monopole
                 + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + half * (t%dipjx(1:3) * A2_oct(1)) &! Dipole - x
                 + A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3)) &
                 + A1_oct(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + two * A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + half * (t%dipjy(1:3) * A2_oct(2)) &! Dipole - y
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_quad(3) * double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4x_dip(1:3)) &
                 + A1_oct(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + half * (t%quadjx(1:3) * A2_sed(1)) &! Quadrupole - xx
                 + A1_mono * (double_cross_product_left(A3xxxx_sed(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_sed(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + three * A1_quad(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + three * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + three * A1_oct(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + three * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + three * A1_oct(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + three * A1_dip(1) * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + three * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + three * A1_quad(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + three * A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + three * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + six * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + six * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + half * (t%quadjy(1:3) * A2_sed(3)) &! Quadrupole - yy
                 + A1_mono * (double_cross_product_left(A3xxyy_sed(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_sed(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + half * (t%quadjxy(1:3) * A2_sed(2)) &! Quadrupole - xy
                 + A1_mono * (double_cross_product_left(A3xxxy_sed(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_sed(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xxx_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3)))

!

      Axx(3) = -t%monoj(3) * Ex_dip(1) &! Monopole
               - t%dipjx(3) * Exx_quad(1) - t%dipjy(3) * Exx_quad(2) &! Dipole
               - t%quadjx(3) * Exxx_oct(1) - t%quadjy(3) * Exyy_oct(1) &! Quadrupole
               - t%quadjxy(3) * Exxx_oct(2)                                                                ! Quadrupole

      Ayy(1:3) = half * (t%monoj(1:3) * A2_quad(2)) &! Monopole
                 + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + half * (t%dipjx(1:3) * A2_oct(3)) &! Dipole - x
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_quad(3) * double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4y_dip(1:3)) &
                 + A1_oct(3) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%dipjx(1:3), A3x_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + half * (t%dipjy(1:3) * A2_oct(4)) &! Dipole - y
                 + A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A3_mono(1:3))) &
                 + A1_quad(2) * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3)) &
                 + A1_oct(3) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + two * A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + half * (t%quadjx(1:3) * A2_sed(3)) &! Quadrupole - xx
                 + A1_mono * (double_cross_product_left(A3xxyy_sed(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_sed(3) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + half * (t%quadjy(1:3) * A2_sed(5)) &! Quadrupole - yy
                 + A1_mono * (double_cross_product_left(A3yyyy_sed(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(4) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(4) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + A1_sed(5) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_oct(4) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_oct(4) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_dip(2) * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + six * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + six * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + half * (t%quadjxy(1:3) * A2_sed(4)) &! Quadrupole - xy
                 + A1_mono * (double_cross_product_left(A3xyyy_sed(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_sed(4) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_oct(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3)))

      Ayy(3) = -t%monoj(3) * Ey_dip(2) &! Monopole
               - t%dipjx(3) * Eyy_quad(1) - t%dipjy(3) * Eyy_quad(2) &! Dipole
               - t%quadjx(3) * Exxy_oct(2) - t%quadjy(3) * Eyyy_oct(2) &! Quadrupole
               - t%quadjxy(3) * Exyy_oct(2)

      Axy(1:3) = half * (t%monoj(1:3) * A2_quad(3)) &! Monopole
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A4_mono(1:3))) &
                 + half * (t%dipjx(1:3) * A2_oct(2)) &! Dipole - x
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_quad(3) * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3)) &
                 + A1_oct(2) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4_mono(1:3)) &
                 + A1_mono * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A4x_dip(1:3))) &
                 + half * (t%dipjy(1:3) * A2_oct(3)) &! Dipole - y
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4x_dip(1:3)) &
                 + A1_oct(3) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4_mono(1:3)) &
                 + two * A1_mono * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A4y_dip(1:3))) &
                 + half * (t%quadjx(1:3) * A2_sed(2)) &! Quadrupole - xx
                 + A1_mono * (double_cross_product_left(A3xxxy_sed(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_sed(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_quad(3) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + A1_oct(1) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4_mono(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A4x_dip(1:3))) &
                 + half * (t%quadjy(1:3) * A2_sed(4)) &! Quadrupole - yy
                 + A1_mono * (double_cross_product_left(A3xyyy_sed(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(4) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + A1_oct(4) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + A1_sed(4) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + three * A1_oct(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_dip(2) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_dip(2) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4x_dip(1:3))) &
                 + three * A1_quad(3) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4_mono(1:3))) &
                 + three * A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + three * A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + six * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + six * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A4y_dip(1:3))) &
                 + half * (t%quadjxy(1:3) * A2_sed(3)) &! Quadrupole - xy
                 + A1_mono * (double_cross_product_left(A3xxyy_sed(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_sed(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_quad(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_mono * (double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_dip(1) * (double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + A1_quad(2) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + A1_oct(3) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4_mono(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4x_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(1) * (double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_mono * (double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(1) * (double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_oct(2) * (double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_dip(2) * (double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A4y_dip(1:3))) &
                 + two * A1_quad(3) * (double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A4y_dip(1:3)))

      Axy(3) = -t%monoj(3) * Ex_dip(2) &! Monopole
               - t%dipjx(3) * Exx_quad(2) - t%dipjy(3) * Eyy_quad(1) &! Dipole
               - t%quadjx(3) * Exxx_oct(2) - t%quadjy(3) * Exyy_oct(2) &! Quadrupole
               - t%quadjxy(3) * Exxy_oct(2)                                                                 ! Quadrupole

      Axx(1:3) = half * Axx(1:3)
      Ayy(1:3) = half * Ayy(1:3)
      Axy(1:3) = half * Axy(1:3)
!
!

      J(3) = zero
      Jirr(3) = zero
      B = Btmp(3)
      A(3) = zero
      Ax(3) = zero
      Ay(3) = zero
      Axx(3) = zero
      Axy(3) = zero
      Ayy(3) = zero

   end subroutine calc_force_darwin_2D

!    subroutine calc_force_darwin_2D3V(t, d, d2,eps2, phi, E, A, J, Jirr, B, Ex, Ey, Ax, Ay, Axx, Axy , Ayy )
   subroutine calc_force_darwin_2D3V(t, d, d2, eps2, phi, E, A, J, Jirr, B, Ax, Ay)
      use module_tool, only: cross_product, double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
!      real(kind_physics), intent(out) ::  E(1:2),A(1:3),J(1:3),Jirr(1:3),B(1:3),phi,Ax(1:3),Ay(1:3),Axx(1:3),Ayy(1:3),Axy(1:3),&
!                                          Ex(1:2),Ey(1:2)
      real(kind_physics), intent(out) ::  E(1:2), A(1:3), J(1:3), Jirr(1:3), B(1:3), phi, Ax(1:3), Ay(1:3)

      real(kind_physics) :: dx, dy, rd2, rd4, rd6, rd8, dx2, dy2, dx3, dy3, r2, logR2e, over2, over4, over6, over8, &
                            phi_mono, phi_dip(1:2), phi_quad(1:3), rho_mono, rho_dip(1:2), rho_quad(1:3), &
                            E_mono(1:3), Ex_dip(1:3), Ey_dip(1:3), Exx_quad(1:3), Eyy_quad(1:3), logTmp, &
                            A1_mono, A1_dip(1:2), A1_quad(1:3), A2_mono, A2_dip(1:2), A2_quad(1:3), &
                            dist(1:3), Jir_mono, Jir_dip(1:2), Jir_quad(1:3), over, Exy_quad(1:3), &
                            A3_mono(1:3), A3x_dip(1:3), A3y_dip(1:3), A3xx_quad(1:3), A3yy_quad(1:3), &
                            A3xy_quad(1:3), A1_oct(1:4), A3xxx_oct(1:3), A3xxy_oct(1:3), &
                            A3xyy_oct(1:3), A3yyy_oct(1:3), Btmp(1:3), Exxx_oct(1:3), Exxy_oct(1:3), &
                            Exyy_oct(1:3), Eyyy_oct(1:3), &!,A4_mono(1:3),A4x_dip(1:3),A4y_dip(1:3),      &
                            A1_sed(1:5), A3xxxx_sed(1:3), A3xxxy_sed(1:3), A3xxyy_sed(1:3), &
                            A3xyyy_sed(1:3), A3yyyy_sed(1:3), over5, over7, over9, dx4, dy4, dx5, dy5, over10, &
                            A2_oct(1:4), A2_sed(1:5)

      dx = d(1)
      dy = d(2)
      dist(1:2) = d(1:2)
      dist(3) = zero
      r2 = dx**2 + dy**2
      over = sqrt(one / r2)
      over2 = one / r2
      over4 = over2 * over2
      over5 = over**5
      over6 = over4 * over2
      over7 = over**7
      over8 = over6 * over2
      over9 = over**9
      over10 = over8 * over2
      rd2 = one / d2
      rd4 = rd2 * rd2
      rd6 = rd4 * rd2
      rd8 = rd6 * rd2

      dx2 = dx * dx
      dy2 = dy * dy
      dx3 = dx2 * dx
      dy3 = dy2 * dy
      dx4 = dx3 * dx
      dy4 = dy3 * dy
      dx5 = dx4 * dx
      dy5 = dy4 * dy

      logTmp = log(r2 / eps2 + one)
      logR2e = eps2 * over2 * logTmp

      !!! Rho/J must be multiply by eps2/pi
      rho_mono = rd4
      rho_dip(1:2) = -four * d(1:2) * rd6
      rho_quad(1) = twentyfour * dx2 * rd8 - four * rd6
      rho_quad(2) = twentyfour * dy2 * rd8 - four * rd6
      rho_quad(3) = twentyfour * dy * dx * rd8

      !!! Jirr must be multiply by 2
      Jir_mono = rd2
      Jir_dip(1:2) = two * d(1:2) * rd4
      Jir_quad(1) = two * rd4 * (four * dx2 * rd2 - one)
      Jir_quad(2) = two * rd4 * (four * dy2 * rd2 - one)
      Jir_quad(3) = eight * dx * dy * rd6

      phi_mono = -log(d2)
      phi_dip = -two * d(1:2) * rd2
      phi_quad(1) = four * dx2 * rd4 - two * rd2
      phi_quad(2) = four * dy2 * rd4 - two * rd2
      phi_quad(3) = four * dy * dx * rd4

      E_mono(1:2) = -phi_dip
      E_mono(3) = zero
      Ex_dip(1:3) = -(/phi_quad(1), phi_quad(3), zero/)
      Ey_dip(1:3) = -(/phi_quad(3), phi_quad(2), zero/)

      Exx_quad(1) = sixteen * dx3 * rd6 - twelve * dx * rd4
      Exx_quad(2) = sixteen * dx2 * dy * rd6 - four * dy * rd4  !! This is E_xy(1)
      Exx_quad(3) = zero
      Eyy_quad(1) = sixteen * dx * dy2 * rd6 - four * dx * rd4  !! This is E_xy(2)
      Eyy_quad(2) = sixteen * dy3 * rd6 - twelve * dy * rd4
      Eyy_quad(3) = zero

      Exy_quad(1) = Exx_quad(2)
      Exy_quad(2) = Eyy_quad(1)
      Exy_quad(3) = zero

      Exxx_oct(1) = ninetysix * dx2 * rd6 - twelve * rd4 - ninetysix * dx4 * rd8
      Exxx_oct(2) = fortyeight * dx * dy * rd6 - ninetysix * dx3 * dy * rd8
      Exxx_oct(3) = zero

      Exxy_oct(1) = fortyeight * dx * dy * rd6 - ninetysix * dx3 * dy * rd8
      Exxy_oct(2) = sixteen * dx2 * rd6 - four * rd4 + sixteen * dy2 * rd6 - ninety * dx2 * dy2 * rd8
      Exxy_oct(3) = zero

      Exyy_oct(1) = sixteen * dx2 * rd6 - four * rd4 + sixteen * dy2 * rd6 - ninety * dx2 * dy2 * rd8
      Exyy_oct(2) = fortyeight * dx * dy * rd6 - ninetysix * dx * dy3 * rd8
      Exyy_oct(3) = zero

      Eyyy_oct(1) = fortyeight * dx * dy * rd6 - ninetysix * dx * dy3 * rd8
      Eyyy_oct(2) = ninetysix * dy2 * rd6 - twelve * rd4 - ninetysix * dy4 * rd8
      Eyyy_oct(3) = zero

      A1_mono = logR2e - one
      A1_dip(1:2) = two * eps2 * d(1:2) * rd2 * over2 - two * over2 * logR2e * d(1:2)
      A1_quad(1) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dx2 * rd2 * over2 &
                                         - two * dx2 * rd4 + four * logTmp * dx2 * over4)
      A1_quad(2) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dy2 * rd2 * over2 &
                                         - two * dy2 * rd4 + four * logTmp * dy2 * over4)
      A1_quad(3) = four * over2 * dx * dy * (two * logR2e * over2 - eps2 * rd4 - two * eps2 * rd2 * over2)

      A1_oct(1) = twelve * eps2 * dx * over2 * (four * dx2 * rd2 * over4 - two * rd2 * over2 - rd4 + &
                                                two * logTmp * over4 + two * dx2 * rd4 * over2 - four * dx2 * logTmp * over6) + sixteen * eps2 * dx3 * rd6 * over2 ! xxx

      A1_oct(2) = four * eps2 * dy * over2 * (two * logTmp * over4 - rd4 - two * rd2 * over2 + twelve * dx2 * rd2 * over4 + &
                                              six * dx2 * rd4 * over2 + four * dx2 * rd6 - twelve * dx2 * logTmp * over6)                                ! xxy

      A1_oct(3) = four * eps2 * dx * over2 * (two * logTmp * over4 - rd4 - two * rd2 * over2 + twelve * dy2 * rd2 * over4 + &
                                              six * dy2 * rd4 * over2 + four * dy2 * rd6 - twelve * dy2 * logTmp * over6)                                ! xxy

      A1_oct(4) = twelve * eps2 * dy * over2 * (four * dy2 * rd2 * over4 - two * rd2 * over2 - rd4 + &
                                                two * logTmp * over4 + two * dy2 * rd4 * over2 - four * dy2 * logTmp * over6) + sixteen * eps2 * dy3 * rd6 * over2 ! yyy

      A1_sed(1) = twentyfour * twelve * dx2 * over6 * rd2 * eps2 - twentyfour * eps2 * over4 * rd2 - twentyfour * sixteen * eps2 * dx4 * over8 * rd2 & !xxxx
                  - twelve * eps2 * over2 * rd4 + twentyfour * eps2 * over6 * logTmp - twentyfour * twelve * dx2 * eps2 * over8 * logTmp &
                  + twentyfour * sixteen * eps2 * dx4 * over10 * logTmp + twelve**2 * dx2 * over4 * eps2 * rd4 + ninetysix * dx2 * eps2 * over2 * rd6 &
                  - twentyfour * eight * dx4 * eps2 * rd4 * over6 - two * eight**2 * dx4 * eps2 * rd6 * over4 - ninetysix * dx4 * eps2 * rd8 * over2

      A1_sed(2) = twelve**2 * dx * dy * eps2 * rd2 * over6 - twentyfour * sixteen * eps2 * dx3 * dy * rd2 * over8 - twelve**2 * dx * dy * eps2 * over8 * logTmp & !xxxy
                  - twentyfour * eight * dx3 * dy * eps2 * rd4 * over6 - two * eight**2 * dx3 * dy * eps2 * rd6 * over4 - ninetysix * dx3 * dy * eps2 * rd8 * over2 &
                  + twentyfour * sixteen * eps2 * dx3 * dy * logTmp * over10 + nine * eight * dx * dy * eps2 * rd4 * over4 + fortyeight * dx * dy * eps2 * rd6 * over2

      A1_sed(3) = fortyeight * eps2 * rd2 * over4 - eight * eps2 * rd2 * over4 - four * eps2 * rd4 * over2 + eight * eps2 * logTmp * over6 & !xxyy
                  - fortyeight * eps2 * logTmp * over6 - twentyfour * sixteen * dx2 * dy2 * eps2 * rd2 * over8 + twentyfour * eps2 * over2 * rd4 &
                  + sixteen * eps2 * rd6 - twentyfour * eight * dx2 * dy2 * eps2 * rd4 * over6 - two * eight**2 * dx2 * dy2 * eps2 * rd6 * over4 &
                  - ninetysix * dx2 * dy2 * eps2 * rd8 * over2 + twentyfour * sixteen * dx2 * dy2 * eps2 * logTmp * over10

      A1_sed(4) = twelve**2 * dx * dy * eps2 * rd2 * over6 - twentyfour * sixteen * eps2 * dx * dy3 * rd2 * over8 - twelve**2 * dx * dy * eps2 * over8 * logTmp & !xyyy
                  - twentyfour * eight * dx * dy3 * eps2 * rd4 * over6 - two * eight**2 * dx * dy3 * eps2 * rd6 * over4 - ninetysix * dx * dy3 * eps2 * rd8 * over2 &
                  + twentyfour * sixteen * eps2 * dx * dy3 * logTmp * over10 + nine * eight * dx * dy * eps2 * rd4 * over4 + fortyeight * dx * dy * eps2 * rd6 * over2

      A1_sed(5) = twentyfour * twelve * dy2 * over6 * rd2 * eps2 - twentyfour * eps2 * over4 * rd2 - twentyfour * sixteen * eps2 * dy4 * over8 * rd2 & !yyyy
                  - twelve * eps2 * over2 * rd4 + twentyfour * eps2 * over6 * logTmp - twentyfour * twelve * dy2 * eps2 * over8 * logTmp &
                  + twentyfour * sixteen * eps2 * dy4 * over10 * logTmp + twelve**2 * dy2 * over4 * eps2 * rd4 + ninetysix * eps2 * dy2 * over2 * rd6 &
                  - twentyfour * eight * dy4 * eps2 * rd4 * over6 - two * eight**2 * dy4 * eps2 * rd6 * over4 - ninetysix * dy4 * eps2 * rd8 * over2

      A2_mono = -A1_mono + one + phi_mono + log(eps2)
      A2_dip(1:2) = -A1_dip(1:2) + phi_dip(1:2)
      A2_quad(1:3) = -A1_quad(1:3) + phi_quad(1:3)
      A2_oct(1) = -A1_oct(1) - Exx_quad(1)
      A2_oct(2) = -A1_oct(2) - Exy_quad(1)
      A2_oct(3) = -A1_oct(2) - Eyy_quad(1)
      A2_oct(4) = -A1_oct(3) - Eyy_quad(2)
      A2_sed(1) = -A1_sed(1) - Exxx_oct(1)
      A2_sed(2) = -A1_sed(2) - Exxy_oct(1)
      A2_sed(3) = -A1_sed(3) - Exyy_oct(1)
      A2_sed(4) = -A1_sed(4) - Eyyy_oct(1)
      A2_sed(5) = -A1_sed(5) - Eyyy_oct(2)

      A3_mono(1:3) = over * dist(1:3)
      A3x_dip(1) = over - dx2 * over * over2
      A3x_dip(2) = -dx * dy * over * over2
      A3x_dip(3) = zero
      A3y_dip(1) = -dx * dy * over * over2
      A3y_dip(2) = over - dy2 * over * over2
      A3y_dip(3) = zero
      A3xx_quad(1) = -three * over * over2 * dx * (one - dx2 * over2)
      A3xx_quad(2) = -over * over2 * dy * (one - three * dx2 * over2)   ! This is equal to A3_xy(1)
      A3xx_quad(3) = zero
      A3yy_quad(1) = -over * over2 * dx * (one - three * dy2 * over2)     ! This is equal to A3_xy(2)
      A3yy_quad(2) = -three * over * over2 * dy * (one - dy2 * over2)
      A3yy_quad(3) = zero
      A3xy_quad(1) = A3xx_quad(2)
      A3xy_quad(2) = A3yy_quad(1)
      A3xy_quad(3) = zero
      A3xxx_oct(1) = three * over**3 * (six * dx2 * over2 - one - five * dx4 * over4)
      A3xxx_oct(2) = three * over**5 * dx * dy * (three - five * dx2 * over2)
      A3xxx_oct(3) = zero

      A3xxy_oct(1) = three * over**5 * dx * dy * (three - five * dx2 * over2)
      A3xxy_oct(2) = over**3 * (two - fifteen * dx2 * dy2 * over4)
      A3xxy_oct(3) = zero

      A3xyy_oct(1) = over**3 * (two - fifteen * dx2 * dy2 * over4)
      A3xyy_oct(2) = three * over**5 * dx * dy * (three - five * dy2 * over2)
      A3xyy_oct(3) = zero

      A3yyy_oct(1) = three * over**5 * dx * dy * (three - five * dy2 * over2)
      A3yyy_oct(2) = three * over**3 * (six * dy2 * over2 - one - five * dy4 * over4)
      A3yyy_oct(3) = zero

      A3xxxx_sed(1) = fortyfive * dx * over5 - hundredfifty * dx3 * over7 + hundredfive * dx5 * over9
      A3xxxx_sed(2) = nine * dy * over5 - ninety * dx2 * dy * over7 + hundredfive * dx4 * dy * over9
      A3xxxx_sed(3) = zero

      A3xxxy_sed(1) = nine * dy * over5 - ninety * dx2 * dy * over7 + hundredfive * dx4 * dy * over9
      A3xxxy_sed(2) = nine * dx * over5 - fifteen * dx3 * over7 + hundredfive * dx3 * dy2 * over9 - fortyfive * dx * dy2 * over7
      A3xxxy_sed(3) = zero

      A3xxyy_sed(1) = nine * dx * over5 - fifteen * dx3 * over7 + hundredfive * dx3 * dy2 * over9 - fortyfive * dx * dy2 * over7
      A3xxyy_sed(2) = nine * dy * over5 - fifteen * dy3 * over7 + hundredfive * dx2 * dy3 * over9 - fortyfive * dx2 * dy * over7
      A3xxyy_sed(3) = zero

      A3xyyy_sed(1) = nine * dy * over5 - fifteen * dy3 * over7 + hundredfive * dx2 * dy3 * over9 - fortyfive * dx2 * dy * over7
      A3xyyy_sed(2) = nine * dx * over5 - ninety * dx * dy2 * over7 + hundredfive * dx * dy4 * over9
      A3xyyy_sed(3) = zero

      A3yyyy_sed(1) = nine * dx * over5 - ninety * dx * dy2 * over7 + hundredfive * dx * dy4 * over9
      A3yyyy_sed(2) = fortyfive * dy * over5 - hundredfifty * dy3 * over7 + hundredfive * dy5 * over9
      A3yyyy_sed(3) = zero

!      A4_mono(1:3)      = dist(1:3)
!      A4x_dip(1)        = one
!      A4x_dip(2)        = zero
!      A4x_dip(3)        = zero
!      A4y_dip(1)        = zero
!      A4y_dip(2)        = one
!      A4y_dip(3)        = zero

      phi = (t%charge * phi_mono) &! Monopole
            + (dot_product(t%dip(1:2), phi_dip(1:2))) &! Dipole
            + (dot_product(t%quad(1:2), phi_quad(1:2))) &! Quadrupole
            + (phi_quad(3) * t%xyquad)                                                                 ! Quadrupole
      phi = half * phi

      E(1:2) = (t%charge * E_mono(1:2)) &! Monopole
               + (t%dip(1) * Ex_dip(1:2) + t%dip(2) * Ey_dip(1:2)) &! Dipole
               + (t%quad(1) * Exx_quad(1:2) + t%quad(2) * Eyy_quad(1:2)) &! Quadrupole
               + (t%xyquad * Exy_quad(1:2))                                                               ! Quadrupole

      E(1:2) = half * E(1:2)

!      Ex(1:2)    = ( t%charge*Ex_dip(1:2)   )                                                              &! Monopole
!                  +( t%dip(1)*Exx_quad(1:2)  + t%dip(2)*Exy_quad(1:2)  )                                   &! Dipole
!                  +( t%quad(1)*Exxx_oct(1:2) + t%quad(2)*Exyy_oct(1:2) )                                   &! Quadrupole
!                  +( t%xyquad*Exxy_oct(1:2) )                                                              ! Quadrupole
!
!      Ex(1:2)    = half*Ex(1:2)
!
!
!      Ey(1:2)    = ( t%charge*Ey_dip(1:2)   )                                                              &! Monopole
!                  +( t%dip(1)*Exy_quad(1:2)  + t%dip(2)*Eyy_quad(1:2)  )                                   &! Dipole
!                  +( t%quad(1)*Exxy_oct(1:2) + t%quad(2)*Eyyy_oct(1:2) )                                   &! Quadrupole
!                  +( t%xyquad*Exyy_oct(1:2) )                                                               ! Quadrupole
!
!      Ey(1:2)    = half*Ey(1:2)

      J(1:3) = (t%monoj(1:3) * rho_mono) &! Monopole
               + (rho_dip(1) * t%dipjx(1:3) + rho_dip(2) * t%dipjy(1:3)) &! Dipole
               + (t%quadjx(1:3) * rho_quad(1) + t%quadjy(1:3) * rho_quad(2)) &! Quadrupole
               + (t%quadjxy(1:3) * rho_quad(3))                                                          ! Quadrupole

      J(1:3) = half * eps2 / pi * J(1:3)

      Jirr(1:3) = (t%monoj(1:3) * (three * half * Jir_mono - eps2 * rho_mono) &
                   + double_cross_product_left(E_mono(1:3), E_mono(1:3), t%monoj(1:3))) &! Monopole
                  + (t%dipjx(1:3) * (three * half * Jir_dip(1) - eps2 * rho_dip(1))) &! Dipole - x
                  + (double_cross_product_left(Ex_dip(1:3), E_mono(1:3), t%dipjx(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Ex_dip(1:3), t%dipjx(1:3))) &
                  + (t%dipjy(1:3) * (three * half * Jir_dip(2) - eps2 * rho_dip(2))) &! Dipole - y
                  + (double_cross_product_left(Ey_dip(1:3), E_mono(1:3), t%dipjy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Ey_dip(1:3), t%dipjy(1:3))) &
                  + (t%quadjx(1:3) * (three * half * Jir_quad(1) - eps2 * rho_quad(1))) &! Quadrupole - xx
                  + (double_cross_product_left(Exx_quad(1:3), E_mono(1:3), t%quadjx(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Exx_quad(1:3), t%quadjx(1:3))) &
                  + (two * double_cross_product_left(Ex_dip(1:3), Ex_dip(1:3), t%quadjx(1:3))) &
                  + (t%quadjy(1:3) * (three * half * Jir_quad(2) - eps2 * rho_quad(2))) &! Quadrupole - yy
                  + (double_cross_product_left(Eyy_quad(1:3), E_mono(1:3), t%quadjy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Eyy_quad(1:3), t%quadjy(1:3))) &
                  + (two * double_cross_product_left(Ey_dip(1:3), Ey_dip(1:3), t%quadjy(1:3))) &
                  + (t%quadjxy(1:3) * (three * half * Jir_quad(3) - eps2 * rho_quad(3))) &! Quadrupole - xy
                  + (double_cross_product_left(Exy_quad(1:3), E_mono(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(E_mono(1:3), Exy_quad(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(Ex_dip(1:3), Ey_dip(1:3), t%quadjxy(1:3))) &
                  + (double_cross_product_left(Ey_dip(1:3), Ex_dip(1:3), t%quadjxy(1:3)))

      Jirr(1:3) = half * oneoverpi * Jirr(1:3)

      Btmp(1:3) = cross_product(t%monoj(1:3), E_mono(1:3)) &! Monopole
                  + (cross_product(t%dipjx(1:3), Ex_dip(1:3))) &! Dipole x
                  + (cross_product(t%dipjy(1:3), Ey_dip(1:3))) &! Dipole y
                  + (cross_product(t%quadjx(1:3), Exx_quad(1:3))) &! Quadrupole xx
                  + (cross_product(t%quadjy(1:3), Eyy_quad(1:3))) &! Quadrupole yy
                  + (cross_product(t%quadjxy(1:3), Exy_quad(1:3)))                                          ! Quadrupole xy

      B(1:3) = half * Btmp(1:3)

      A(1:3) = half * t%monoj(1:3) * A2_mono &
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)) &! Monopole
               + half * t%dipjx(1:3) * A2_dip(1) &! Dipole - x
               + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3_mono(1:3)) &! I order  I term
               + A1_mono * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)) &! I order  II term
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)) &! I order  III term
               + half * t%dipjy(1:3) * A2_dip(2) &! Dipole - y
               + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3_mono(1:3)) &! I order  I term
               + A1_mono * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)) &! I order  II term
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)) &! I order  III term
               + half * t%quadjx(1:3) * A2_quad(1) &! Quadrupole - xx
               + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3_mono(1:3)) &!II order  I term
               + A1_dip(1) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3_mono(1:3)) &!II order  II-a term
               + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3x_dip(1:3)) &!II order  II-b term
               + A1_mono * double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)) &!II order  III-a term
               + two * A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3x_dip(1:3)) &!II order  III-b term
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &!II order  III-c term
               + half * t%quadjy(1:3) * A2_quad(2) &! Quadrupole - yy
               + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II       I term
               + A1_dip(2) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order  II-a term
               + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order  II-b term
               + A1_mono * double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order  III-a term
               + two * A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order  III-b term
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3yy_quad(1:3)) &! II order  III-c term
               + half * t%quadjxy(1:3) * A2_quad(3) &! Quadrupole - xy
               + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II - a term
               + A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II - b term
               + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II - c term
               + A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! III - a term
               + A1_mono * double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! III - b term
               + A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! III - c term
               + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! IV - a term
               + A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! IV - b term
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3xy_quad(1:3))     ! IV - c term

      A(3) = t%monoj(3) * phi_mono &! Monopole
             + t%dipjx(3) * phi_dip(1) + t%dipjy(3) * phi_dip(2) &! Dipole
             + t%quadjx(3) * phi_quad(1) + t%quadjy(3) * phi_quad(2) &! Quadrupole
             + t%quadjxy(3) * phi_quad(3)                                                               ! Quadrupole

      A(1:3) = half * A(1:3)

      Ax(1:3) = half * t%monoj(1:3) * A2_dip(1) &! Monopole
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)) &
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)) &
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3x_dip(1:3)) &
                + half * t%dipjx(1:3) * A2_quad(1) &! Dipole - x I term
                + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3_mono(1:3)) &! I order - II term
                + two * A1_dip(1) * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)) &! I order - III - a term
                + two * A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)) &! I order - III - b term
                + A1_mono * double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A3_mono(1:3)) &! I order - IV  - a term
                + two * A1_mono * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)) &! I order - IV  - b term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3xx_quad(1:3)) &! I order - IV  - c term
                + half * t%dipjy(1:3) * A2_quad(3) &! Dipole - y
                + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3_mono(1:3)) &! II term
                + A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A3_mono(1:3)) &
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3x_dip(1:3)) &
                + A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)) &! III term
                + A1_mono * double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)) &
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A3x_dip(1:3)) &
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)) &! IV term
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)) &
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3xy_quad(1:3)) &
                + half * t%quadjx(1:3) * A2_oct(1) &! Quadrupole - xx
                + A1_oct(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - I a term
                + three * A1_quad(1) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - II a term
                + three * A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - II b term
                + three * A1_dip(1) * double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - III a term
                + six * A1_dip(1) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - III b term
                + three * A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &! II order - III c term
                + A1_mono * double_cross_product_left(A3xxx_oct(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - IV a term
                + three * A1_mono * double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - IV b term
                + three * A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &! II order - IV c term
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &! II order - IV d term
                + half * t%quadjy(1:3) * A2_oct(3) &! Quadrupole - yy
                + A1_oct(3) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - I a term
                + A1_quad(2) * double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - II a term
                + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3x_dip(1:3)) &! II order - II b term
                + two * A1_quad(3) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - III a term
                + two * A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - III b term
                + two * A1_dip(2) * double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - IV  a term
                + two * A1_dip(2) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3x_dip(1:3)) &! II order - IV  b term
                + two * A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - IV  c term
                + two * A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3xy_quad(1:3)) &! II order - IV  d term
                + A1_dip(1) * double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - V  a term
                + two * A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - V  b term
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3yy_quad(1:3)) &! II order - V  c term
                + A1_mono * double_cross_product_left(A3xyy_oct(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - VI  a term
                + A1_mono * double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3x_dip(1:3)) &! II order - VI  b term
                + two * A1_mono * double_cross_product_left(A3xy_quad(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - VI  c term
                + two * A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3xy_quad(1:3)) &! II order - VI  d term
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjy(1:3), A3yy_quad(1:3)) &! II order - VI  e term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3xyy_oct(1:3)) &! II order - VI  f term
                + half * t%quadjxy(1:3) * A2_oct(2) &! Quadrupole - xy
                + A1_oct(2) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - I a term
                + A1_quad(1) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - II a term
                + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - II b term
                + two * A1_quad(3) * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - III a term
                + two * A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - III b term
                + two * A1_dip(1) * double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - IV a term
                + two * A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - IV b term
                + A1_dip(2) * double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - IV a term
                + two * A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - IV b term
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3xx_quad(1:3)) &! II order - IV c term
                + A1_mono * double_cross_product_left(A3xxy_oct(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - V a term
                + A1_mono * double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - V b term
                + two * A1_mono * double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - V c term
                + two * A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3xy_quad(1:3)) &! II order - V d term
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3xx_quad(1:3)) &! II order - V e term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3xxy_oct(1:3))     ! II order - V f term

      Ax(3) = t%monoj(3) * phi_dip(1) &! Monopole
              + t%dipjx(3) * phi_quad(1) + t%dipjy(3) * phi_quad(3) &! Dipole
              - t%quadjx(3) * Exx_quad(1) - t%quadjy(3) * Eyy_quad(1) &! Quadrupole
              - t%quadjxy(3) * Exx_quad(2)                                                                 ! Quadrupole

      Ay(1:3) = half * t%monoj(1:3) * A2_dip(2) &
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)) &! Monopole
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)) &! Monopole
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)) &! Monopole
                + half * t%dipjx(1:3) * A2_quad(3) &! Dipole - x I term
                + A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3_mono(1:3)) &! II term
                + A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A3_mono(1:3)) &! II term
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3y_dip(1:3)) &! II term
                + A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)) &! III term
                + A1_mono * double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A3_mono(1:3)) &! III term
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%dipjx(1:3), A3y_dip(1:3)) &! III term
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)) &! IV term
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)) &! IV term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjx(1:3), A3xy_quad(1:3)) &! IV term
                + half * t%dipjy(1:3) * A2_quad(2) &! Dipole - y
                + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3_mono(1:3)) &! I order - I term
                + two * A1_dip(2) * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)) &! I order - II - a term
                + two * A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)) &! I order - II - b term
                + A1_mono * double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)) &! I order - III - a term
                + two * A1_mono * double_cross_product_left(A3y_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)) &! I order - III - b term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%dipjy(1:3), A3yy_quad(1:3)) &! I order - III - c term
                + half * t%quadjx(1:3) * A2_oct(2) &! Quadrupole - xx
                + A1_oct(2) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - I term
                + A1_quad(1) * double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - II - a term
                + A1_quad(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3y_dip(1:3)) &! II order - II - b term
                + two * A1_quad(3) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - III- a term
                + two * A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - III- b term
                + two * A1_dip(1) * double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - IV- a term
                + two * A1_dip(1) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3y_dip(1:3)) &! II order - IV- b term
                + two * A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - IV- c term
                + two * A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3xy_quad(1:3)) &! II order - IV- d term
                + A1_dip(2) * double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - V- a term
                + two * A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - V- b term
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &! II order - V- c term
                + A1_mono * double_cross_product_left(A3xxy_oct(1:3), t%quadjx(1:3), A3_mono(1:3)) &! II order - VI- a term
                + A1_mono * double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3y_dip(1:3)) &! II order - VI- b term
                + two * A1_mono * double_cross_product_left(A3xy_quad(1:3), t%quadjx(1:3), A3x_dip(1:3)) &! II order - VI- c term
                + two * A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjx(1:3), A3xy_quad(1:3)) &! II order - VI- d term
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjx(1:3), A3xx_quad(1:3)) &! II order - VI- e term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjx(1:3), A3xxy_oct(1:3)) &! II order - VI- f term
                + half * t%quadjy(1:3) * A2_oct(4) &! Quadrupole - yy
                + A1_oct(4) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - I- a term
                + three * A1_quad(2) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - II- a term
                + three * A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - II- b term
                + three * A1_dip(2) * double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - III- a term
                + six * A1_dip(2) * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - III- b term
                + three * A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3yy_quad(1:3)) &! II order - III- c term
                + A1_mono * double_cross_product_left(A3yyy_oct(1:3), t%quadjy(1:3), A3_mono(1:3)) &! II order - IV- a term
                + three * A1_mono * double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3y_dip(1:3)) &! II order - IV- b term
                + three * A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjy(1:3), A3yy_quad(1:3)) &! II order - IV- c term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjy(1:3), A3yyy_oct(1:3)) &! II order - IV- d term
                + half * t%quadjxy(1:3) * A2_oct(3) &! Quadrupole - xy
                + A1_oct(3) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - I- a term
                + A1_quad(2) * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - II- a term
                + A1_quad(2) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - II- b term
                + two * A1_quad(3) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - III- a term
                + two * A1_quad(3) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - III- b term
                + two * A1_dip(2) * double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - III- c term
                + two * A1_dip(2) * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - III- d term
                + two * A1_dip(2) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - III- e term
                + two * A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3xy_quad(1:3)) &! II order - III- f term
                + A1_dip(1) * double_cross_product_left(A3xx_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - IV- a term
                + two * A1_dip(1) * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - IV- b term
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3yy_quad(1:3)) &! II order - IV- c term
                + A1_mono * double_cross_product_left(A3xyy_oct(1:3), t%quadjxy(1:3), A3_mono(1:3)) &! II order - V- a term
                + A1_mono * double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A3x_dip(1:3)) &! II order - V- a term
                + two * A1_mono * double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3y_dip(1:3)) &! II order - V- a term
                + two * A1_mono * double_cross_product_left(A3y_dip(1:3), t%quadjxy(1:3), A3xy_quad(1:3)) &! II order - V- a term
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%quadjxy(1:3), A3yy_quad(1:3)) &! II order - V- a term
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%quadjxy(1:3), A3xyy_oct(1:3))     ! II order - V- a term

!
      Ay(3) = t%monoj(3) * phi_dip(2) &! Monopole
              + t%dipjx(3) * phi_quad(3) + t%dipjy(3) * phi_quad(2) &! Dipole
              - t%quadjx(3) * Exx_quad(2) - t%quadjy(3) * Eyy_quad(2) &! Quadrupole
              - t%quadjxy(3) * Eyy_quad(1)                                                                 ! Quadrupole

      Ax(1:3) = half * Ax(1:3)
      Ay(1:3) = half * Ay(1:3)

!
!      Axx(1:3)   =  half*t%monoj(1:3)*A2_quad(1)                                                            &! Monopole
!                 + A1_quad(1)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - I- a term
!
!              +two*A1_dip(1)*double_cross_product_left( A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - II- a term
!              +two*A1_dip(1)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3x_dip(1:3)  )         &! 0 order - II- b term
!
!                 + A1_mono*double_cross_product_left( A3xx_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - III- a term
!              +two*A1_mono*double_cross_product_left( A3x_dip(1:3)  , t%monoj(1:3), A3x_dip(1:3)  )         &! 0 order - III- b term
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3xx_quad(1:3))         &! 0 order - III- c term
!
!                 + half*  t%dipjx(1:3)*A2_oct(1)                                                            &! Dipole - x I term
!                 + A1_oct(1)* double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! I order - I- a term
!
!            +three*A1_quad(1)*double_cross_product_left( A3x_dip(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! I order - II- a term
!            +three*A1_quad(1)*double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - II- b term
!
!            +three*A1_dip(1)*double_cross_product_left( A3xx_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - III- a term
!            +six*  A1_dip(1)*double_cross_product_left( A3x_dip(1:3)  , t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - III- b term
!            +three*A1_dip(1)*double_cross_product_left( A3_mono(1:3)  , t%dipjx(1:3), A3xx_quad(1:3))       &! I order - III- c term
!
!                 + A1_mono* double_cross_product_left(  A3xxx_oct(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - IV- a term
!            +three*A1_mono* double_cross_product_left(  A3xx_quad(1:3), t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - IV- b term
!            +three*A1_mono* double_cross_product_left(  A3x_dip(1:3)  , t%dipjx(1:3), A3xx_quad(1:3))       &! I order - IV- c term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjx(1:3), A3xxx_oct(1:3))       &! I order - IV- d term
!
!                 + half*  t%dipjy(1:3)*A2_oct(2)                                                            &! Dipole - y
!                 + A1_oct(2) *double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! I order - I- a term
!
!                 + A1_quad(1)* double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - II- a term
!                 + A1_quad(1)* double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - II- b term
!
!              +two*A1_quad(3)*double_cross_product_left( A3x_dip(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! I order - III- a term
!              +two*A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3x_dip(1:3)  )       &! I order - III- b term
!
!              +two*A1_dip(1)* double_cross_product_left(A3xy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - IV- a term
!              +two*A1_dip(1)* double_cross_product_left(  A3y_dip(1:3), t%dipjy(1:3), A3x_dip(1:3)  )       &! I order - IV- b term
!              +two*A1_dip(1)* double_cross_product_left(  A3x_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - IV- b term
!              +two*A1_dip(1)* double_cross_product_left(  A3_mono(1:3), t%dipjy(1:3), A3xy_quad(1:3))       &! I order - IV- c term
!
!                 + A1_dip(2)*double_cross_product_left( A3xx_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - V- a term
!              +two*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)  , t%dipjy(1:3), A3x_dip(1:3)  )       &! I order - V- b term
!                 + A1_dip(2)*double_cross_product_left( A3_mono(1:3)  , t%dipjy(1:3), A3xx_quad(1:3))       &! I order - V- c term
!
!                 + A1_mono* double_cross_product_left(  A3xxy_oct(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3xx_quad(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - VI- b term
!              +two*A1_mono* double_cross_product_left(  A3xy_quad(1:3), t%dipjy(1:3), A3x_dip(1:3)  )       &! I order - VI- c term
!              +two*A1_mono* double_cross_product_left(  A3x_dip(1:3)  , t%dipjy(1:3), A3xy_quad(1:3))       &! I order - VI- d term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)  , t%dipjy(1:3), A3xx_quad(1:3))       &! I order - VI- e term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjy(1:3), A3xxy_oct(1:3))       &! I order - VI- f term
!
!                 + half*  t%quadjx(1:3)*A2_sed(1)                                                           &! Quadrupole - xx
!                 + A1_sed(1) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!             +four*A1_oct(1) *double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - II- a term
!             +four*A1_oct(1) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - II- b term
!
!              +six*A1_quad(1)*double_cross_product_left( A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - III- a term
!           +twelve*A1_quad(1)*double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - III- b term
!              +six*A1_quad(1)*double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - III- c term
!
!             +four*A1_dip(1)*double_cross_product_left( A3xxx_oct(1:3) , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!           +twelve*A1_dip(1)*double_cross_product_left( A3xx_quad(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - IV- b term
!           +twelve*A1_dip(1)*double_cross_product_left( A3x_dip(1:3)   , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - IV- c term
!             +four*A1_dip(1)*double_cross_product_left( A3_mono(1:3)   , t%quadjx(1:3), A3xxx_oct(1:3)  )   &! II order - IV- d term
!
!                 + A1_mono* double_cross_product_left(  A3xxxx_sed(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - V- a term
!             +four*A1_mono* double_cross_product_left(  A3xxx_oct(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - V- b term
!              +six*A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - V- c term
!             +four*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjx(1:3), A3xxx_oct(1:3)  )   &! II order - V- d term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjx(1:3), A3xxxx_sed(1:3))    &! II order - V- e term
!
!                 + half*  t%quadjy(1:3)*A2_sed(3)                                                           &! Quadrupole - yy
!                 + A1_sed(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!              +two*A1_oct(2) *double_cross_product_left( A3y_dip(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - II- a term
!              +two*A1_oct(2) *double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - II- b term
!
!                 + A1_quad(1)* double_cross_product_left(A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - III- a term
!              +two*A1_quad(1)* double_cross_product_left(A3y_dip(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - III- b term
!                 + A1_quad(1)* double_cross_product_left(A3_mono(1:3)  , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - III- c term
!
!              +two*A1_oct(3)*double_cross_product_left(  A3x_dip(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!              +two*A1_oct(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - IV- b term
!
!             +four*A1_quad(3)*double_cross_product_left( A3xy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - V- a term
!             +four*A1_quad(3)*double_cross_product_left( A3y_dip(1:3)  , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - V- b term
!             +four*A1_quad(3)*double_cross_product_left(A3x_dip(1:3)   , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - V- c term
!             +four*A1_quad(3)*double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - V- d term
!
!              +two*A1_dip(1)* double_cross_product_left(A3xyy_oct(1:3) , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VI- a term
!              +two*A1_dip(1)* double_cross_product_left(A3yy_quad(1:3) , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - VI- b term
!             +four*A1_dip(1)* double_cross_product_left(A3xy_quad(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - VI- c term
!             +four*A1_dip(1)* double_cross_product_left(A3y_dip(1:3)   , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - VI- d term
!              +two*A1_dip(1)* double_cross_product_left(A3x_dip(1:3)   , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - VI- e term
!              +two*A1_dip(1)* double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3xyy_oct(1:3))     &! II order - VI- f term
!
!                 + A1_quad(2)*double_cross_product_left( A3xx_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VII- a term
!              +two*A1_quad(2)*double_cross_product_left( A3x_dip(1:3)  , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - VII- b term
!                 + A1_quad(2)*double_cross_product_left( A3_mono(1:3)  , t%quadjy(1:3), A3xx_quad(1:3))     &! II order - VII- c term
!
!              +two*A1_dip(2)*double_cross_product_left( A3xxy_oct(1:3) , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VIII- a term
!              +two*A1_dip(2)*double_cross_product_left( A3xx_quad(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - VIII- b term
!             +four*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)   , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - VIII- c term
!             +four*A1_dip(2)*double_cross_product_left( A3xy_quad(1:3) , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - VIII- d term
!              +two*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjy(1:3), A3xx_quad(1:3))     &! II order - VIII- e term
!              +two*A1_dip(2)* double_cross_product_left(  A3_mono(1:3) , t%quadjy(1:3), A3xxy_oct(1:3))     &! II order - VIII- f term
!
!
!                 + A1_mono* double_cross_product_left(  A3xxyy_sed(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - IX- a term
!              +two*A1_mono* double_cross_product_left(  A3xxy_oct(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - IX- b term
!                 + A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - IX- c term
!
!              +two*A1_mono* double_cross_product_left(  A3xyy_oct(1:3) , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - IX- d term
!             +four*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - IX- e term
!              +two*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjy(1:3), A3xyy_oct(1:3))     &! II order - IX- f term
!
!                 + A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjy(1:3), A3xx_quad(1:3))     &! II order - IX- g term
!              +two*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjy(1:3), A3xxy_oct(1:3))     &! II order - IX- h term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjy(1:3), A3xxyy_sed(1:3))    &! II order - IX- i term
!
!
!
!                 + half*  t%quadjxy(1:3)*A2_sed(2)                                                          &! Quadrupole - xy
!
!                 + A2_sed(2)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - I- a term
!
!                 + A1_oct(1) * double_cross_product_left( A3y_dip(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - II- a term
!                 + A1_oct(1) * double_cross_product_left( A3_mono(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - II- b term
!
!            +three*A1_oct(2)*double_cross_product_left(  A3x_dip(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - III- a term
!            +three*A1_oct(2)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - III- b term
!
!            +three*A1_quad(1)* double_cross_product_left(A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - IV- a term
!            +three*A1_quad(1)* double_cross_product_left( A3x_dip(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - IV- b term
!            +three*A1_quad(1)* double_cross_product_left( A3y_dip(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - IV- c term
!            +three*A1_quad(1)* double_cross_product_left( A3_mono(1:3) , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - IV- d term
!
!            +three*A1_quad(3)*double_cross_product_left( A3xx_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - V- a term
!              +six*A1_quad(3)*double_cross_product_left( A3x_dip(1:3)  , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - V- b term
!            +three*A1_quad(3)*double_cross_product_left( A3_mono(1:3)  , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - V- c term
!
!            +three*A1_dip(1)* double_cross_product_left(  A3xxy_oct(1:3),t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VI- a term
!            +three*A1_dip(1)* double_cross_product_left(A3xx_quad(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VI- b term
!              +six*A1_dip(1)* double_cross_product_left(  A3xy_quad(1:3),t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VI- c term
!              +six*A1_dip(1)* double_cross_product_left(  A3x_dip(1:3),t%quadjxy(1:3), A3xy_quad(1:3)  )    &! II order - VI- d term
!            +three*A1_dip(1)* double_cross_product_left(  A3y_dip(1:3) , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - VI- e term
!            +three*A1_dip(1)* double_cross_product_left(  A3_mono(1:3) , t%quadjxy(1:3), A3xxy_oct(1:3))    &! II order - VI- f term
!
!                 + A1_dip(2)*double_cross_product_left( A3xxx_oct(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VII- a term
!            +three*A1_dip(2)*double_cross_product_left( A3xx_quad(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VII- b term
!            +three*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)   , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - VII- c term
!                 + A1_dip(2)*double_cross_product_left( A3_mono(1:3)   , t%quadjxy(1:3), A3xxx_oct(1:3))    &! II order - VII- d term
!
!                 + A1_mono* double_cross_product_left(  A3xxxy_sed(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VII- a term
!                 + A1_mono* double_cross_product_left(  A3xxx_oct(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VII- b term
!            +three*A1_mono* double_cross_product_left(  A3xxy_oct(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VII- c term
!            +three*A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - VII- d term
!            +three*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - VII- e term
!            +three*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjxy(1:3), A3xxy_oct(1:3))    &! II order - VII- f term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjxy(1:3), A3xxx_oct(1:3))    &! II order - VII- g term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjxy(1:3), A3xxxy_sed(1:3))    ! II order - VII- h term
!
!
!      Axx(3)     =   t%monoj(3)*phi_quad(1)                                                                  &! Monopole
!                  -  t%dipjx(3)*Exx_quad(1)   - t%dipjy(3)*Exx_quad(2)                                       &! Dipole
!                  -  t%quadjx(3)*Exxx_oct(1)  - t%quadjy(3)*Exyy_oct(1)                                      &! Quadrupole
!                  -  t%quadjxy(3)*Exxy_oct(1)                                                                 ! Quadrupole
!
!
!
!      Ayy(1:3)   =  half*t%monoj(1:3)*A2_quad(2)                                                            &! Monopole
!                 + A1_quad(2)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - I- a term
!
!              +two*A1_dip(2)*double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - II- a term
!              +two*A1_dip(2)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)  )         &! 0 order - II- b term
!
!                 + A1_mono*double_cross_product_left( A3yy_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! 0 order - II- a term
!              +two*A1_mono*double_cross_product_left( A3y_dip(1:3)  , t%monoj(1:3), A3y_dip(1:3)  )         &! 0 order - II- a term
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3yy_quad(1:3))         &! 0 order - II- a term
!
!                 + half*  t%dipjx(1:3)*A2_oct(3)                                                            &! Dipole - x I term
!                 + A1_oct(3)* double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! I order - I- a term
!
!                 + A1_quad(2)* double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - II- a term
!                 + A1_quad(2)* double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - II- b term
!
!              +two*A1_quad(3)*double_cross_product_left( A3y_dip(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! I order - III- a term
!              +two*A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3y_dip(1:3)  )       &! I order - III- b term
!
!              +two*A1_dip(2)* double_cross_product_left(A3xy_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - IV- a term
!              +two*A1_dip(2)* double_cross_product_left(A3x_dip(1:3)  , t%dipjx(1:3), A3y_dip(1:3)  )       &! I order - IV- b term
!              +two*A1_dip(2)* double_cross_product_left(  A3y_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - IV- c term
!              +two*A1_dip(2)* double_cross_product_left(  A3_mono(1:3), t%dipjx(1:3), A3xy_quad(1:3))       &! I order - IV- d term
!
!
!                 + A1_dip(1)*double_cross_product_left( A3yy_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - V- a term
!              +two*A1_dip(1)*double_cross_product_left( A3y_dip(1:3)  , t%dipjx(1:3), A3y_dip(1:3)  )       &! I order - V- b term
!                 + A1_dip(1)*double_cross_product_left( A3_mono(1:3)  , t%dipjx(1:3), A3yy_quad(1:3))       &! I order - V- c term
!
!
!                 + A1_mono* double_cross_product_left(  A3xyy_oct(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! I order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3yy_quad(1:3), t%dipjx(1:3), A3x_dip(1:3)  )       &! I order - VI- b term
!              +two*A1_mono* double_cross_product_left(  A3xy_quad(1:3), t%dipjx(1:3), A3y_dip(1:3)  )       &! I order - VI- c term
!              +two*A1_mono* double_cross_product_left(  A3y_dip(1:3)  , t%dipjx(1:3), A3xy_quad(1:3))       &! I order - VI- d term
!                 + A1_mono* double_cross_product_left(  A3x_dip(1:3)  , t%dipjx(1:3), A3yy_quad(1:3))       &! I order - VI- e term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjx(1:3), A3xyy_oct(1:3))       &! I order - VI- f term
!
!
!                 + half*  t%dipjy(1:3)*A2_oct(4)                                                            &! Dipole - y
!                 + A1_oct(4) *double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! I order - I- a term
!
!            +three*A1_quad(2)*double_cross_product_left( A3y_dip(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! I order - II- a term
!            +three*A1_quad(2)*double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - II- b term
!
!            +three*A1_dip(2)*double_cross_product_left( A3yy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - III- a term
!            +six*  A1_dip(2)*double_cross_product_left( A3y_dip(1:3)  , t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - III- b term
!            +three*A1_dip(2)*double_cross_product_left( A3_mono(1:3)  , t%dipjy(1:3), A3yy_quad(1:3))       &! I order - III- c term
!
!
!                 + A1_mono* double_cross_product_left(  A3yyy_oct(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! I order - IV- a term
!            +three*A1_mono* double_cross_product_left(  A3yy_quad(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! I order - IV- b term
!            +three*A1_mono* double_cross_product_left(  A3y_dip(1:3)  , t%dipjy(1:3), A3yy_quad(1:3))       &! I order - IV- c term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjy(1:3), A3yyy_oct(1:3))       &! I order - IV- d term
!
!
!                 + half*  t%quadjx(1:3)*A2_sed(3)                                                           &! Quadrupole - xx
!                 + A1_sed(3) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!              +two*A1_oct(2) *double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - II- a term
!              +two*A1_oct(2) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - II- b term
!
!              +two*A1_quad(1)*double_cross_product_left( A3yy_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - III- a term
!             +four*A1_quad(1)*double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - III- b term
!              +two*A1_quad(1)*double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3yy_quad(1:3))     &
!
!              +two*A1_oct(3) *double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!              +two*A1_oct(3) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - IV- b term
!
!              +two*A1_oct(2) *double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - V- a term
!              +two*A1_oct(2) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - V- b term
!
!             +four*A1_quad(3)*double_cross_product_left( A3xy_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - VI- b term
!             +four*A1_quad(3)*double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - VI- b term
!             +four*A1_quad(3)*double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - VI- b term
!             +four*A1_quad(3)*double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - VI- b term
!
!              +two*A1_dip(1)*double_cross_product_left( A3xyy_oct(1:3) , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - VII- a term
!              +two*A1_dip(1)* double_cross_product_left(A3yy_quad(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - VII- b term
!             +four*A1_dip(1)*double_cross_product_left( A3xy_quad(1:3) , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - VII- c term
!             +four*A1_dip(1)*double_cross_product_left( A3y_dip(1:3)   , t%quadjx(1:3), A3xy_quad(1:3)  )   &! II order - VII- d term
!              +two*A1_dip(1)*double_cross_product_left( A3x_dip(1:3)   , t%quadjx(1:3), A3yy_quad(1:3))     &! II order - VII- e term
!              +two*A1_dip(1)*double_cross_product_left( A3_mono(1:3)   , t%quadjx(1:3), A3xyy_oct(1:3)  )   &! II order - VII- f term
!
!                 + A1_quad(2)* double_cross_product_left(A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - VIII- a term
!              +two*A1_quad(2)* double_cross_product_left( A3x_dip(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - VIII- b term
!                 + A1_quad(2)* double_cross_product_left( A3_mono(1:3) , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - VIII- c term
!
!              +two*A1_dip(2)* double_cross_product_left( A3xxy_oct(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - IX- a term
!              +two*A1_dip(2)* double_cross_product_left( A3xx_quad(1:3), t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - IX- b term
!             +four*A1_dip(2)* double_cross_product_left(  A3xy_quad(1:3), t%quadjx(1:3), A3x_dip(1:3) )     &! II order - IX- c term
!             +four*A1_dip(2)* double_cross_product_left(  A3x_dip(1:3) , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - IX- d term
!              +two*A1_dip(2)* double_cross_product_left(  A3y_dip(1:3) , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - IX- e term
!              +two*A1_dip(2)* double_cross_product_left(  A3_mono(1:3) , t%quadjx(1:3), A3xxy_oct(1:3))     &! II order - IX- f term
!
!                 + A1_mono* double_cross_product_left(  A3xxyy_sed(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - X- a term
!              +two*A1_mono* double_cross_product_left(  A3xxy_oct(1:3) , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - X- b term
!                 + A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjx(1:3), A3yy_quad(1:3))     &! II order - X- c term
!              +two*A1_mono* double_cross_product_left(  A3xyy_oct(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - X- d term
!             +four*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - X- e term
!              +two*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjx(1:3), A3xyy_oct(1:3))     &! II order - X- f term
!
!                 + A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - XI- a term
!              +two*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjx(1:3), A3xxy_oct(1:3))     &! II order - XI- b term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjx(1:3), A3xxyy_sed(1:3))    &! II order - XI- c term
!
!
!
!                 + half*  t%quadjy(1:3)*A2_sed(5)                                                           &! Quadrupole - yy
!                 + A1_sed(5)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!             +four*A1_oct(4)*double_cross_product_left(  A3y_dip(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - II- a term
!             +four*A1_oct(4)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - II- b term
!
!              +six*A1_quad(2)*double_cross_product_left( A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - III- a term
!           +twelve*A1_quad(2)*double_cross_product_left( A3y_dip(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - III- a term
!              +six*A1_quad(2)*double_cross_product_left( A3_mono(1:3)  , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - III- a term
!
!             +four*A1_dip(2)*double_cross_product_left( A3yyy_oct(1:3) , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!           +twelve*A1_dip(2)*double_cross_product_left( A3yy_quad(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - IV- a term
!           +twelve*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - IV- a term
!             +four*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - IV- a term
!             +four*A1_dip(2)*double_cross_product_left( A3_mono(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - IV- a term
!
!
!                 + A1_mono* double_cross_product_left(  A3yyyy_sed(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - V- a term
!             +four*A1_mono* double_cross_product_left(  A3yyy_oct(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - V- a term
!             +six* A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - V- a term
!             +four*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - V- a term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - V- a term
!
!                 + A1_dip(2)* double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjy(1:3), A3yyyy_sed(1:3))    &! II order - VI- a term
!
!                 + half*  t%quadjxy(1:3)*A2_sed(4)                                                          &! Quadrupole - xy
!                 + A2_sed(4)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - I- a term
!
!                 + A1_oct(4) *double_cross_product_left( A3x_dip(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - II- a term
!                 + A1_oct(4) *double_cross_product_left(A3_mono(1:3)   , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - II- a term
!
!            +three*A1_oct(3)*double_cross_product_left(  A3y_dip(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - III- a term
!            +three*A1_oct(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - III- b term
!
!            +three*A1_quad(2)*double_cross_product_left( A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - IV- a term
!            +three*A1_quad(2)*double_cross_product_left( A3x_dip(1:3)  , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - IV- b term
!            +three*A1_quad(2)*double_cross_product_left(A3y_dip(1:3)   , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - IV- c term
!            +three*A1_quad(2)*double_cross_product_left(A3_mono(1:3)   , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - IV- d term
!
!            +three*A1_quad(3)*double_cross_product_left( A3yy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - V- a term
!            +six*  A1_quad(3)*double_cross_product_left( A3y_dip(1:3)  , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - V- b term
!            +three*A1_quad(3)*double_cross_product_left( A3_mono(1:3)  , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - V- c term
!
!            +three*A1_dip(2)*double_cross_product_left( A3xyy_oct(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VI- a term
!            +three*A1_dip(2)*double_cross_product_left( A3yy_quad(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VI- b term
!              +six*A1_dip(2)*double_cross_product_left( A3xy_quad(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VI- c term
!              +six*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - VI- d term
!            +three*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)   , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - VI- e term
!            +three*A1_dip(2)*double_cross_product_left( A3_mono(1:3)   , t%quadjxy(1:3), A3xyy_oct(1:3))    &! II order - VI- f term
!
!                 + A1_dip(1)* double_cross_product_left( A3yyy_oct(1:3) ,t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VII- a term
!            +three*A1_dip(1)* double_cross_product_left(  A3yy_quad(1:3),t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VII- b term
!            +three*A1_dip(1)* double_cross_product_left(  A3y_dip(1:3)  , t%quadjxy(1:3), A3yy_quad(1:3))   &! II order - VII- c term
!                 + A1_dip(1)* double_cross_product_left(  A3_mono(1:3) , t%quadjxy(1:3), A3yyy_oct(1:3))    &! II order - VII- d term
!
!
!
!                 + A1_mono* double_cross_product_left(  A3xyyy_sed(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VIII- a term
!                 + A1_mono* double_cross_product_left(  A3yyy_oct(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VIII- b term
!            +three*A1_mono* double_cross_product_left(  A3xyy_oct(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VIII- c term
!            +three*A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - VIII- d term
!            +three*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - VIII- e term
!            +three*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjxy(1:3), A3xyy_oct(1:3))    &! II order - VIII- f term
!                 + A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjxy(1:3), A3yyy_oct(1:3))    &! II order - VIII- g term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjxy(1:3), A3xyyy_sed(1:3))    ! II order - VIII- h term
!
!
!      Ayy(3)     =   t%monoj(3)*phi_quad(2)                                                                  &! Monopole
!                  -  t%dipjx(3)*Eyy_quad(1)   - t%dipjy(3)*Eyy_quad(2)                                       &! Dipole
!                  -  t%quadjx(3)*Exxy_oct(2)  - t%quadjy(3)*Eyyy_oct(2)                                      &! Quadrupole
!                  -  t%quadjxy(3)*Exyy_oct(2)                                                                 ! Quadrupole
!
!
!      Axy(1:3)   =  half*t%monoj(1:3)*A2_quad(3)                                                            &! Monopole
!                 + A1_quad(3)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! II order - I- a term
!
!                 + A1_dip(1)*double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! II order - II- a term
!                 + A1_dip(1)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)  )         &! II order - II- b term
!
!                 + A1_dip(2)*double_cross_product_left(A3_mono(1:3) , t%monoj(1:3), A3x_dip(1:3)  )         &! II order - III- a term
!                 + A1_dip(2)*double_cross_product_left(A3x_dip(1:3) , t%monoj(1:3), A3_mono(1:3)  )         &! II order - III- b term
!
!                 + A1_mono*double_cross_product_left( A3xy_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &! II order - IV- a term
!                 + A1_mono*double_cross_product_left( A3x_dip(1:3)  , t%monoj(1:3), A3y_dip(1:3)  )         &! II order - IV- b term
!                 + A1_mono*double_cross_product_left( A3y_dip(1:3)  , t%monoj(1:3), A3x_dip(1:3)  )         &! II order - IV- c term
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3xy_quad(1:3))         &! II order - IV- d term
!
!                 + half*  t%dipjx(1:3)*A2_oct(2)                                                            &! Dipole - x I term
!
!                 + A1_oct(2)* double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! II order - I- a term
!
!                 + A1_quad(1)*double_cross_product_left( A3y_dip(1:3) , t%dipjx(1:3), A3_mono(1:3)  )       &! II order - II- a term
!                 + A1_quad(1)*double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3y_dip(1:3)  )       &! II order - II- b term
!
!              +two*A1_quad(3)*double_cross_product_left(A3x_dip(1:3)  , t%dipjx(1:3), A3_mono(1:3)  )       &! II order - III- a term
!              +two*A1_quad(3)*double_cross_product_left(A3_mono(1:3)  , t%dipjx(1:3), A3x_dip(1:3)  )       &! II order - III- b term
!
!              +two*A1_dip(1)*double_cross_product_left( A3xy_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! II order - IV- a term
!              +two*A1_dip(1)*double_cross_product_left( A3x_dip(1:3)  , t%dipjx(1:3), A3y_dip(1:3)  )       &! II order - IV- b term
!              +two*A1_dip(1)*double_cross_product_left( A3y_dip(1:3)  , t%dipjx(1:3), A3x_dip(1:3)  )       &! II order - IV- c term
!              +two*A1_dip(1)*double_cross_product_left( A3_mono(1:3)  , t%dipjx(1:3), A3xy_quad(1:3))       &! II order - IV- d term
!
!                 + A1_dip(2)* double_cross_product_left(A3xx_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! II order - V- a term
!              +two*A1_dip(2)* double_cross_product_left(A3x_dip(1:3)  , t%dipjx(1:3), A3x_dip(1:3)  )       &! II order - V- b term
!                 + A1_dip(2)* double_cross_product_left(A3_mono(1:3)  , t%dipjx(1:3), A3xx_quad(1:3))       &! II order - V- c term
!
!                 + A1_mono* double_cross_product_left(  A3xxy_oct(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &! II order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3xx_quad(1:3), t%dipjx(1:3), A3y_dip(1:3)  )       &! II order - VI- b term
!              +two*A1_mono* double_cross_product_left(  A3xy_quad(1:3), t%dipjx(1:3), A3x_dip(1:3)  )       &! II order - VI- c term
!              +two*A1_mono* double_cross_product_left(  A3x_dip(1:3)  , t%dipjx(1:3), A3xy_quad(1:3))       &! II order - VI- d term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)  , t%dipjx(1:3), A3xx_quad(1:3))       &! II order - VI- e term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjx(1:3), A3xxy_oct(1:3))       &! II order - VI- f term
!
!
!                 + half*  t%dipjy(1:3)*A2_oct(3)                                                            &! Dipole - y
!                 + A1_oct(3) *double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! I order - I- a term
!
!                 + A1_quad(2)*double_cross_product_left(A3x_dip(1:3)  , t%dipjy(1:3), A3_mono(1:3)  )       &! II order - II- a term
!                 + A1_quad(2)*double_cross_product_left(A3_mono(1:3)  , t%dipjy(1:3), A3x_dip(1:3)  )       &! II order - II- a term
!
!              +two*A1_quad(3)*double_cross_product_left( A3y_dip(1:3) , t%dipjy(1:3), A3_mono(1:3)  )       &! II order - III- a term
!              +two*A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3y_dip(1:3)  )       &! II order - III- b term
!
!
!              +two*A1_dip(2)*double_cross_product_left( A3xy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! II order - IV- a term
!              +two*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)  , t%dipjy(1:3), A3y_dip(1:3)  )       &! II order - IV- b term
!              +two*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)  , t%dipjy(1:3), A3x_dip(1:3)  )       &! II order - IV- c term
!              +two*A1_dip(2)*double_cross_product_left( A3_mono(1:3)  , t%dipjy(1:3), A3xy_quad(1:3))       &! II order - IV- D term
!
!                 + A1_dip(1)* double_cross_product_left(A3yy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! II order - V- a term
!              +two*A1_dip(1)* double_cross_product_left(  A3y_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! II order - V- b term
!                 + A1_dip(1)* double_cross_product_left(  A3_mono(1:3), t%dipjy(1:3), A3yy_quad(1:3))       &! II order - V- c term
!
!                 + A1_mono* double_cross_product_left(  A3xyy_oct(1:3), t%dipjy(1:3), A3_mono(1:3)  )       &! II order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3yy_quad(1:3), t%dipjy(1:3), A3x_dip(1:3)  )       &! II order - VI- b term
!              +two*A1_mono* double_cross_product_left(  A3xy_quad(1:3), t%dipjy(1:3), A3y_dip(1:3)  )       &! II order - VI- c term
!              +two*A1_mono* double_cross_product_left(  A3y_dip(1:3)  , t%dipjy(1:3), A3xy_quad(1:3))       &! II order - VI- d term
!                 + A1_mono* double_cross_product_left(  A3x_dip(1:3)  , t%dipjy(1:3), A3yy_quad(1:3))       &! II order - VI- e term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)  , t%dipjy(1:3), A3xyy_oct(1:3))       &! II order - VI- f term
!
!
!                 + half*  t%quadjx(1:3)*A2_sed(2)                                                           &! Quadrupole - xx
!                 + A1_sed(2) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!                 + A1_oct(1) *double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - II- a term
!                 + A1_oct(1) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - II- b term
!
!            +three*A1_oct(2) *double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - III- a term
!            +three*A1_oct(2) *double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - III- b term
!
!            +three*A1_quad(1)*double_cross_product_left( A3xy_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - III- a term
!            +three*A1_quad(1)*double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - III- b term
!            +three*A1_quad(1)*double_cross_product_left( A3y_dip(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - III- c term
!            +three*A1_quad(1)*double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - III- d term
!
!            +three*A1_quad(3)*double_cross_product_left( A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!            +six*  A1_quad(3)*double_cross_product_left( A3x_dip(1:3)  , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - IV- b term
!            +three*A1_quad(3)*double_cross_product_left( A3_mono(1:3)  , t%quadjx(1:3), A3xx_quad(1:3))     &! II order - IV- c term
!
!                 + A1_dip(1)*double_cross_product_left( A3xxy_oct(1:3) , t%quadjx(1:3), A3_mono(1:3)  )     &! II order - V- a term
!                 + A1_dip(1)*double_cross_product_left( A3xx_quad(1:3) , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - V- b term
!              +two*A1_dip(1)*double_cross_product_left( A3xy_quad(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - V- c term
!              +two*A1_dip(1)*double_cross_product_left( A3x_dip(1:3)   , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - V- d term
!                 + A1_dip(1)*double_cross_product_left( A3y_dip(1:3)   , t%quadjx(1:3), A3xx_quad(1:3)  )   &! II order - V- e term
!                 + A1_dip(1)*double_cross_product_left( A3_mono(1:3)   , t%quadjx(1:3), A3xxy_oct(1:3)  )   &! II order - V- f term
!
!
!                 + A1_mono* double_cross_product_left(  A3xxxy_sed(1:3), t%quadjx(1:3), A3_mono(1:3)  )     &! II order - VI- a term
!                 + A1_mono* double_cross_product_left(  A3xxx_oct(1:3) , t%quadjx(1:3), A3y_dip(1:3)  )     &! II order - VI- b term
!            +three*A1_mono* double_cross_product_left(  A3xxy_oct(1:3) , t%quadjx(1:3), A3x_dip(1:3)  )     &! II order - VI- c term
!            +three*A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjx(1:3), A3xy_quad(1:3))     &! II order - VI- d term
!            +three*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjx(1:3), A3xx_quad(1:3)  )   &! II order - VI- e term
!            +three*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjx(1:3), A3xxy_oct(1:3)  )   &! II order - VI- f term
!                 + A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjx(1:3), A3xxx_oct(1:3))     &! II order - VI- g term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjx(1:3), A3xxxy_sed(1:3))    &! II order - VI- h term
!
!                 + half*  t%quadjy(1:3)*A2_sed(4)                                                           &! Quadrupole - yy
!                 + A1_sed(4)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - I- a term
!
!                 + A1_oct(4) *double_cross_product_left( A3x_dip(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - II- a term
!                 + A1_oct(4) *double_cross_product_left( A3_mono(1:3)  , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - II- b term
!
!            +three*A1_oct(3)*double_cross_product_left(  A3y_dip(1:3)  , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - III- a term
!            +three*A1_oct(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - III- b term
!
!
!            +three*A1_quad(2)*double_cross_product_left( A3xy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - IV- a term
!            +three*A1_quad(2)*double_cross_product_left( A3x_dip(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - IV- b term
!            +three*A1_quad(2)*double_cross_product_left( A3y_dip(1:3)  , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - IV- c term
!            +three*A1_quad(2)*double_cross_product_left( A3_mono(1:3)  , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - IV- d term
!
!            +three*A1_quad(3)*double_cross_product_left( A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - V- a term
!            +six*  A1_quad(3)*double_cross_product_left( A3y_dip(1:3)  , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - V- b term
!            +three*A1_quad(3)*double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - V- c term
!
!            +three*A1_dip(2)*double_cross_product_left( A3xyy_oct(1:3) , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VI- a term
!            +three*A1_dip(2)*double_cross_product_left( A3yy_quad(1:3) , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - VI- b term
!            +six*  A1_dip(2)*double_cross_product_left( A3xy_quad(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - VI- c term
!            +six*  A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - VI- d term
!            +three*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)   , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - VI- e term
!            +three*A1_dip(2)* double_cross_product_left(  A3_mono(1:3) , t%quadjy(1:3), A3xyy_oct(1:3))     &! II order - VI- f term
!
!                 + A1_dip(1)* double_cross_product_left(A3yyy_oct(1:3) , t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VII- a term
!              +two*A1_dip(1)* double_cross_product_left(A3yy_quad(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - VII- b term
!                 + A1_dip(1)* double_cross_product_left(A3_mono(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - VII- c term
!
!                 + A1_mono* double_cross_product_left(  A3xyyy_sed(1:3), t%quadjy(1:3), A3_mono(1:3)  )     &! II order - VIII- a term
!                 + A1_mono* double_cross_product_left(  A3yyy_oct(1:3) , t%quadjy(1:3), A3x_dip(1:3)  )     &! II order - VIII- b term
!            +three*A1_mono* double_cross_product_left(  A3xyy_oct(1:3) , t%quadjy(1:3), A3y_dip(1:3)  )     &! II order - VIII- c term
!            +three*A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjy(1:3), A3xy_quad(1:3))     &! II order - VIII- d term
!            +three*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjy(1:3), A3yy_quad(1:3))     &! II order - VIII- e term
!            +three*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjy(1:3), A3xyy_oct(1:3))     &! II order - VIII- f term
!                 + A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjy(1:3), A3yyy_oct(1:3))     &! II order - VIII- g term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjy(1:3), A3xyyy_sed(1:3))    &! II order - VIII- h term
!
!
!                 + half*  t%quadjxy(1:3)*A2_sed(3)                                                          &! Quadrupole - xy
!                 + A2_sed(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - I- a term
!
!              +two*A1_oct(2)*double_cross_product_left(  A3y_dip(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - II- a term
!              +two*A1_oct(2)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - II- b term
!
!                 + A1_quad(1)* double_cross_product_left(A3yy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - III- a term
!              +two*A1_quad(1)* double_cross_product_left( A3y_dip(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - III- b term
!                 + A1_quad(1)* double_cross_product_left( A3_mono(1:3) , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - III- c term
!
!              +two*A1_oct(3)*double_cross_product_left(  A3x_dip(1:3)  , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - IV- a term
!              +two*A1_oct(3)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - IV- b term
!
!             +four*A1_quad(3)*double_cross_product_left( A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - V- a term
!             +four*A1_quad(3)*double_cross_product_left( A3x_dip(1:3)  , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - V- b term
!             +four*A1_quad(3)*double_cross_product_left( A3y_dip(1:3)  , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - V- c term
!             +four*A1_quad(3)*double_cross_product_left( A3_mono(1:3)  , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - V- d term
!
!              +two*A1_dip(1)* double_cross_product_left(  A3xyy_oct(1:3),t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VI- a term
!              +two*A1_dip(1)* double_cross_product_left( A3yy_quad(1:3), t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VI- b term
!             +four*A1_dip(1)* double_cross_product_left(  A3xy_quad(1:3),t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VI- c term
!             +four*A1_dip(1)* double_cross_product_left(  A3y_dip(1:3) , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - VI- d term
!              +two*A1_dip(1)*double_cross_product_left(   A3x_dip(1:3) , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - VI- e term
!              +two*A1_dip(1)*double_cross_product_left(  A3_mono(1:3)  , t%quadjxy(1:3), A3xyy_oct(1:3))    &! II order - VI- f term
!
!
!                 + A1_quad(2)*double_cross_product_left(A3xx_quad(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VII- a term
!              +two*A1_quad(2)*double_cross_product_left(A3x_dip(1:3)   , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VII- b term
!                 + A1_quad(2)*double_cross_product_left( A3_mono(1:3)  , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - VII- c term
!
!
!              +two*A1_dip(2)*double_cross_product_left( A3xxy_oct(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - VIII- a term
!              +two*A1_dip(2)*double_cross_product_left( A3xx_quad(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - VIII- b term
!             +four*A1_dip(2)*double_cross_product_left( A3xy_quad(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - VIII- c term
!             +four*A1_dip(2)*double_cross_product_left( A3x_dip(1:3)   , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - VIII- d term
!              +two*A1_dip(2)*double_cross_product_left( A3y_dip(1:3)   , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - VIII- e term
!              +two*A1_dip(2)*double_cross_product_left( A3_mono(1:3)   , t%quadjxy(1:3), A3xxy_oct(1:3))    &! II order - VIII- f term
!
!                 + A1_mono* double_cross_product_left(  A3xxyy_sed(1:3), t%quadjxy(1:3), A3_mono(1:3)  )    &! II order - IX- a term
!              +two*A1_mono* double_cross_product_left(  A3xxy_oct(1:3) , t%quadjxy(1:3), A3y_dip(1:3)  )    &! II order - IX- b term
!                 + A1_mono* double_cross_product_left(  A3xx_quad(1:3) , t%quadjxy(1:3), A3yy_quad(1:3))    &! II order - IX- c term
!              +two*A1_mono* double_cross_product_left(  A3xyy_oct(1:3) , t%quadjxy(1:3), A3x_dip(1:3)  )    &! II order - IX- d term
!             +four*A1_mono* double_cross_product_left(  A3xy_quad(1:3) , t%quadjxy(1:3), A3xy_quad(1:3))    &! II order - IX- e term
!              +two*A1_mono* double_cross_product_left(  A3x_dip(1:3)   , t%quadjxy(1:3), A3xyy_oct(1:3))    &! II order - IX- f term
!                 + A1_mono* double_cross_product_left(  A3yy_quad(1:3) , t%quadjxy(1:3), A3xx_quad(1:3))    &! II order - IX- g term
!              +two*A1_mono* double_cross_product_left(  A3y_dip(1:3)   , t%quadjxy(1:3), A3xxy_oct(1:3))    &! II order - IX- h term
!                 + A1_mono* double_cross_product_left(  A3_mono(1:3)   , t%quadjxy(1:3), A3xxyy_sed(1:3))    ! II order - IX- i term
!
!      Axy(3)     =   t%monoj(3)*phi_quad(3)                                                                  &! Monopole
!                  -  t%dipjx(3)*Exx_quad(2)   - t%dipjy(3)*Eyy_quad(1)                                       &! Dipole
!                  -  t%quadjx(3)*Exxx_oct(2)  - t%quadjy(3)*Eyyy_oct(1)                                      &! Quadrupole
!                  -  t%quadjxy(3)*Exxy_oct(2)                                                                 ! Quadrupole
!
!
!      Axx(1:3)   = half*Axx(1:3)
!      Ayy(1:3)   = half*Ayy(1:3)
!      Axy(1:3)   = half*Axy(1:3)

      A(1:2) = -A(1:2)
      Ax(1:2) = -Ax(1:2)
      Ay(1:2) = -Ay(1:2)
!      Axx(1:2)   = -Axx(1:2)
!      Axy(1:2)   = -Axy(1:2)
!      Ayy(1:2)   = -Ayy(1:2)

   end subroutine calc_force_darwin_2D3V

   subroutine calc_force_darwin_3D(t, d, d2, eps2, phi, exyz, Axyz, Jxyz, Jirrxyz, Bxyz, dAx, dAy, dAz)
      use module_tool, only: cross_product, double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  ::  d(1:3), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(1:3), Axyz(1:3), Jxyz(1:3), Jirrxyz(1:3), Bxyz(1:3), phi, dAx(1:3), dAy(1:3), dAz(1:3)

      real(kind_physics) :: dx, dy, dz, r2, r, rd, rEps, rd2, rd3, rd4, rd5, rd7, rd8, rd9, over, over2, over3, over4, over5, over7, over6, over8, over9, &
                            dx2, dy2, dz2, dx3, dy3, dz3, rd6, rd10, rd12, rd14, &
                            rho_mono, rho_dip, rho_quad(1:4), &                                     !density/current density
                            phi_mono, phi_dip, phi_quad(1:4), &                                     !potential
                            E_mono(1:3), E_dipx(1:3), E_dipy(1:3), E_dipz(1:3), E_quadxx(1:3), E_quadyy(1:3), E_quadzz(1:3), & !Eirr
                            E_quadxy(1:3), E_quadyz(1:3), E_quadzx(1:3), &
                            fmono, fpmono, A1_mono, A1_dip, A1_quad, A1_oct, &                         !A
                            A2_mono, A2_dip, A2_quad, A2_oct, f_dip(1:3), f_quad(1:6), f_oct(1:10), &
                            fp_dip(1:3), fp_quad(1:6), fp_oct(1:10), Ix(1:3), Iy(1:3), Iz(1:3), fppmono, fpppmono, fppppmono, &
                            DxDxJ(1:3), DxDxJ_x(1:3), I_xxDxJ_x(1:3), DxI_xxJ_x(1:3), &
                            DxDxJ_y(1:3), I_yxDxJ_y(1:3), DxI_yxJ_y(1:3), &
                            DxDxJ_z(1:3), I_zxDxJ_z(1:3), DxI_zxJ_z(1:3), &
                            DxDxJ_xx(1:3), I_xxDxJ_xx(1:3), DxI_xxJ_xx(1:3), &
                            DxDxJ_yy(1:3), I_yxDxJ_yy(1:3), DxI_yxJ_yy(1:3), &
                            I_xxI_xxJ_xx(1:3), I_yxI_xxJ_yy(1:3)

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r2 = dx**2 + dy**2 + dz**2
      r = sqrt(r2)
      rEps = sqrt(d2)
      rd2 = one / d2 ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd = sqrt(rd2)
      rd3 = rd**3
      rd4 = rd2**2
      rd5 = rd**5
      rd6 = rd2**3
      rd7 = rd**7
      rd8 = rd2**4
      rd9 = rd**9
      rd10 = rd2**5
      rd12 = rd2**6
      rd14 = rd2**7
      over = one / r
      over2 = over * over
      over3 = over2 * over
      over4 = over3 * over
      over5 = over4 * over
      over6 = over5 * over
      over7 = over6 * over
      over8 = over7 * over
      over9 = over8 * over

      dx2 = dx * dx
      dy2 = dy * dy
      dz2 = dz * dz
      dx3 = dx * dx2
      dy3 = dy * dy2
      dz3 = dz * dz2

      Ix(1:3) = (/one, zero, zero/)
      Iy(1:3) = (/zero, one, zero/)
      Iz(1:3) = (/zero, zero, one/)

      rho_mono = rd5
      rho_dip = -five * rd7
      rho_quad(1) = five * rd7 * (seven * dx2 * rd2 - one) !xx
      rho_quad(2) = five * rd7 * (seven * dy2 * rd2 - one) !yy
      rho_quad(3) = five * rd7 * (seven * dz2 * rd2 - one) !zz

      rho_quad(4) = five * seven * rd9         !mixed terms

      phi_mono = rd
      phi_dip = -rd3
      phi_quad(1) = three * dx2 * rd5 - rd3  !xx
      phi_quad(2) = three * dy2 * rd5 - rd3  !yy
      phi_quad(3) = three * dz2 * rd5 - rd3  !zz

      phi_quad(4) = three * rd5         !mixed terms

      E_mono = rd3 * d(1:3)
      E_dipx = (/rd3 - three * rd5 * dx2, -three * dx * dy * rd5, -three * dx * dz * rd5/)
      E_dipy = (/-three * dx * dy * rd5, rd3 - three * rd5 * dy2, -three * dy * dz * rd5/)
      E_dipz = (/-three * dx * dz * rd5, -three * dy * dz * rd5, rd3 - three * rd5 * dz2/)
      E_quadxx = (/fifteen * rd7 * dx3 - nine * dx * rd5, fifteen * rd7 * dx2 * dy - three * dy * rd5, fifteen * rd7 * dx2 * dz - three * dz * rd5/)
      E_quadyy = (/fifteen * rd7 * dx * dy2 - three * dx * rd5, fifteen * rd7 * dy3 - nine * dy * rd5, fifteen * rd7 * dy2 * dz - three * dz * rd5/)
      E_quadzz = (/fifteen * rd7 * dx * dz2 - three * dx * rd5, fifteen * rd7 * dz2 * dy - three * dy * rd5, fifteen * rd7 * dz3 - nine * dz * rd5/)
      E_quadxy = (/fifteen * rd7 * dx2 * dy - three * dy * rd5, fifteen * rd7 * dx * dy2 - three * dx * rd5, fifteen * dx * dy * dz * rd7/)
      E_quadyz = (/E_quadxy(3), fifteen * rd7 * dy2 * dz - three * dz * rd5, fifteen * rd7 * dz2 * dy - three * dy * rd5/)
      E_quadzx = (/fifteen * rd7 * dx2 * dz - three * dz * rd5, E_quadxy(3), fifteen * rd7 * dz2 * dx - three * dx * rd5/)

      Jxyz = rho_mono * t%monoj(1:3) + rho_dip * (dx * t%dipjx(1:3) + dy * t%dipjy(1:3) + dz * t%dipjz(1:3)) &
             + rho_quad(1) * t%quadjx(1:3) + rho_quad(2) * t%quadjy(1:3) + rho_quad(3) * t%quadjz(1:3) &
             + rho_quad(4) * (dx * dy * t%quadjxy(1:3) + dz * dy * t%quadjyz(1:3) + dx * dz * t%quadjzx(1:3))

      Jxyz = three / four / pi * eps2 * Jxyz

      phi = (t%charge * phi_mono) &! Monopole
            + phi_dip * (dx * t%dip(1) + dy * t%dip(2) + dz * t%dip(3)) &! Dipole
            + phi_quad(1) * t%quad(1) + phi_quad(2) * t%quad(2) + phi_quad(3) * t%quad(3) &! Quadrupole
            + phi_quad(4) * (dx * dy * t%xyquad + dz * dy * t%yzquad + dx * dz * t%zxquad)

      Exyz = t%charge * E_mono &! Monopole
             + t%dip(1) * E_dipx + t%dip(2) * E_dipy + t%dip(3) * E_dipz &! Dipole
             + t%quad(1) * E_quadxx + t%quad(2) * E_quadyy + t%quad(3) * E_quadzz &! Quadrupole
             + t%xyquad * E_quadxy + t%yzquad * E_quadyz + t%zxquad * E_quadzx

      Bxyz = cross_product(t%monoj(1:3), E_mono(1:3)) &! Monopole
             + cross_product(t%dipjx(1:3), E_dipx(1:3)) &! Dipole x
             + cross_product(t%dipjy(1:3), E_dipy(1:3)) &! Dipole y
             + cross_product(t%dipjz(1:3), E_dipz(1:3)) &! Dipole z
             + cross_product(t%quadjx(1:3), E_quadxx(1:3)) &! Quadrupole xx
             + cross_product(t%quadjy(1:3), E_quadyy(1:3)) &! Quadrupole yy
             + cross_product(t%quadjz(1:3), E_quadzz(1:3)) &! Quadrupole zz
             + cross_product(t%quadjxy(1:3), E_quadxy(1:3)) &! Quadrupole xy
             + cross_product(t%quadjyz(1:3), E_quadyz(1:3)) &! Quadrupole yz
             + cross_product(t%quadjzx(1:3), E_quadzx(1:3))                                ! Quadrupole zx

      A1_mono = over2
      A1_dip = -two * over4                  ! first  order derivative
      A1_quad = eight * over6                   ! second order derivative
      A1_oct = -fortyeight * over8              ! third order derivative

      A2_mono = A1_mono * rd
      A2_dip = -over4 * rd3 * (two * d2 + r2)                                                      ! first  order derivative
      A2_quad = four * (two * d2 + r2) * over6 * rd3 - six * over4 * rd3 + three * (two * d2 + r2) * over4 * rd5       ! second order derivative
      A2_oct = fortyeight * over6 * rd3 + six**2 * over4 * rd5 - twentyfour * (two * d2 + r2) * over6 * rd5 &
               - fifteen * (two * d2 + r2) * over4 * rd7 - twentyfour * (two * d2 + r2) * over8 * rd3    ! third order derivative

      fmono = -half * rEps * over2 + half * eps2 * over3 * log(rEps + r)
      fpmono = -three * over2 * fmono - over2 * rd
      fppmono = -(three * A1_mono * fpmono + three * fmono * A1_dip + A2_dip)
      fpppmono = -(three * A1_mono * fppmono + six * A1_dip * fpmono + three * A1_quad * fmono + A2_quad)
      fppppmono = -(three * A1_mono * fpppmono + nine * A1_dip * fppmono + nine * A1_quad * fpmono + three * A1_oct * fmono + A2_oct)

      f_dip(1:3) = fpmono * (/dx, dy, dz/)
      fp_dip(1:3) = fppmono * (/dx, dy, dz/)

      f_quad(1:3) = fppmono * (/dx2, dy2, dz2/) + fpmono                            !diagonal
      f_quad(4:6) = fppmono * (/dx * dy, dy * dz, dz * dx/)                               !mixed terms

      fp_quad(1:3) = fpppmono * (/dx2, dy2, dz2/) + fppmono                         !diagonal
      fp_quad(4:6) = fpppmono * (/dx * dy, dy * dz, dz * dx/)                              !mixed terms

      f_oct(1:3) = fpppmono * (/dx3, dy3, dz3/) + three * fppmono * (/dx, dy, dz/)
      f_oct(4:5) = fpppmono * (/dx2 * dy, dx2 * dz/) + fppmono * (/dy, dz/)
      f_oct(6:7) = fpppmono * (/dy2 * dx, dy2 * dz/) + fppmono * (/dx, dz/)
      f_oct(8:9) = fpppmono * (/dz2 * dx, dz2 * dy/) + fppmono * (/dx, dy/)
      f_oct(10) = fpppmono * dx * dy * dz

      fp_oct(1:3) = fppppmono * (/dx3, dy3, dz3/) + three * fpppmono * (/dx, dy, dz/)
      fp_oct(4:5) = fppppmono * (/dx2 * dy, dx2 * dz/) + fpppmono * (/dy, dz/)
      fp_oct(6:7) = fppppmono * (/dy2 * dx, dy2 * dz/) + fpppmono * (/dx, dz/)
      fp_oct(8:9) = fppppmono * (/dz2 * dx, dz2 * dy/) + fpppmono * (/dx, dy/)
      fp_oct(10) = fppppmono * dx * dy * dz

      DxDxJ = double_cross_product_left(d(1:3), d(1:3), t%monoj(1:3))
      DxDxJ_x = double_cross_product_left(d(1:3), d(1:3), t%dipjx(1:3))
      I_xxDxJ_x = double_cross_product_left(Ix(1:3), d(1:3), t%dipjx(1:3))
      DxI_xxJ_x = double_cross_product_left(d(1:3), Ix(1:3), t%dipjx(1:3))
      DxDxJ_y = double_cross_product_left(d(1:3), d(1:3), t%dipjy(1:3))
      I_yxDxJ_y = double_cross_product_left(Iy(1:3), d(1:3), t%dipjy(1:3))
      DxI_yxJ_y = double_cross_product_left(d(1:3), Iy(1:3), t%dipjy(1:3))
      DxDxJ_z = double_cross_product_left(d(1:3), d(1:3), t%dipjz(1:3))
      I_zxDxJ_z = double_cross_product_left(Iz(1:3), d(1:3), t%dipjz(1:3))
      DxI_zxJ_z = double_cross_product_left(d(1:3), Iz(1:3), t%dipjz(1:3))
      DxDxJ_xx = double_cross_product_left(d(1:3), d(1:3), t%quadjx(1:3))
      I_xxDxJ_xx = double_cross_product_left(Ix(1:3), d(1:3), t%quadjx(1:3))
      DxI_xxJ_xx = double_cross_product_left(d(1:3), Ix(1:3), t%quadjx(1:3))
      I_xxI_xxJ_xx = double_cross_product_left(Ix(1:3), Ix(1:3), t%quadjx(1:3))
      DxDxJ_yy = double_cross_product_left(d(1:3), d(1:3), t%quadjy(1:3))
      I_yxDxJ_yy = double_cross_product_left(Iy(1:3), d(1:3), t%quadjy(1:3))
      I_yxI_xxJ_yy = double_cross_product_left(d(1:3), Iy(1:3), t%quadjy(1:3))

      Axyz = fpmono * DxDxJ - two * fmono * t%monoj(1:3) &  ! Monopole
             + fp_dip(1) * DxDxJ_x &  ! Dipole - x a
             + fpmono * I_xxDxJ_x &  !b
             + fpmono * DxI_xxJ_x &  !c
             - two * f_dip(1) * t%dipjx(1:3) &  !d
             + fp_dip(2) * DxDxJ_y &  ! Dipole - y a
             + fpmono * I_yxDxJ_y &  !b
             + fpmono * DxI_yxJ_y &  !c
             - two * f_dip(2) * t%dipjy(1:3) &  !d
             + fp_dip(3) * DxDxJ_z &  ! Dipole - z a
             + fpmono * I_zxDxJ_z &  !b
             + fpmono * DxI_zxJ_z &  !c
             - two * f_dip(3) * t%dipjz(1:3) &  !d
             + fp_quad(1) * DxDxJ_xx & ! Quadrupole - xx a
             + fp_dip(1) * I_xxDxJ_xx * two & ! b
             + fp_dip(1) * DxI_xxJ_xx * two & ! c
             + fpmono * I_xxI_xxJ_xx * two & ! d
             - two * f_quad(1) * t%quadjx(1:3) & ! e
             + fp_quad(2) * DxDxJ_yy & ! Quadrupole - yy a
             + fp_dip(2) * I_yxDxJ_yy * two & ! b
             + fp_dip(2) * I_yxI_xxJ_yy * two & ! c
             + fpmono * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two & ! d
             - two * f_quad(2) * t%quadjy(1:3) & ! e
             + fp_quad(3) * double_cross_product_left(d(1:3), d(1:3), t%quadjz(1:3)) & ! Quadrupole - zz a
             + fp_dip(3) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjz(1:3)) * two & ! b
             + fp_dip(3) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! c
             + fpmono * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! d
             - two * f_quad(3) * t%quadjz(1:3) & ! e
             + fp_quad(4) * double_cross_product_left(d(1:3), d(1:3), t%quadjxy(1:3)) & ! Quadrupole - xy a
             + fp_dip(1) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjxy(1:3)) & ! b
             + fp_dip(1) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
             + fp_dip(2) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjxy(1:3)) & ! d
             + fp_dip(2) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
             + fpmono * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! f
             + fpmono * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! g
             - two * f_quad(4) * t%quadjxy(1:3) & ! h
             + fp_quad(5) * double_cross_product_left(d(1:3), d(1:3), t%quadjyz(1:3)) & ! Quadrupole - yz a
             + fp_dip(3) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjyz(1:3)) & ! b
             + fp_dip(3) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
             + fp_dip(2) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjyz(1:3)) & ! d
             + fp_dip(2) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
             + fpmono * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! f
             + fpmono * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! g
             - two * f_quad(5) * t%quadjyz(1:3) & ! h
             + fp_quad(6) * double_cross_product_left(d(1:3), d(1:3), t%quadjzx(1:3)) & ! Quadrupole - zx a
             + fp_dip(3) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjzx(1:3)) & ! b
             + fp_dip(3) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
             + fp_dip(1) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjzx(1:3)) & ! d
             + fp_dip(1) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
             + fpmono * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! f
             + fpmono * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! g
             - two * f_quad(6) * t%quadjzx(1:3)                                                                           ! h

      dAx = fp_dip(1) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Ix(1:3), d(1:3), t%monoj(1:3)) &
            + fpmono * double_cross_product_left(d(1:3), Ix(1:3), t%monoj(1:3)) &
            - two * f_dip(1) * t%monoj(1:3) &
            + fp_quad(1) * DxDxJ_x &  ! Dipole - x a
            + fp_dip(1) * I_xxDxJ_x &
            + fp_dip(1) * DxI_xxJ_x * three &
            !                      +fp_dip(1)*double_cross_product_left( Ix(1:3) , d(1:3), t%dipjx(1:3)  )                               &  !b
            + fpmono * double_cross_product_left(Ix(1:3), Ix(1:3), t%dipjx(1:3)) &
            !                     +fp_dip(1)*double_cross_product_left( d(1:3) , Ix(1:3), t%dipjx(1:3)  )                                &  !c
            + fpmono * double_cross_product_left(Ix(1:3), Ix(1:3), t%dipjx(1:3)) &  !c
            - two * f_quad(1) * t%dipjx(1:3) &  !d
            + fp_quad(4) * DxDxJ_y &  ! Dipole - y a
            + fp_dip(2) * double_cross_product_left(Ix(1:3), d(1:3), t%dipjy(1:3)) &  !
            + fp_dip(2) * double_cross_product_left(d(1:3), Ix(1:3), t%dipjy(1:3)) &  !
            + fp_dip(1) * I_yxDxJ_y &  !b
            + fpmono * double_cross_product_left(Iy(1:3), Ix(1:3), t%dipjy(1:3)) &  !b
            + fp_dip(1) * DxI_yxJ_y &  !c
            + fpmono * double_cross_product_left(Ix(1:3), Iy(1:3), t%dipjy(1:3)) &  !c
            - two * f_quad(4) * t%dipjy(1:3) &  !d
            + fp_quad(6) * DxDxJ_z &  ! Dipole - z a
            + fp_dip(3) * double_cross_product_left(Ix(1:3), d(1:3), t%dipjz(1:3)) &
            + fp_dip(3) * double_cross_product_left(d(1:3), Ix(1:3), t%dipjz(1:3)) &
            + fp_dip(1) * I_zxDxJ_z &  !b
            + fpmono * double_cross_product_left(Iz(1:3), Ix(1:3), t%dipjz(1:3)) &  !
            + fp_dip(1) * DxI_zxJ_z &  !c
            + fpmono * double_cross_product_left(Ix(1:3), Iz(1:3), t%dipjz(1:3)) &  !
            - two * f_quad(6) * t%dipjz(1:3) &  !d
            + fp_oct(1) * DxDxJ_xx & ! Quadrupole - xx a
            !                     +fp_quad(1)*double_cross_product_left( Ix(1:3) , d(1:3), t%quadjx(1:3)  )                               & !
            !                     +fp_quad(1)*double_cross_product_left( d(1:3) , Ix(1:3), t%quadjx(1:3)  )                               & !
            + fp_quad(1) * I_xxDxJ_xx * three &! b
            + fp_dip(1) * I_xxI_xxJ_xx * six &! b
            + fp_quad(1) * DxI_xxJ_xx * three &! c
            !                     +fp_dip(1)*double_cross_product_left( Ix(1:3) , Ix(1:3), t%quadjx(1:3)  )*two                           &! c
            !                     +fp_dip(1)*double_cross_product_left( Ix(1:3) , Ix(1:3), t%quadjx(1:3)  )*two                           &! d
            - two * f_oct(1) * t%quadjx(1:3) &! e
            + fp_oct(6) * DxDxJ_yy & ! Quadrupole - yy a
            + fp_quad(2) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjy(1:3)) & !
            + fp_quad(2) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjy(1:3)) & !
            + fp_quad(4) * I_yxDxJ_yy * two &! b
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjy(1:3)) * two &! b
            + fp_quad(4) * I_yxI_xxJ_yy * two &! c
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjy(1:3)) * two &! c
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two &! d
            - two * f_oct(6) * t%quadjy(1:3) &! e
            + fp_oct(8) * double_cross_product_left(d(1:3), d(1:3), t%quadjz(1:3)) & ! Quadrupole - zz a
            + fp_quad(3) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjz(1:3)) & !
            + fp_quad(3) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjz(1:3)) & !
            + fp_quad(6) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjz(1:3)) * two &! b
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjz(1:3)) * two &! b
            + fp_quad(6) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjz(1:3)) * two &! c
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjz(1:3)) * two &! c
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two &! d
            - two * f_oct(8) * t%quadjz(1:3) &! e
            + fp_oct(4) * double_cross_product_left(d(1:3), d(1:3), t%quadjxy(1:3)) & ! Quadrupole - xy a
            + fp_quad(4) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjxy(1:3)) & !
            + fp_quad(4) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjxy(1:3)) & !
            + fp_quad(1) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjxy(1:3)) & ! b
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! b
            + fp_quad(1) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
            + fp_quad(4) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjxy(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Ix(1:3), t%quadjxy(1:3)) & ! d
            + fp_quad(4) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! f
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! g
            - two * f_oct(4) * t%quadjxy(1:3) & ! h
            + fp_oct(10) * double_cross_product_left(d(1:3), d(1:3), t%quadjyz(1:3)) & ! Quadrupole - yz a
            + fp_quad(5) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjyz(1:3)) &
            + fp_quad(5) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjyz(1:3)) &
            + fp_quad(6) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjyz(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjyz(1:3)) & ! b
            + fp_quad(6) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_quad(4) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjyz(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjyz(1:3)) & ! d
            + fp_quad(4) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! f
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! g
            - two * f_oct(10) * t%quadjyz(1:3) & ! h
            + fp_oct(5) * double_cross_product_left(d(1:3), d(1:3), t%quadjzx(1:3)) & ! Quadrupole - zx a
            + fp_quad(6) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjzx(1:3)) & !
            + fp_quad(6) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjzx(1:3)) & !
            + fp_quad(6) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjzx(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Ix(1:3), t%quadjzx(1:3)) & ! b
            + fp_quad(6) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_quad(1) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjzx(1:3)) & ! d
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! d
            + fp_quad(1) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! f
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! g
            - two * f_oct(5) * t%quadjzx(1:3)                                                                           ! h

      dAy = fp_dip(2) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Iy(1:3), d(1:3), t%monoj(1:3)) &  ! Monopole
            + fpmono * double_cross_product_left(d(1:3), Iy(1:3), t%monoj(1:3)) - two * f_dip(2) * t%monoj(1:3) &
            + fp_quad(4) * DxDxJ_x &  ! Dipole - x a
            + fp_dip(1) * double_cross_product_left(Iy(1:3), d(1:3), t%dipjx(1:3)) &  !
            + fp_dip(1) * double_cross_product_left(d(1:3), Iy(1:3), t%dipjx(1:3)) &  !
            + fp_dip(2) * I_xxDxJ_x &  !b
            + fpmono * double_cross_product_left(Ix(1:3), Iy(1:3), t%dipjx(1:3)) &  !b
            + fp_dip(2) * DxI_xxJ_x &  !c
            + fpmono * double_cross_product_left(Iy(1:3), Ix(1:3), t%dipjx(1:3)) &  !c
            - two * f_quad(4) * t%dipjx(1:3) &  !d
            + fp_quad(2) * DxDxJ_y &  ! Dipole - y a
            + fp_dip(2) * I_yxDxJ_y * two &  !
            + fp_dip(2) * DxI_yxJ_y * two &  !
            !                     +fp_dip(2)*double_cross_product_left( Iy(1:3) , d(1:3), t%dipjy(1:3)  )                                 &  !b
            + fpmono * double_cross_product_left(Iy(1:3), Iy(1:3), t%dipjy(1:3)) &  !b
            !                     +fp_dip(2)*double_cross_product_left( d(1:3) , Iy(1:3), t%dipjy(1:3)  )                                 &  !c
            + fpmono * double_cross_product_left(Iy(1:3), Iy(1:3), t%dipjy(1:3)) &  !c
            - two * f_quad(2) * t%dipjy(1:3) &  !d
            + fp_quad(5) * DxDxJ_z &  ! Dipole - z a
            + fp_dip(3) * double_cross_product_left(Iy(1:3), d(1:3), t%dipjz(1:3)) &  ! Dipole - z a
            + fp_dip(3) * double_cross_product_left(d(1:3), Iy(1:3), t%dipjz(1:3)) &  ! Dipole - z a
            + fp_dip(2) * I_zxDxJ_z &  !b
            + fpmono * double_cross_product_left(Iz(1:3), Iy(1:3), t%dipjz(1:3)) &  !b
            + fp_dip(2) * DxI_zxJ_z &  !c
            + fpmono * double_cross_product_left(Iy(1:3), Iz(1:3), t%dipjz(1:3)) &  !c
            - two * f_quad(5) * t%dipjz(1:3) &  !d
            + fp_oct(4) * DxDxJ_xx & ! Quadrupole - xx a
            + fp_quad(1) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjx(1:3)) & !
            + fp_quad(1) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjx(1:3)) & !
            + fp_quad(4) * I_xxDxJ_xx * two & !b
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjx(1:3)) * two & !b
            + fp_quad(4) * DxI_xxJ_xx * two & !c
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjx(1:3)) * two & !c
            + fp_dip(2) * I_xxI_xxJ_xx * two & !d
            - two * f_oct(4) * t%quadjx(1:3) & !e
            + fp_oct(2) * DxDxJ_yy & ! Quadrupole - yy a
            !                     +fp_quad(2)*double_cross_product_left( Iy(1:3) , d(1:3), t%quadjy(1:3)  )                               & !
            !                     +fp_quad(2)*double_cross_product_left( d(1:3) , Iy(1:3), t%quadjy(1:3)  )                               & !
            + fp_quad(2) * I_yxDxJ_yy * three & !b
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two & !b
            + fp_quad(2) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjy(1:3)) * three & !c
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two & !c
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two & !d
            - two * f_oct(2) * t%quadjy(1:3) & !e
            + fp_oct(9) * double_cross_product_left(d(1:3), d(1:3), t%quadjz(1:3)) & ! Quadrupole - zz a
            + fp_quad(3) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjz(1:3)) & !
            + fp_quad(3) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjz(1:3)) & !
            + fp_quad(5) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjz(1:3)) * two & !b
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjz(1:3)) * two & !b
            + fp_quad(5) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjz(1:3)) * two & !c
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjz(1:3)) * two & !c
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two & !d
            - two * f_oct(9) * t%quadjz(1:3) & !e
            + fp_oct(6) * double_cross_product_left(d(1:3), d(1:3), t%quadjxy(1:3)) & ! Quadrupole - xy a
            + fp_quad(4) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjxy(1:3)) & !
            + fp_quad(4) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjxy(1:3)) & !
            + fp_quad(4) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjxy(1:3)) & ! b
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjxy(1:3)) & ! b
            + fp_quad(4) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjxy(1:3)) & ! d
            + fp_quad(2) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjxy(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! d
            + fp_quad(2) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! f
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! g
            - two * f_oct(6) * t%quadjxy(1:3) & ! h
            + fp_oct(7) * double_cross_product_left(d(1:3), d(1:3), t%quadjyz(1:3)) & ! Quadrupole - yz a
            + fp_quad(5) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjyz(1:3)) & !
            + fp_quad(5) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjyz(1:3)) & !
            + fp_quad(5) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjyz(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjyz(1:3)) & ! b
            + fp_quad(5) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_quad(2) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjyz(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! d
            + fp_quad(2) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! f
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! g
            - two * f_oct(7) * t%quadjyz(1:3) & ! h
            + fp_oct(10) * double_cross_product_left(d(1:3), d(1:3), t%quadjzx(1:3)) & ! Quadrupole - zx a
            + fp_quad(6) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjzx(1:3)) & !
            + fp_quad(6) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjzx(1:3)) & !
            + fp_quad(5) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjzx(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjzx(1:3)) & ! b
            + fp_quad(5) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_quad(4) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjzx(1:3)) & ! d
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjzx(1:3)) & ! d
            + fp_quad(4) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! f
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! g
            - two * f_oct(10) * t%quadjzx(1:3)                                                                            ! h

      dAz = fp_dip(3) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Iz(1:3), d(1:3), t%monoj(1:3)) &  !
            + fpmono * double_cross_product_left(d(1:3), Iz(1:3), t%monoj(1:3)) - two * f_dip(3) * t%monoj(1:3) &  !
            + fp_quad(6) * DxDxJ_x &  ! Dipole - x a
            + fp_dip(1) * double_cross_product_left(Iz(1:3), d(1:3), t%dipjx(1:3)) &  !
            + fp_dip(1) * double_cross_product_left(d(1:3), Iz(1:3), t%dipjx(1:3)) &  !
            + fp_dip(3) * I_xxDxJ_x &  !b
            + fpmono * double_cross_product_left(Ix(1:3), Iz(1:3), t%dipjx(1:3)) &  !b
            + fp_dip(3) * DxI_xxJ_x &  !c
            + fpmono * double_cross_product_left(Iz(1:3), Ix(1:3), t%dipjx(1:3)) &  !c
            - two * f_quad(6) * t%dipjx(1:3) &  !d
            + fp_quad(5) * DxDxJ_y &  ! Dipole - y a
            + fp_dip(2) * double_cross_product_left(Iz(1:3), d(1:3), t%dipjy(1:3)) &  !
            + fp_dip(2) * double_cross_product_left(d(1:3), Iz(1:3), t%dipjy(1:3)) &  !
            + fp_dip(3) * I_yxDxJ_y &  !b
            + fpmono * double_cross_product_left(Iy(1:3), Iz(1:3), t%dipjy(1:3)) &  !b
            + fp_dip(3) * DxI_yxJ_y &  !c
            + fpmono * double_cross_product_left(Iz(1:3), Iy(1:3), t%dipjy(1:3)) &  !c
            - two * f_quad(5) * t%dipjy(1:3) &  !d
            + fp_quad(3) * DxDxJ_z &  ! Dipole - z a
            + fp_dip(3) * I_zxDxJ_z * two &  !
            + fp_dip(3) * DxI_zxJ_z * two &  !
            !                     +fp_dip(3)*double_cross_product_left( Iz(1:3) , d(1:3), t%dipjz(1:3)  )                                 &  !b
            + fpmono * double_cross_product_left(Iz(1:3), Iz(1:3), t%dipjz(1:3)) &  !b
            !                     +fp_dip(3)*double_cross_product_left( d(1:3) , Iz(1:3), t%dipjz(1:3)  )                                 &  !c
            + fpmono * double_cross_product_left(Iz(1:3), Iz(1:3), t%dipjz(1:3)) &  !c
            - two * f_quad(3) * t%dipjz(1:3) &  !d
            + fp_oct(5) * DxDxJ_xx & ! Quadrupole - xx a
            + fp_quad(1) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjx(1:3)) & !
            + fp_quad(1) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjx(1:3)) & !
            + fp_quad(6) * I_xxDxJ_xx * two & ! b
            + fp_dip(1) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjx(1:3)) * two & ! b
            + fp_quad(6) * DxI_xxJ_xx * two & ! c
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjx(1:3)) * two & ! c
            + fp_dip(3) * I_xxI_xxJ_xx * two & ! d
            - two * f_oct(5) * t%quadjx(1:3) & ! e
            + fp_oct(7) * DxDxJ_yy & ! Quadrupole - yy a
            + fp_quad(2) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjy(1:3)) & !
            + fp_quad(2) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjy(1:3)) & !
            + fp_quad(5) * I_yxDxJ_yy * two & ! b
            + fp_dip(2) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjy(1:3)) * two & ! b
            + fp_quad(5) * I_yxI_xxJ_yy * two & ! c
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjy(1:3)) * two & ! c
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iy(1:3), t%quadjy(1:3)) * two & ! d
            - two * f_oct(7) * t%quadjy(1:3) & ! e
            + fp_oct(3) * double_cross_product_left(d(1:3), d(1:3), t%quadjz(1:3)) & ! Quadrupole - zz a
            + fp_quad(3) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjz(1:3)) & !
            + fp_quad(3) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjz(1:3)) & !
            + fp_quad(3) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjz(1:3)) * two & ! b
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! b
            + fp_quad(3) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! c
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! c
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjz(1:3)) * two & ! d
            - two * f_oct(3) * t%quadjz(1:3) & ! e
            + fp_oct(10) * double_cross_product_left(d(1:3), d(1:3), t%quadjxy(1:3)) & ! Quadrupole - xy a
            + fp_quad(4) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjxy(1:3)) & !
            + fp_quad(4) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjxy(1:3)) & !
            + fp_quad(6) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjxy(1:3)) & ! b
            + fp_dip(1) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjxy(1:3)) & ! b
            + fp_quad(6) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjxy(1:3)) & ! c
            + fp_quad(5) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjxy(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjxy(1:3)) & ! d
            + fp_quad(5) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjxy(1:3)) & ! e
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iy(1:3), t%quadjxy(1:3)) & ! f
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Ix(1:3), t%quadjxy(1:3)) & ! g
            - two * f_oct(10) * t%quadjxy(1:3) & ! h
            + fp_oct(9) * double_cross_product_left(d(1:3), d(1:3), t%quadjyz(1:3)) & ! Quadrupole - yz a
            + fp_quad(5) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjyz(1:3)) & !
            + fp_quad(5) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjyz(1:3)) & !
            + fp_quad(3) * double_cross_product_left(Iy(1:3), d(1:3), t%quadjyz(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! b
            + fp_quad(3) * double_cross_product_left(d(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! c
            + fp_quad(5) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjyz(1:3)) & ! d
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjyz(1:3)) & ! d
            + fp_quad(5) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(2) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjyz(1:3)) & ! e
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Iy(1:3), t%quadjyz(1:3)) & ! f
            + fp_dip(3) * double_cross_product_left(Iy(1:3), Iz(1:3), t%quadjyz(1:3)) & ! g
            - two * f_oct(9) * t%quadjyz(1:3) & ! h
            + fp_oct(8) * double_cross_product_left(d(1:3), d(1:3), t%quadjzx(1:3)) & ! Quadrupole - zx a
            + fp_quad(6) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjzx(1:3)) & !
            + fp_quad(6) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjzx(1:3)) & !
            + fp_quad(3) * double_cross_product_left(Ix(1:3), d(1:3), t%quadjzx(1:3)) & ! b
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! b
            + fp_quad(3) * double_cross_product_left(d(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! c
            + fp_quad(6) * double_cross_product_left(Iz(1:3), d(1:3), t%quadjzx(1:3)) & ! d
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjzx(1:3)) & ! d
            + fp_quad(6) * double_cross_product_left(d(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(1) * double_cross_product_left(Iz(1:3), Iz(1:3), t%quadjzx(1:3)) & ! e
            + fp_dip(3) * double_cross_product_left(Iz(1:3), Ix(1:3), t%quadjzx(1:3)) & ! f
            + fp_dip(3) * double_cross_product_left(Ix(1:3), Iz(1:3), t%quadjzx(1:3)) & ! g
            - two * f_oct(8) * t%quadjzx(1:3)                                                                             ! h

      phi = half * phi
      exyz = half * exyz
      Jxyz = half * Jxyz
      Bxyz = half * Bxyz
      Jirrxyz = half * Jirrxyz
      Axyz = half * Axyz
      dAx = half * dAx
      dAy = half * dAy
      dAz = half * dAz

   end subroutine calc_force_darwin_3D

!    subroutine calc_force_darwin_2D3V_direct(t, d, d2,eps2, phi, E, A, J, Jirr, B, Ex, Ey, Ax, Ay, Axx, Axy , Ayy )
   subroutine calc_force_darwin_2D3V_direct(t, d, d2, eps2, phi, E, A, J, Jirr, B, Ax, Ay)
      use module_tool, only: cross_product, double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
!      real(kind_physics), intent(out) ::  E(1:2),A(1:3),J(1:3),Jirr(1:3),B(1:3),phi,Ax(1:3),Ay(1:3),Axx(1:3),Ayy(1:3),Axy(1:3),&
!                                          Ex(1:2),Ey(1:2)
      real(kind_physics), intent(out) ::  E(1:2), A(1:3), J(1:3), Jirr(1:3), B(1:3), phi, Ax(1:3), Ay(1:3)

      real(kind_physics) :: dx, dy, rd2, rd4, rd6, rd8, dx2, dy2, dx3, dy3, r2, logR2e, over2, over4, over6, over8, &
                            phi_mono, phi_dip(1:2), rho_mono, rho_dip(1:2), phi_quad(1:3), &
                            E_mono(1:3), Ex_dip(1:3), Ey_dip(1:3), logTmp, &
                            A1_mono, A1_dip(1:2), A1_quad(1:3), A2_mono, A2_dip(1:2), &
                            dist(1:3), Jir_mono, over, Exy_quad(1:3), A2_quad(1:3), &
                            A3_mono(1:3), A3x_dip(1:3), A3y_dip(1:3), A3xx_quad(1:3), A3yy_quad(1:3), &
                            A3xy_quad(1:3), Btmp(1:3)!,A4_mono(1:3),A4x_dip(1:3),A4y_dip(1:3)

      dx = d(1)
      dy = d(2)
      dist(1:2) = d(1:2)
      dist(3) = zero
      r2 = dx**2 + dy**2
      over = sqrt(one / r2)
      over2 = one / r2
      over4 = over2 * over2
!      over5 = over**5
      over6 = over4 * over2
!      over7 = over**7
      over8 = over6 * over2
!      over9 = over**9
!      over10= over8*over2
      rd2 = one / d2
      rd4 = rd2 * rd2
      rd6 = rd4 * rd2
      rd8 = rd6 * rd2

      dx2 = dx * dx
      dy2 = dy * dy
      dx3 = dx2 * dx
      dy3 = dy2 * dy
!      dx4 = dx3*dx
!      dy4 = dy3*dy
!      dx5 = dx4*dx
!      dy5 = dy4*dy

      logTmp = log(r2 / eps2 + one)
      logR2e = eps2 * over2 * logTmp

      !!! Rho/J must be multiply by eps2/pi
      rho_mono = rd4
!      rho_dip(1:2)= -four*d(1:2)*rd6
!      rho_quad(1)=  twentyfour*dx2*rd8 - four*rd6
!      rho_quad(2)=  twentyfour*dy2*rd8 - four*rd6
!      rho_quad(3)=  twentyfour*dy*dx*rd8

      !!! Jirr must be multiply by 2
      Jir_mono = rd2
!      Jir_dip(1:2)  = two*d(1:2)*rd4
!      Jir_quad(1)   = two*rd4*( four*dx2*rd2 - one )
!      Jir_quad(2)   = two*rd4*( four*dy2*rd2 - one )
!      Jir_quad(3)   = eight*dx*dy*rd6

      phi_mono = -log(d2)
      phi_dip = -two * d(1:2) * rd2
      phi_quad(1) = four * dx2 * rd4 - two * rd2
      phi_quad(2) = four * dy2 * rd4 - two * rd2
      phi_quad(3) = four * dy * dx * rd4

      E_mono(1:2) = -phi_dip
      E_mono(3) = zero
      Ex_dip(1:3) = -(/phi_quad(1), phi_quad(3), zero/)
      Ey_dip(1:3) = -(/phi_quad(3), phi_quad(2), zero/)

!      Exx_quad(1)=  sixteen*dx3*rd6    - twelve*dx*rd4
!      Exx_quad(2)=  sixteen*dx2*dy*rd6 - four*dy*rd4  !! Tshis is E_xy(1)
!      Exx_quad(3)=  zero
!      Eyy_quad(1)=  sixteen*dx*dy2*rd6 - four*dx*rd4  !! This is E_xy(2)
!      Eyy_quad(2)=  sixteen*dy3*rd6    - twelve*dy*rd4
!      Eyy_quad(3)=  zero
!
!      Exy_quad(1)=  Exx_quad(2)
!      Exy_quad(2)=  Eyy_quad(1)
!      Exy_quad(3)=  zero
!
!      Exxx_oct(1)= ninetysix*dx2*rd6   - twelve*rd4 - ninetysix*dx2**2*rd8
!      Exxx_oct(2)= fortyeight*dx*dy*rd6             - ninetysix*dx3*dy*rd8
!      Exxx_oct(3)= zero
!
!      Exxy_oct(1)= fortyeight*dx*dy*rd6             - ninetysix*dx3*dy*rd8
!      Exxy_oct(2)= sixteen*dx2*rd6 - four*rd4 + sixteen*dy2*rd6 - ninety*dx2*dy2*rd8
!      Exxy_oct(3)= zero
!
!      Exyy_oct(1)= sixteen*dx2*rd6 - four*rd4 + sixteen*dy2*rd6 - ninety*dx2*dy2*rd8
!      Exyy_oct(2)= fortyeight*dx*dy*rd6             - ninetysix*dx*dy3*rd8
!      Exyy_oct(3)= zero
!
!      Eyyy_oct(1)= fortyeight*dx*dy*rd6             - ninetysix*dx*dy3*rd8
!      Eyyy_oct(2)= ninetysix*dy2*rd6   - twelve*rd4 - ninetysix*dy2**2*rd8
!      Eyyy_oct(3)= zero

      A1_mono = logR2e - one
      A1_dip(1:2) = two * eps2 * d(1:2) * rd2 * over2 - two * over2 * logR2e * d(1:2)
      A1_quad(1) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dx2 * rd2 * over2 &
                                         - two * dx2 * rd4 + four * logTmp * dx2 * over4)
      A1_quad(2) = two * eps2 * over2 * (rd2 - logTmp * over2 - four * dy2 * rd2 * over2 &
                                         - two * dy2 * rd4 + four * logTmp * dy2 * over4)
      A1_quad(3) = four * over2 * dx * dy * (two * logR2e * over2 - eps2 * rd4 - two * eps2 * rd2 * over2)
!
!      A1_oct(1)  =  twelve*eps2*dx*over2*( four*dx2*rd2*over4 - two*rd2*over2 - rd4 +&
!                    two*logTmp*over4 + two*dx2*rd4*over2 - four*dx2*logTmp*over6  ) + sixteen*eps2*dx3*rd6*over2 ! xxx
!
!      A1_oct(2)  =  four*eps2*dy*over2*( two*logTmp*over4 - rd4 - two*rd2*over2 + twelve*dx2*rd2*over4 +&
!                    six*dx2*rd4*over2 + four*dx2*rd6 - twelve*dx2*logTmp*over6  )                                ! xxy
!
!      A1_oct(3)  =  four*eps2*dx*over2*( two*logTmp*over4 - rd4 - two*rd2*over2 + twelve*dy2*rd2*over4 +&
!                    six*dy2*rd4*over2 + four*dy2*rd6 - twelve*dy2*logTmp*over6  )                                ! xxy
!
!      A1_oct(4)  =  twelve*eps2*dy*over2*( four*dy2*rd2*over4 - two*rd2*over2 - rd4 +&
!                    two*logTmp*over4 + two*dy2*rd4*over2 - four*dy2*logTmp*over6  ) + sixteen*eps2*dy3*rd6*over2 ! yyy
!
!      A1_sed(1)  =  twentyfour*twelve*dx2*over2*rd2*eps2 - twentyfour*eps2*over4*rd2 - twentyfour*sixteen*eps2*dx4*over8*rd2 +& !xxxx
!                    -twelve*eps2*over2*rd4 + twentyfour*eps2*over6*logTmp - twentyfour*twelve*dx2*eps2*over8*logTmp          +&
!                    +twentyfour*sixteen*eps2*dx4*over5*logTmp + twelve**2*dx2*over4*eps2*rd4 + ninetysix*dx2*over2*rd6       +&
!                    -twentyfour*eight*dx4*eps2*rd4*over6 - two*eight**2*dx4*eps2*rd6*over4 - ninetysix*dx4*eps2*rd8*over2
!
!      A1_sed(2)  =  twelve**2*dx*dy*eps2*rd2*over6 - twentyfour*sixteen*eps2*dx3*dy*rd2*over8 - twelve**2*dx*dy*eps2*over8*logTmp             +& !xxxy
!                    - twentyfour*eight*dx3*dy*eps2*rd6*over4 - ninetysix*dx3*dy*eps2*rd8*over2 + twentyfour*sixteen*eps2*dx3*dy*logTmp*over10 + &
!                    + nine*eight*dx*dy*eps2*rd4*over4 + fortyeight*dx*dy*eps2*rd6*over2
!
!      A1_sed(3)  =  fortyeight*eps2*rd2*over4 - eight*eps2*rd2*over4  - four*eps2*rd4*over2 +eight*eps2*logTmp*over6                          +& !xxyy
!                     - fortyeight*eps2*logTmp*over6 - twentyfour*sixteen*dx2*dy2*eps2*rd2*over8 + twentyfour*eps2*over2*rd4                   +&
!                     + sixteen*eps2*over2*rd6 - twentyfour*eight*dx2*dy2*eps2*rd4*over6 - two*eight**2*dx2*dy2*eps2*rd6*over4                 +&
!                     - ninetysix*dx2*dy2*eps2*rd8*over2 + twentyfour*sixteen*dx2*dy2*eps2*logTmp*over10
!
!      A1_sed(4)  =  twelve**2*dx*dy*eps2*rd2*over6 - twentyfour*sixteen*eps2*dx*dy3*rd2*over8 - twelve**2*dx*dy*eps2*over8*logTmp             +& !xyyy
!                    - twentyfour*eight*dx*dy3*eps2*rd6*over4 - ninetysix*dx*dy3*eps2*rd8*over2 + twentyfour*sixteen*eps2*dx*dy3*logTmp*over10 + &
!                    + nine*eight*dx*dy*eps2*rd4*over4 + fortyeight*dx*dy*eps2*rd6*over2
!
!      A1_sed(5)  =  twentyfour*twelve*dy2*over2*rd2*eps2 - twentyfour*eps2*over4*rd2 - twentyfour*sixteen*eps2*dy4*over8*rd2 +& !yyyy
!                    -twelve*eps2*over2*rd4 + twentyfour*eps2*over6*logTmp - twentyfour*twelve*dy2*eps2*over8*logTmp          +&
!                    +twentyfour*sixteen*eps2*dy4*over5*logTmp + twelve**2*dy2*over4*eps2*rd4 + ninetysix*dy2*over2*rd6       +&
!                    -twentyfour*eight*dy4*eps2*rd4*over6 - two*eight**2*dy4*eps2*rd6*over4 - ninetysix*dy4*eps2*rd8*over2
!
      A2_mono = -A1_mono + one + phi_mono + log(eps2)
      A2_dip(1:2) = -A1_dip(1:2) + phi_dip(1:2)
      A2_quad(1:3) = -A1_quad(1:3) + phi_quad(1:3)
!      A2_oct(1)   = -A1_oct(1)     - Exx_quad(1)
!      A2_oct(2)   = -A1_oct(2)     - Exy_quad(1)
!      A2_oct(3)   = -A1_oct(2)     - Eyy_quad(1)
!      A2_oct(4)   = -A1_oct(3)     - Eyy_quad(2)
!      A2_sed(1)   = -A1_sed(1)     - Exxx_oct(1)
!      A2_sed(2)   = -A1_sed(2)     - Exxy_oct(1)
!      A2_sed(3)   = -A1_sed(3)     - Exyy_oct(1)
!      A2_sed(4)   = -A1_sed(4)     - Eyyy_oct(1)
!      A2_sed(5)   = -A1_sed(5)     - Eyyy_oct(2)

      A3_mono(1:3) = over * dist(1:3)
      A3x_dip(1) = over - dx2 * over * over2
      A3x_dip(2) = -dx * dy * over * over2
      A3x_dip(3) = zero
      A3y_dip(1) = -dx * dy * over * over2
      A3y_dip(2) = over - dy2 * over * over2
      A3y_dip(3) = zero
      A3xx_quad(1) = -three * over * over2 * dx * (one - dx2 * over2)
      A3xx_quad(2) = -over * over2 * dy * (one - three * dx2 * over2)   ! This is equal to A3_xy(1)
      A3xx_quad(3) = zero
      A3yy_quad(1) = -over * over2 * dx * (one - three * dy2 * over2)     ! This is equal to A3_xy(2)
      A3yy_quad(2) = -three * over * over2 * dy * (one - dy2 * over2)
      A3yy_quad(3) = zero
      A3xy_quad(1) = A3xx_quad(2)
      A3xy_quad(2) = A3yy_quad(1)
      A3xy_quad(3) = zero
!      A3xxx_oct(1)  = three*over**3*( six*dx2*over2 - one - five*dx2**2*over4  )
!      A3xxx_oct(2)  = three*over**5*dx*dy*( three - five*dx2*over2  )
!      A3xxx_oct(3)  = zero
!
!      A3xxy_oct(1)  = three*over**5*dx*dy*( three - five*dx2*over2  )
!      A3xxy_oct(2)  = over**3*( two - fifteen*dx2*dy2*over4  )
!      A3xxy_oct(3)  = zero
!
!      A3xyy_oct(1)  = over**3*( two - fifteen*dx2*dy2*over4  )
!      A3xyy_oct(2)  = three*over**5*dx*dy*( three - five*dy2*over2  )
!      A3xyy_oct(3)  = zero
!
!      A3yyy_oct(1)  = three*over**5*dx*dy*( three - five*dy2*over2  )
!      A3yyy_oct(2)  = three*over**3*( six*dy2*over2 - one - five*dy2**2*over4  )
!      A3yyy_oct(3)  = zero
!
!      A3xxxx_sed(1) = fortyfive*dx*over5 - hundredfifty*dx3*over7 + hundredfive*dx5*over9
!      A3xxxx_sed(2) = nine*dy*over5      - ninety*dx2*dy*over7    + hundredfive*dx4*dy*over9
!      A3xxxx_sed(3) = zero
!
!      A3xxxy_sed(1) = nine*dy*over5 - ninety*dx2*dy*over7   + hundredfive*dx4*dy*over9
!      A3xxxy_sed(2) = nine*dx*over5 - fifteen*dx3*over7     + hundredfive*dx3*dy2*over9 - fortyfive*dx*dy2*over7
!      A3xxxy_sed(3) = zero
!
!      A3xxyy_sed(1) = nine*dx*over5 - fifteen*dx3*over7     + hundredfive*dx3*dy2*over9 - fortyfive*dx*dy2*over7
!      A3xxyy_sed(2) = nine*dy*over5 - fifteen*dy3*over7     + hundredfive*dx2*dy3*over9 - fortyfive*dx2*dy*over7
!      A3xxyy_sed(3) = zero
!
!      A3xyyy_sed(1) = nine*dy*over5 - fifteen*dy3*over7     + hundredfive*dx2*dy3*over9 - fortyfive*dx2*dy*over7
!      A3xyyy_sed(2) = nine*dx*over5 - ninety*dx*dy2*over7   + hundredfive*dx*dy4*over9
!      A3xyyy_sed(3) = zero
!
!      A3yyyy_sed(1) = nine*dx*over5      - ninety*dx*dy2*over7    + hundredfive*dx*dy4*over9
!      A3yyyy_sed(2) = fortyfive*dy*over5 - hundredfifty*dy3*over7 + hundredfive*dy5*over9
!      A3yyyy_sed(3) = zero

!      A4_mono(1:3)      = dist(1:3)
!      A4x_dip(1)        = one
!      A4x_dip(2)        = zero
!      A4x_dip(3)        = zero
!      A4y_dip(1)        = zero
!      A4y_dip(2)        = one
!      A4y_dip(3)        = zero

      phi = (t%charge * phi_mono)                                                                   !& Monopole
      phi = half * phi

      E(1:2) = (t%charge * E_mono(1:2))                                                             !& Monopole

      E(1:2) = half * E(1:2)

!      Ex(1:2)    = ( t%charge*Ex_dip(1:2)   )                                                      !& Monopole
!
!      Ex(1:2)    = half*Ex(1:2)
!
!
!      Ey(1:2)    = ( t%charge*Ey_dip(1:2)   )                                                      !& Monopole
!
!      Ey(1:2)    = half*Ey(1:2)

      J(1:3) = (t%monoj(1:3) * rho_mono)                                                            !& Monopole

      J(1:3) = half * eps2 / pi * J(1:3)

      Jirr(1:3) = (t%monoj(1:3) * (three * half * Jir_mono - eps2 * rho_mono) &
                   + double_cross_product_left(E_mono(1:3), E_mono(1:3), t%monoj(1:3)))             !& Monopole
      Jirr(1:3) = half * oneoverpi * Jirr(1:3)

      Btmp(1:3) = cross_product(t%monoj(1:3), E_mono(1:3))                                          !& Monopole
      B(1:3) = half * Btmp(1:3)

      A(1:3) = half * t%monoj(1:3) * A2_mono &
               + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3))      !& Monopole
!
      A(3) = t%monoj(3) * phi_mono                                                                  !& Monopole

      A(1:3) = half * A(1:3)

      Ax(1:3) = half * t%monoj(1:3) * A2_dip(1) &                                                   !& Monopole
                + A1_mono * double_cross_product_left(A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)) &
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3x_dip(1:3)) &
                + A1_dip(1) * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3))

      Ax(3) = t%monoj(3) * phi_dip(1)                                                               !& Monopole

      Ay(1:3) = half * t%monoj(1:3) * A2_dip(2) &
                + A1_mono * double_cross_product_left(A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)) &   !& Monopole
                + A1_mono * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)) &   !& Monopole
                + A1_dip(2) * double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3))   !& Monopole

!
      Ay(3) = t%monoj(3) * phi_dip(2)                                                               !& Monopole

      Ax(1:3) = half * Ax(1:3)
      Ay(1:3) = half * Ay(1:3)

!
!
!
!      Axx(1:3)   =  half*t%monoj(1:3)*A2_quad(1)                                                            &! Monopole
!                 + A1_quad(1)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(1)*double_cross_product_left( A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(1)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3x_dip(1:3)  )         &
!
!                 + A1_dip(1)*double_cross_product_left(A3x_dip(1:3) , t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3xx_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3x_dip(1:3)  , t%monoj(1:3), A3x_dip(1:3)  )         &
!
!                 + A1_dip(1)*double_cross_product_left(A3_mono(1:3) , t%monoj(1:3), A3x_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3x_dip(1:3)  , t%monoj(1:3), A3x_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3xx_quad(1:3))
!
!      Axx(3)     =   t%monoj(3)*phi_quad(1)                                                                   ! Monopole
!
!
!
!      Ayy(1:3)   =  half*t%monoj(1:3)*A2_quad(2)                                                            &! Monopole
!                 + A1_quad(2)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(2)*double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(2)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)  )         &
!
!                 + A1_dip(2)*double_cross_product_left(A3y_dip(1:3) , t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3yy_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3y_dip(1:3)  , t%monoj(1:3), A3y_dip(1:3)  )         &
!
!                 + A1_dip(2)*double_cross_product_left(A3_mono(1:3) , t%monoj(1:3), A3y_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3y_dip(1:3)  , t%monoj(1:3), A3y_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3yy_quad(1:3))
!
!!
!      Ayy(3)     =   t%monoj(3)*phi_quad(2)                                                                   ! Monopole
!
!
!      Axy(1:3)   =  half*t%monoj(1:3)*A2_quad(3)                                                            &! Monopole
!                 + A1_quad(3)*double_cross_product_left(A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(1)*double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_dip(1)*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3y_dip(1:3)  )         &
!
!                 + A1_dip(2)*double_cross_product_left(A3x_dip(1:3) , t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3xy_quad(1:3), t%monoj(1:3), A3_mono(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3x_dip(1:3)  , t%monoj(1:3), A3y_dip(1:3)  )         &
!
!                 + A1_dip(2)*double_cross_product_left(A3_mono(1:3) , t%monoj(1:3), A3x_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3y_dip(1:3)  , t%monoj(1:3), A3x_dip(1:3)  )         &
!                 + A1_mono*double_cross_product_left( A3_mono(1:3)  , t%monoj(1:3), A3xy_quad(1:3))
!
!      Axy(3)     =   t%monoj(3)*phi_quad(3)                                                                   ! Monopole
!
!
!      Axx(1:3)   = half*Axx(1:3)
!      Ayy(1:3)   = half*Ayy(1:3)
!      Axy(1:3)   = half*Axy(1:3)
!

      A(1:2) = -A(1:2)
      Ax(1:2) = -Ax(1:2)
      Ay(1:2) = -Ay(1:2)
!      Axx(1:2)   = -Axx(1:2)
!      Axy(1:2)   = -Axy(1:2)
!      Ayy(1:2)   = -Ayy(1:2)

   end subroutine calc_force_darwin_2D3V_direct

   subroutine calc_force_darwin_3D_direct(t, d, d2, eps2, phi, Exyz, Axyz, Jxyz, Jirrxyz, Bxyz, dAx, dAy, dAz)
      use module_tool, only: cross_product, double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(1:3), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
!      real(kind_physics), intent(out) ::  E(1:2),A(1:3),J(1:3),Jirr(1:3),B(1:3),phi,Ax(1:3),Ay(1:3),Axx(1:3),Ayy(1:3),Axy(1:3),&
!                                          Ex(1:2),Ey(1:2)
      real(kind_physics), intent(out) ::  Exyz(1:3), Axyz(1:3), Jxyz(1:3), Jirrxyz(1:3), Bxyz(1:3), phi, dAx(1:3), dAy(1:3), dAz(1:3)

      real(kind_physics) :: dx, dy, dz, r2, r, rd, rEps, rd2, over, over2, over3, over4, &!over5,over7,over6,over8,over9, &
                            dx2, dy2, dz2, rd3, rd5, &!,dx3,dy3,dz3rd6,rd10,rd12,rd14 ,rd4,rd7,rd8,rd9                 , &
                            rho_mono, &!,rho_dip,rho_quad(1:4),&                                     !density/current density
                            phi_mono, &!phi_dip,phi_quad(1:4),&                                     !potential
                            E_mono(1:3), &!,E_dipx(1:3),E_dipy(1:3),E_dipz(1:3),E_quadxx(1:3),E_quadyy(1:3),E_quadzz(1:3),& !Eirr
                            !                            E_quadxy(1:3),E_quadyz(1:3),E_quadzx(1:3),&
                            fmono, fpmono, A1_mono, A1_dip, A1_quad, A1_oct, &                         !A
                            A2_mono, A2_dip, f_dip(1:3), fp_dip(1:3), Ix(1:3), Iy(1:3), Iz(1:3), fppmono, DxDxJ(1:3)

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r2 = dx**2 + dy**2 + dz**2
      r = sqrt(r2)
      rEps = sqrt(d2)
      rd2 = one / d2 ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd = sqrt(rd2)
      rd3 = rd**3
!      rd4   = rd2**2
      rd5 = rd**5
!      rd6   = rd2**3
!      rd7   = rd**7
!      rd8   = rd2**4
!      rd9   = rd**9
!      rd10  = rd2**5
!      rd12  = rd2**6
!      rd14  = rd2**7
      over = one / r
      over2 = over * over
      over3 = over2 * over
      over4 = over3 * over
!      over5 = over4*over
!      over6 = over5*over
!      over7 = over6*over
!      over8 = over7*over
!      over9 = over8*over

      dx2 = dx * dx
      dy2 = dy * dy
      dz2 = dz * dz
!      dx3 = dx*dx2
!      dy3 = dy*dy2
!      dz3 = dz*dz2

      Ix(1:3) = (/one, zero, zero/)
      Iy(1:3) = (/zero, one, zero/)
      Iz(1:3) = (/zero, zero, one/)

      rho_mono = rd5

      phi_mono = rd
!      phi_dip        = -rd3
!      phi_quad(1)    = three*dx2*rd5 - rd3  !xx
!      phi_quad(2)    = three*dy2*rd5 - rd3  !yy
!      phi_quad(3)    = three*dz2*rd5 - rd3  !zz
!
!      phi_quad(4)    = three*rd5         !mixed terms

      E_mono = rd3 * d(1:3)
!      E_dipx       = (/ rd3 - three*rd5*dx2 , - three*dx*dy*rd5  , - three*dx*dz*rd5 /)
!      E_dipy       = (/ - three*dx*dy*rd5   ,rd3 - three*rd5*dy2 , - three*dy*dz*rd5 /)
!      E_dipz       = (/ - three*dx*dz*rd5   , - three*dy*dz*rd5  ,rd3 - three*rd5*dz2/)
!      E_quadxx     = (/ fifteen*rd7*dx3 - nine*dx*rd5 , fifteen*rd7*dx2*dy - three*dy*rd5  , fifteen*rd7*dx2*dz - three*dz*rd5   /)
!      E_quadyy     = (/ fifteen*rd7*dx*dy2 - three*dx*rd5 , fifteen*rd7*dy3 - nine*dy*rd5  , fifteen*rd7*dy2*dz - three*dz*rd5   /)
!      E_quadzz     = (/ fifteen*rd7*dx*dz2 - three*dx*rd5 , fifteen*rd7*dz2*dy - three*dy*rd5  , fifteen*rd7*dz3 - nine*dz*rd5   /)
!      E_quadxy     = (/ fifteen*rd7*dx2*dy - three*dy*rd5 , fifteen*rd7*dx*dy2 - three*dx*rd5  , fifteen*dx*dy*dz*rd7   /)
!      E_quadyz     = (/ E_quadxy(3) , fifteen*rd7*dy2*dz - three*dz*rd5 , fifteen*rd7*dz2*dy - three*dy*rd5             /)
!      E_quadzx     = (/  fifteen*rd7*dx2*dz - three*dz*rd5 , E_quadxy(3) , fifteen*rd7*dz2*dx - three*dx*rd5            /)

      Jxyz = three / four / pi * eps2 * rho_mono * t%monoj(1:3)

      phi = (t%charge * phi_mono)

      Exyz = t%charge * E_mono

      Bxyz = cross_product(t%monoj(1:3), E_mono(1:3))

      A1_mono = over2
      A1_dip = -two * over4 ! first order derivative

      A2_mono = A1_mono * rd
      A2_dip = -over2 * rd * (two * d2 + r2) * over2 * rd2 ! first order derivative

      fmono = -half * rEps * over2 + half * eps2 * over3 * log(rEps + r)
      fpmono = -three * over2 * fmono - over2 * rd
      fppmono = -(three * A1_mono * fpmono + three * fmono * A1_dip + A2_dip)

      f_dip(1:3) = fpmono * (/dx, dy, dz/)
      fp_dip = fppmono * (/dx, dy, dz/)

      DxDxJ = double_cross_product_left(d(1:3), d(1:3), t%monoj(1:3))

      Axyz = fpmono * DxDxJ - two * fmono * t%monoj(1:3)             ! Monopole

      dAx = fp_dip(1) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Ix(1:3), d(1:3), t%monoj(1:3)) &
            + fpmono * double_cross_product_left(d(1:3), Ix(1:3), t%monoj(1:3)) &
            - two * f_dip(1) * t%monoj(1:3)

      dAy = fp_dip(2) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Iy(1:3), d(1:3), t%monoj(1:3)) &
            + fpmono * double_cross_product_left(d(1:3), Iy(1:3), t%monoj(1:3)) &
            - two * f_dip(2) * t%monoj(1:3)

      dAz = fp_dip(3) * DxDxJ &  ! Monopole
            + fpmono * double_cross_product_left(Iz(1:3), d(1:3), t%monoj(1:3)) &
            + fpmono * double_cross_product_left(d(1:3), Iz(1:3), t%monoj(1:3)) &
            - two * f_dip(3) * t%monoj(1:3)

      phi = half * phi
      exyz = half * exyz
      Jxyz = half * Jxyz
      Bxyz = half * Bxyz
      Jirrxyz = half * Jirrxyz
      Axyz = half * Axyz
      dAx = half * dAx
      dAy = half * dAy
      dAz = half * dAz

   end subroutine calc_force_darwin_3D_direct

   !>
   !> Calculates 3D Kelbg interaction of particle p with particle inode
   !> that is shifted by the lattice vector vbox
   !> results are returned in exyz, phi
   !>
   subroutine calc_force_kelbg_3D_direct(particle, t, d, dist2, kelbg_invsqrttemp, exyz, phi)
      implicit none

      type(t_particle), intent(inout) :: particle
      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi
      real(kind_physics), intent(in) :: kelbg_invsqrttemp
      real(kind_physics) :: rd, r, rd3
      real(kind_physics), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kind_physics) :: ome, rol, lambda, q, fprefac

      q = t%charge

      ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
      if (particle%data%q * q .lt. 0.) then
         ! e-i or i-e interaction
         lambda = 1.00027227_8 * kelbg_invsqrttemp
      else
         if (q .gt. 0.) then
            ! i-i interaction
            lambda = 0.03300355_8 * kelbg_invsqrttemp
         else
            ! e-e interaction
            lambda = 1.41421356_8 * kelbg_invsqrttemp
         end if
      end if

      r = sqrt(dist2)
      rd = one / r
      rd3 = rd * rd * rd
      rol = r / lambda        !< "r over lambda"
      ome = 1 - exp(-rol * rol) !< "one minus exp(stuff)"

      ! potential
      phi = q * rd * (ome + sqrtpi * rol * (1 - erf(rol)))
      !  forces
      fprefac = q * rd3 * ome
      exyz = fprefac * d
   end subroutine calc_force_kelbg_3D_direct
end module

