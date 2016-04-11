! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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

    public calc_force_coulomb_3D
    public calc_force_coulomb_3D_direct
    public calc_force_coulomb_2D
    public calc_force_coulomb_2D_direct
    public calc_force_LJ
    public calc_force_kelbg_3D_direct
    public calc_force_darwin_2D
    public calc_force_darwin_2D3V
    public calc_force_darwin_2D_direct
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

      real(kind_physics) :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6, m2rd5, m5rd7, &
        pre1, pre2x, pre2y, pre2z, preQ1, preQ2, preQ3

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r  = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd = one/r
      rd2 = rd *rd
      rd3 = rd *rd2
      rd5 = rd3*rd2
      rd7 = rd5*rd2

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dx3 = dx*dx2
      dy3 = dy*dy2
      dz3 = dz*dz2

      fd1 = three*dx2*rd5 - rd3
      fd2 = three*dy2*rd5 - rd3
      fd3 = three*dz2*rd5 - rd3
      fd4 = three*dx*dy*rd5
      fd5 = three*dy*dz*rd5
      fd6 = three*dx*dz*rd5

      m2rd5 = two*rd5
      m5rd7 = five*rd7
      pre1 = m5rd7*dx*dy*dz
      pre2x = m5rd7*dx2 - rd5
      pre2y = m5rd7*dy2 - rd5
      pre2z = m5rd7*dz2 - rd5
      preQ1 = pre2x*t%quad(1)
      preQ2 = pre2y*t%quad(2)
      preQ3 = pre2z*t%quad(3)

      phi = t%charge*rd                                             &  !  monopole term
            - (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
            + fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3)         &  !  quadrupole
            + fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

      exyz(1) = t%charge*dx*rd3                                     &  ! monopole term
                - (fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3))      &  ! dipole term
                + three * (                                         &  ! quadrupole term
                   dx * (                                           &
                       ( pre2x - m2rd5 )*t%quad(1)                  &
                     + preQ2                                        &
                     + preQ3                                        &
                   )                                                &
                   + dy*pre2x*t%xyquad                              &
                   + dz*pre2x*t%zxquad                              &
                   + pre1*t%yzquad                                  &
                  )

      exyz(2) = t%charge*dy*rd3                                     &
                - (fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3))      &
                + three * (                                         &
                   dy * (                                           &
                       ( pre2y - m2rd5 )*t%quad(2)                  &
                     + preQ1                                        &
                     + preQ3                                        &
                   )                                                &
                   + dx*pre2y*t%xyquad                              &
                   + dz*pre2y*t%yzquad                              &
                   + pre1*t%zxquad                                  &
                  )

      exyz(3) = t%charge*dz*rd3                                     &
                - (fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1))      &
                + three * (                                         &
                   dz * (                                           &
                     + ( pre2z - m2rd5 )*t%quad(3)                  &
                     + preQ2                                        &
                     + preQ1                                        &
                   )                                                &
                   + dx*pre2z*t%zxquad                              &
                   + dy*pre2z*t%yzquad                              &
                   + pre1*t%xyquad                                  &
                  )
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
      real(kind_physics), intent(out) ::  exy(1:2),phi

      real(kind_physics) :: dx,dy,rd2,rd4,rd6,dx2,dy2,dx3,dy3

      dx = d(1)
      dy = d(2)

      rd2 = one/d2
      rd4 = rd2*rd2
      rd6 = rd4*rd2

      dx2 = dx *dx
      dy2 = dy *dy
      dx3 = dx2*dx
      dy3 = dy2*dy

      phi = - half * t%charge * log(d2)            & !  monopole term
           - (dx * t%dip(1) + dy * t%dip(2)) * rd2 & !  dipole
           + t%quad(1) * (two * dx2 * rd4 - rd2)   & ! quadrupole
           + t%quad(2) * (two * dy2 * rd4 - rd2)   &
           + two * t%xyquad * dx * dy * rd4

      exy(1) = t%charge * dx * rd2                                 & ! monopole
           - t%dip(1) * (two * dx2 * rd4 - rd2)                    & ! dipole
           - t%dip(2) * two * dx * dy * rd4                        &
           + t%quad(1) * (eight * dx3 * rd6 - six * dx * rd4)      & ! quadrupole
           + t%quad(2) * (eight * dx * dy2 * rd6 - two * dx * rd4) &
           + t%xyquad * (eight * dx2 * dy * rd6 - two * dy * rd4)

      exy(2) = t%charge * dy * rd2                                 & ! monopole
           - t%dip(2) * (two * dy2 * rd4 - rd2)                    & ! dipole
           - t%dip(1) * two * dx * dy * rd4                        &
           + t%quad(2) * (eight * dy3 * rd6 - six * dy * rd4)      & ! quadrupole
           + t%quad(1) * (eight * dy * dx2 * rd6 - two * dy * rd4) &
           + t%xyquad * (eight * dy2 * dx * rd6 - two * dx * rd4)

    end subroutine calc_force_coulomb_2D


    subroutine calc_force_darwin_2D(t, d, d2,eps2, exy, Axy, Jxy, Jirrxy, Bz, phi)
      use module_globals, only: vtilde  
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(1:2),Axy(1:2),Jxy(1:2),Jirrxy(1:2),Bz,phi

      real(kind_physics) :: dx,dy,rd2,rd4,rd6,rd8,dx2,dy2,dx3,dy3,r2,logR2e,over2,over4,over6,over8
      real(kind_physics) :: Ax1,Ax2,Ay1,Ay2,Ay3,A0,A1,A2,A3,A4,A5,A6,A7,A8,vdotxM,vdotxDx,vdotxDy
      real(kind_physics) :: vdotxQxx,vdotxQxy,vdotxQyy, vcrossxM, vcrossxDx,vcrossxDy,vcrossxQx,vcrossxQxy,vcrossxQy
      real(kind_physics) :: rho0,rho1,rho2,rho3,phix1,phiy1,phix2,phiy2,phixy,E1,E2,E3,E4,Jirr_dipx(1:2),Jirr_dipy(1:2)
      real(kind_physics) :: Jirr_quadx(1:2),Jirr_quady(1:2),Jirr_quadxy(1:2)



      dx = d(1)
      dy = d(2)

      r2 =  dx**2 + dy**2
      over2 = one/r2
      over4 = over2*over2
      over6 = over4*over2
      over8 = over6*over2
      rd2 = one/d2
      rd4 = rd2*rd2
      rd6 = rd4*rd2
      rd8 = rd6*rd2

      dx2 = dx *dx
      dy2 = dy *dy
      dx3 = dx2*dx
      dy3 = dy2*dy

      phix1 = - two*rd2*dx
      phix2 = two*( - rd2 + two*rd4*dx2 )
      phiy1 = - two*rd2*dy
      phiy2 = two*( - rd2 + two*rd4*dy2 )
      phixy = four*rd4*dx*dy

      phi = - t%charge * log(d2)                                          &     !  monopole term
             + phix1*t%dip(1) + phiy1*t%dip(2)                            &     !  dipole
             + phix2*t%quad(1) + phiy2*t%quad(2) + phixy*t%xyquad               !  quadrupole



      E1     =  two*rd4*( - six + eight*rd2*dx2  )
      E2     =  two*rd4*( - six + eight*rd2*dy2  )
      E3     =  two*rd4*( - two + eight*rd2*dx2  )
      E4     =  two*rd4*( - two + eight*rd2*dy2  )


      exy(1) = - phix1*t%charge                                           &     ! monopole
                - phix2*t%dip(1) - phixy*t%dip(2)                         &     ! dipole
                + dx*E1*t%quad(1) + dx*E4*t%quad(2) + dy*E3*t%xyquad            ! quadrupole


      exy(2) = - phiy1*t%charge                                           &     ! monopole
                - phixy*t%dip(1) - phiy2*t%dip(2)                         &     ! dipole
                + dy*E3*t%quad(1) + dy*E2*t%quad(2) + dx*E4*t%xyquad            ! quadrupole


      rho0   = eps2/pi
      rho1   = four*rho0*rd6*( six*dx2*rd2 - one )
      rho2   = four*rho0*rd6*( six*dy2*rd2 - one )
      rho3   = twentyfour*rho0*rd8*dx*dy


!                !!!!!!!!!rho
      Jxy    = rho0*rd4*t%charge                                 &     ! monopole
             - four*rho0*rd6*( t%dip(1)*dx + t%dip(2)*dy )       &     ! dipole
             + rho1*t%quad(1) + rho2*t%quad(2) + rho3*t%xyquad         ! quadrupole

!      Jxy(1) = rho0*rd4*t%monoj(1)                                       &     ! monopole
!              - four*rho0*rd6*( t%dipjx(1)*dx + t%dipjy(1)*dy  )         &     ! dipole
!              + rho1*t%quadjx(1) + rho2*t%quadjy(1) + rho3*t%quadjxy(1)        ! quadrupole
!
!
!      Jxy(2) = rho0*rd4*t%monoj(2)                                       &     ! monopole
!              - four*rho0*rd6*( t%dipjx(2)*dx + t%dipjy(2)*dy  )         &     ! dipole
!              + rho1*t%quadjx(2) + rho2*t%quadjy(2) + rho3*t%quadjxy(2)        ! quadrupole


      vcrossxM   = t%monoj(1)*dy - t%monoj(2)*dx
      vcrossxDx  = t%dipjx(1)*dy - t%dipjx(2)*dx
      vcrossxDy  = t%dipjy(1)*dy - t%dipjy(2)*dx
      vcrossxQx  = t%quadjx(1)*dy - t%quadjx(2)*dx
      vcrossxQy  = t%quadjy(1)*dy - t%quadjy(2)*dx
      vcrossxQxy = t%quadjxy(1)*dy - t%quadjxy(2)*dx


      Bz     = two*rd2 *vcrossxM                                                                     &   ! - monopole
              - two *rd2*t%dipjx(2) - four *rd4*dx*vcrossxDx                                         &   ! - x dipole
              + two *rd2*t%dipjy(1) - four *rd4*dy*vcrossxDy                                         &   ! - y dipole
              - four *rd4*( vcrossxQx - two*t%quadjx(2)*dx ) + four*four*rd6 *dx2*vcrossxQx          &   ! - xx  quadrupole
              - four *rd4*( two*t%quadjy(1)*dy + vcrossxQy ) + four*four*rd6 *dy2*vcrossxQy          &   ! - yy  quadrupole
              - four *rd4*( t%quadjxy(1)*dx - t%quadjxy(2)*dy ) + four*four *rd6*dx*dy*vcrossxQxy        ! - xy quadrupole


      logR2e = log( one + r2/eps2 )

      vdotxM   = t%monoj(1)*dx + t%monoj(2)*dy
      vdotxDx  = t%dipjx(1)*dx + t%dipjx(2)*dy
      vdotxDy  = t%dipjy(1)*dx + t%dipjy(2)*dy
      vdotxQxx = t%quadjx(1)*dx + t%quadjx(2)*dy
      vdotxQyy = t%quadjy(1)*dx + t%quadjy(2)*dy
      vdotxQxy = t%quadjxy(1)*dx + t%quadjxy(2)*dy


      Jirr_dipx(1)    = eight*dx2*vdotxDx*rd6 - two*vdotxDx*rd4 - four*t%dipjx(1)*dx*rd4
      Jirr_dipx(2)    = eight*dx*dy*vdotxDx*rd6 - two*t%dipjx(2)*dx*rd4 - two*t%dipjx(1)*dy*rd4
      Jirr_dipy(1)    = eight*dx*dy*vdotxDy*rd6 - two*t%dipjy(2)*dx*rd4 - two*t%dipjy(1)*dy*rd4
      Jirr_dipy(2)    = eight*dy2*vdotxDy*rd6 - two*vdotxDy*rd4 - four*t%dipjy(2)*dy*rd4

      Jirr_quadx(1)   = twentyfour*dx*vdotxQxx*rd6 - six*t%quadjx(1)*rd4 - two*twentyfour*dx3*vdotxQxx*rd8 + twentyfour*t%quadjx(1)*dx2*rd6
      Jirr_quadx(2)   = eight*dy*vdotxQxx*rd6 - two*t%quadjx(2)*rd4 + eight*t%quadjx(2)*dx2*rd6 - two*twentyfour*dx2*dy*vdotxQxx*rd8 + four**2*t%quadjx(1)*dx*dy*rd6
      Jirr_quady(1)   = eight*dx*vdotxQyy*rd6 - two*t%quadjy(1)*rd4 + eight*t%quadjy(1)*dy2*rd6 - two*twentyfour*dx*dy2*vdotxQyy*rd8 + four**2*t%quadjy(2)*dx*dy*rd6
      Jirr_quady(2)   = twentyfour*dy*vdotxQyy*rd6 - six*t%quadjy(2)*rd4 - two*twentyfour*dy3*vdotxQyy*rd8 + twentyfour*t%quadjy(2)*dy2*rd6
      Jirr_quadxy(1)  = eight*dy*vdotxQxy*rd6 - two*t%quadjxy(2)*rd4 + eight*t%quadjxy(2)*dx2*rd6 - two*twentyfour*dx2*dy*vdotxQxy*rd8 + four**2*t%quadjxy(1)*dx*dy*rd6
      Jirr_quadxy(2)  = eight*dx*vdotxQxy*rd6 - two*t%quadjxy(1)*rd4 + eight*t%quadjxy(1)*dy2*rd6 - two*twentyfour*dx*dy2*vdotxQxy*rd8 + four**2*t%quadjxy(2)*dx*dy*rd6

      Jirrxy(1:2)   =  ( t%monoj(1:2)*rd2 - two*d(1:2)*rd4*vdotxM )                        &! Monopole
                         + (  Jirr_dipx(1:2) + Jirr_dipy(1:2)  )                           &! Dipole
                         + ( Jirr_quadx(1:2) + Jirr_quady(1:2) + Jirr_quadxy(1:2) )         ! Quadrupole

      Jirrxy(1:2)   = two*Jirrxy(1:2)

      A0       = eps2*over2*logR2e  - log(d2)
      A1       = eps2*over2*logR2e - one
      A2       = rd2 - over2*logR2e
      A3       = two*over2*rd2 + rd4 - two*over4*logR2e
      A4       = eps2*over2*rd2 - rd2 - eps2*over4*logR2e
      A5       = two*eps2*over2*logR2e - one - eps2*rd2
      A6       = rd4 - eps2*over2*rd4 + two*over6*eps2*logR2e
      A7       = rd4 + two*over6*eps2*logR2e
      A8       = -eps2*over4*rd2 - eps2*over2*rd4 + rd4 + two*over6*eps2*logR2e
      Ax1      = over2*logR2e - rd2 + four*over2*rd2*dx2 + two*rd4*dx2 - four*over4*dx2*logR2e
      Ax2      = rd4 - two*eps2*over4*rd2 + two*over6*eps2*logR2e - eps2*over2*rd4
      Ay1      = over2*logR2e - rd2 + four*over2*rd2*dy2 + two*rd4*dy2 - four*over4*dy2*logR2e
      !Ay2      = eps2*over2*rd2 - four*eps2*dy2*over4*rd2 - rd2 - over4*eps2*logR2e
      Ay3      = eps2*over2*rd2 - four*eps2*dy2*over4*rd2 - two*over2*rd4*eps2*dy2 - rd2 - over4*eps2*logR2e

      Axy(1) = half *t%monoj(1)*A0 - one *dx*over2*A1*vdotxM                                                &    !  monopole
               + one *dx*A4*t%dipjx(1) +  one *dy*A4*t%dipjy(1)                                             &    !  dipole I               V
               + two *over4*dx2*A5*vdotxDx - one *over2*A1*( vdotxDx + t%dipjx(1)*dx )                      &    !  - x dipole II          V
               - two *over4*eps2*dx*dy*vdotxDy*A2 + one *over2*dx*A1*( two*over2*dy*vdotxDy - t%dipjy(2) )  &    !  - y dipole II          V
               + two *dx2*t%quadjx(1)*Ax2 + one *t%quadjx(1)*A4                                             &    !  - xx quadrupole I      V
               + two *t%quadjy(1)*dy2*Ax2 + one *t%quadjy(1)*A4                                             &    !  - yy quadrupole I      V
               + two *t%quadjxy(1)*dx*dy*Ax2                                                                &    !  - xy quadrupole I      V
               + four *over4*t%quadjx(1)*dx2*A1 - four *over4*dx*eps2*A2*( vdotxQxx + t%quadjx(1)*dx )      &    !  - xx quadrupole II     V
               - eight *over6*dx3*A1*vdotxQxx - two *over2*A1*t%quadjx(1)                                   &
               + six *over4*dx*A1*vdotxQxx + eight *over6*dx3*eps2*A2*vdotxQxx                              &
               + two *eps2*over4*dx*Ax1*vdotxQxx                                                            &
               + two *over4*dx*A1*vdotxQyy + two *over4*dx*eps2*vdotxQyy*Ay1                                &    !  - yy quadrupole II		V
               - four *over4*dx*dy*A2*eps2*t%quadjy(2) - eight *over6*dx*dy2*A1*vdotxQyy                    &
               + four *over4*dx*dy*A1*t%quadjy(2) + eight *over6*eps2*dx*dy2*A2*vdotxQyy                    &
               + two *over4*dx2*A1*t%quadjxy(2) - two *over4*eps2*dy*A2*vdotxQxy - over2 *A1*t%quadjxy(2)   &    !  - xy quadrupole II		V
               + two *over4*dy*A1*vdotxQxy +four *over6*eps2*dx2*dy*A2*vdotxQxy                             &
               + four *over4*eps2*dx2*dy*A3*vdotxQxy - two *over4*eps2*dx2*A2*t%quadjxy(2)                  &
               - two *over4*eps2*dx*dy*A2*t%quadjxy(1) - eight *over6*dx2*dy*A1*vdotxQxy                    &
               + two *over4*dx*dy*A1*t%quadjxy(1) + four *over6*eps2*dx2*dy*A2*vdotxQxy
!
      Axy(2) = half *t%monoj(2)*A0   - one *dy*over2*A1*vdotxM                                                 &    !   monopole
               + one *A4*t%dipjx(2)*dx + one *A4*t%dipjy(2)*dy                                                 &    !   dipole I           V
               - two *eps2*dx*dy*over4*A2*vdotxDx + two *over4*dx*dy*A1*vdotxDx - one *over2*dy*A1*t%dipjx(1)  &    !  - x dipole II
               + two *over4*dy2*A5*vdotxDy - over2 *A1*( vdotxDy + t%dipjy(2)*dy)                              &    !  - y dipole II
               + two *dx2*t%quadjx(2)*Ax2 + one *t%quadjx(2)*A4                                                &    !  - xx quadrupole I   V
               + two *dy2*t%quadjy(2)*Ax2 + one *t%quadjy(2)*A4                                                &    !  - yy quadrupole I   V
               + two *t%quadjxy(2)*dx*dy*Ax2                                                                   &    !  - xy quadrupole I   V
               + two *over4*dy*A1*vdotxQxx + two *over4*dy*eps2*Ax1*vdotxQxx                                   &    !  - xx quadrupole II	V
               - four *eps2*over4*dx*dy*A2*t%quadjx(1) - eight *over6*dx2*dy*A1*vdotxQxx                       &
               + four *over4*dx*dy*A1*t%quadjx(1) + eight *over6*dx2*dy*eps2*A2*vdotxQxx                       &
               + four *over4*dy2*A1*t%quadjy(2) - four *over4*dy*eps2*A2*( vdotxQyy + t%quadjy(2)*dy )         &    !  - yy quadrupole II	V
               - eight *over6*dy3*A1*vdotxQyy - two *over2*A1*t%quadjy(2)                                      &
               + six *over4*dy*A1*vdotxQyy + eight *over6*eps2*dy3*A2*vdotxQyy                                 &
               + two *over4*eps2*dy*Ay1*vdotxQyy                                                               &
               + two *over4*dy2*A1*t%quadjxy(1) - two *over4*eps2*dx*A2*vdotxQxy - over2 *A1*t%quadjxy(1)      &    !  - xy quadrupole II	V
               + two *over4*dx*A1*vdotxQxy + four *over6*eps2*dx*dy2*A2*vdotxQxy                               &
               + four *over4*eps2*dx*dy2*A3*vdotxQxy - two *over4*eps2*dx*dy*A2*t%quadjxy(2)                   &
               - two *over4*eps2*dy2*A2*t%quadjxy(1) - eight *over6*dx*dy2*A1*vdotxQxy                         &
               + two *over4*dx*dy*A1*t%quadjxy(2) + four *over6*eps2*dx*dy2*A2*vdotxQxy


      phi           = half*phi
      exy           = half*exy
      Jxy           = half*Jxy
      Bz            = half*Bz
      Jirrxy        = half*Jirrxy

      Axy           = half*Axy
      
      Axy           =  Axy/vtilde  
      Bz            =  Bz/vtilde


    end subroutine calc_force_darwin_2D
    
    subroutine calc_force_darwin_2D3V(t, d, d2,eps2, exy, Axy, Jxy, Jirrxy, Bxy, phi, gradxA, gradyA)
      use module_globals, only: vtilde 
      use module_tool   , only: cross_product,double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(1:2),Axy(1:3),Jxy(1:3),Jirrxy(1:3),Bxy(1:3),phi, gradxA(1:3), gradyA(1:3)

      real(kind_physics) :: dx,dy,rd2,rd4,rd6,rd8,dx2,dy2,dx3,dy3,r2,logR2e,over2,over4,over6,over8, &
                            phi_mono,phi_dip(1:2),phi_quad(1:3),rho_mono,rho_dip(1:2),rho_quad(1:3), &
                            E_mono(1:3),Ex_dip(1:3),Ey_dip(1:3),Exx_quad(1:3),Eyy_quad(1:3),logTmp,  &
                            A1_mono,A1_dip(1:2),A1_quad(1:3),A2_mono,A2_dip(1:2),A2_quad(1:3),       &
                            A3_mono(1:3),A3x_dip(1:3),A3y_dip(1:3),A3xx_quad(1:3),A3yy_quad(1:3),    &   
                            dist(1:3),Jir_mono,Jir_dip(1:2),Jir_quad(1:3),over,Exy_quad(1:3),        &   
                            A3xy_quad(1:3),A1_oct(1:4),A3xxx_oct(1:3),A3xxy_oct(1:3)                 &
                            ,A3xyy_oct(1:3),A3yyy_oct(1:3) 

      dx = d(1)
      dy = d(2)
      dist(1:2) = d(1:2)
      dist(3)   = zero
      r2 =  dx**2 + dy**2
      over  = sqrt(one/r2)
      over2 = one/r2
      over4 = over2*over2
      over6 = over4*over2
      over8 = over6*over2
      rd2 = one/d2
      rd4 = rd2*rd2
      rd6 = rd4*rd2
      rd8 = rd6*rd2

      dx2 = dx *dx
      dy2 = dy *dy
      dx3 = dx2*dx
      dy3 = dy2*dy
      
      logTmp     = log(r2/eps2 + one)
      logR2e     = eps2*over2*logTmp 
      
      !!! Rho/J must be multiply by eps2/pi 
      rho_mono   =  rd4
      rho_dip(1:2)= -four*d(1:2)*rd6
      rho_quad(1)=  twentyfour*dx2*rd8 - four*rd6
      rho_quad(2)=  twentyfour*dy2*rd8 - four*rd6
      rho_quad(3)=  twentyfour*dy*dx*rd8 
      
      !!! Jirr must be multiply by 2 
      Jir_mono      = rd2
      Jir_dip(1:2)  = two*d(1:2)*rd4
      Jir_quad(1)   = two*rd4*( four*dx2*rd2 - one )
      Jir_quad(2)   = two*rd4*( four*dy2*rd2 - one )
      Jir_quad(3)   = eight*dx*dy*rd6
      
      phi_mono   = -log(d2)
      phi_dip    = -two*d(1:2)*rd2
      phi_quad(1)=  four*dx2*rd4 - two*rd2
      phi_quad(2)=  four*dy2*rd4 - two*rd2
      phi_quad(3)=  four*dy*dx*rd4 

      E_mono(1:2)= -phi_dip
      E_mono(3)  =  zero
      Ex_dip(1:3)= -(/phi_quad(1),phi_quad(3),zero /)
      Ey_dip(1:3)= -(/phi_quad(3),phi_quad(2),zero /)
      
      Exx_quad(1)=  sixteen*dx3*rd6    - twelve*dx*rd4 
      Exx_quad(2)=  sixteen*dx2*dy*rd6 - four*dy*rd4  !! This is E_xy(1)
      Exx_quad(3)=  zero
      Eyy_quad(1)=  sixteen*dx*dy2*rd6 - four*dx*rd4  !! This is E_xy(2)  
      Eyy_quad(2)=  sixteen*dy3*rd6    - twelve*dy*rd4 
      Eyy_quad(3)=  zero
      
      Exy_quad(1)=  Exx_quad(2)
      Exy_quad(2)=  Eyy_quad(1)
      Exy_quad(3)=  zero
      
      A1_mono    =  logR2e - one
      A1_dip(1:2)=  two*eps2*d(1:2)*rd2*over2 - two*over2*logR2e*d(1:2)
      A1_quad(1) =  two*eps2*over2*( rd2      - logTmp*over2 - four*dx2*rd2*over2 &
                    - two*dx2*rd4 + four*logTmp*dx2*over4 )
      A1_quad(2) =  two*eps2*over2*( rd2      - logTmp*over2 - four*dy2*rd2*over2 &
                    - two*dy2*rd4 + four*logTmp*dy2*over4 )
      A1_quad(3) =  four*over2*dx*dy*( two*logR2e*over2 - eps2*rd4 - two*eps2*rd2*over2 )
      
      A1_oct(1)  =  twelve*eps2*dx*over2*( four*dx2*rd2*over4 - two*rd2*over2 - rd4 +&
                    two*logTmp*over4 + two*dx2*rd4*over2 - four*dx2*logTmp*over6  ) + sixteen*eps2*dx3*rd6*over2 ! xxx
                    
      A1_oct(2)  =  four*eps2*dy*over2*( two*logTmp*over4 - rd4 - two*rd2*over2 + twelve*dx2*rd2*over4 +&
                    six*dx2*rd4*over2 + four*dx2*rd6 - twelve*dx2*logTmp*over6  )                                ! xxy
                    
      A1_oct(3)  =  four*eps2*dx*over2*( two*logTmp*over4 - rd4 - two*rd2*over2 + twelve*dy2*rd2*over4 +&
                    six*dy2*rd4*over2 + four*dy2*rd6 - twelve*dy2*logTmp*over6  )                                ! xxy
                    
      A1_oct(4)  =  twelve*eps2*dy*over2*( four*dy2*rd2*over4 - two*rd2*over2 - rd4 +&
                    two*logTmp*over4 + two*dy2*rd4*over2 - four*dy2*logTmp*over6  ) + sixteen*eps2*dy3*rd6*over2 ! yyy             
                    
      A2_mono    = -A1_mono + one + phi_mono
      A2_dip(1:2)= -A1_dip(1:2)   + phi_dip(1:2)
      A2_quad(1:3)=-A1_quad(1:3)  + phi_quad(1:3)
      
      A3_mono(1:3)  = over*dist(1:3)
      A3x_dip(1)    = over - dx2*over*over2
      A3x_dip(2)    =      - dx*dy*over*over2
      A3x_dip(3)    = zero
      A3y_dip(1)    =      - dx*dy*over*over2 
      A3y_dip(2)    = over - dy2*over*over2
      A3y_dip(3)    = zero
      A3xx_quad(1)  =-three*over*over2*dx*( one - dx2*over2 )
      A3xx_quad(2)  =-over*over2*dy*( one - three*dx2*over2 )   ! This is equal to A3_xy(1)
      A3xx_quad(3)  = zero
      A3yy_quad(1)  =-over*over2*dx*( one - three*dy2*over2 )     ! This is equal to A3_xy(2)  
      A3yy_quad(2)  =-three*over*over2*dy*( one - dy2*over2 ) 
      A3yy_quad(3)  = zero
      A3xy_quad(1)  = A3xx_quad(2)
      A3xy_quad(2)  = A3yy_quad(1)
      A3xy_quad(3)  = zero
      A3xxx_oct(1)  = three*over**3*( six*dx2*over2 - one - five*dx2**2*over4  )
      A3xxx_oct(2)  = three*over**5*dx*dy*( three - five*dx2*over2  )
      A3xxx_oct(3)  = zero
      
      A3xxy_oct(1)  = three*over**5*dx*dy*( three - five*dx2*over2  )
      A3xxy_oct(2)  = over**3*( two - fifteen*dx2*dy2*over4  )
      A3xxy_oct(3)  = zero
      
      A3xyy_oct(1)  = over**3*( two - fifteen*dx2*dy2*over4  )
      A3xyy_oct(2)  = three*over**5*dx*dy*( three - five*dy2*over2  )
      A3xyy_oct(3)  = zero
      
      A3yyy_oct(1)  = three*over**5*dx*dy*( three - five*dy2*over2  )
      A3yyy_oct(2)  = three*over**3*( six*dy2*over2 - one - five*dy2**2*over4  )
      A3yyy_oct(3)  = zero

      
      phi        = ( t%charge*phi_mono      )                                                              &! Monopole
                  +( dot_product(t%dip(1:2),phi_dip(1:2) ) )     !                                          &! Dipole
!                  +( dot_product(t%quad(1:2),phi_quad(1:2) ) )                                             &! Quadrupole
!                  +( phi_quad(3)*t%xyquad )                                                                 ! Quadrupole
      phi        = half*phi
                   
      exy(1:2)   = ( t%charge*E_mono(1:2)   )                                                              &! Monopole
                  +( t%dip(1)*Ex_dip(1:2) + t%dip(2)*Ey_dip(1:2)  )                         !               &! Dipole
!                  +( t%quad(1)*Exx_quad(1:2) + t%quad(2)*Eyy_quad(1:2) )                                   &! Quadrupole
!                  +( t%xyquad*Exy_quad(1:2) )                                                               ! Quadrupole 
      
      exy(1:2)   = half*exy(1:2)
      
      Jxy(1:3)   = ( t%monoj(1:3)*rho_mono )                                                               &! Monopole
                  +( rho_dip(1)*t%dipjx(1:3) + rho_dip(2)*t%dipjy(1:3)  ) !                                 &! Dipole
!                  +( t%quadjx(1:3)*rho_quad(1) + t%quadjy(1:3)*rho_quad(2) )                               &! Quadrupole
!                  +( t%quadjxy(1:3) *rho_quad(3) )                                                          ! Quadrupole 
                  
      Jxy(1:3)   = half*eps2/pi*Jxy(1:3)
                  
                  
      Jirrxy(1:3)= ( t%monoj(1:3)*(three*half*Jir_mono - eps2*rho_mono)                                     &
                   + double_cross_product_left( E_mono(1:3), E_mono(1:3), t%monoj(1:3)  ) )                 &! Monopole
                 + ( t%dipjx(1:3)*( three*half*Jir_dip(1)   - eps2*rho_dip(1) ) )                           &! Dipole - x 
                 + ( double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)  ) )               &
                 + ( double_cross_product_left( E_mono(1:3), Ex_dip(1:3), t%dipjx(1:3)  ) )                 &
                 + ( t%dipjy(1:3)*( three*half*Jir_dip(2)   - eps2*rho_dip(2) ) )                           &! Dipole - y 
                 + ( double_cross_product_left( Ey_dip(1:3), E_mono(1:3), t%dipjy(1:3)  ) )                 &
                 + ( double_cross_product_left( E_mono(1:3), Ey_dip(1:3), t%dipjy(1:3)  ) )        !         &
!                 + ( t%quadjx(1:3)*( three*half*Jir_quad(1)   - eps2*rho_quad(1) ) )                        &! Quadrupole - xx 
!                 + ( double_cross_product_left( Exx_quad(1:3), E_mono(1:3), t%quadjx(1:3)  ) )              &
!                 + ( double_cross_product_left( E_mono(1:3), Exx_quad(1:3), t%quadjx(1:3)  ) )              &
!                 + ( two*double_cross_product_left( Ex_dip(1:3), Ex_dip(1:3), t%quadjx(1:3)  ) )            &
!                 + ( t%quadjx(1:3)*( three*half*Jir_quad(2)   - eps2*rho_quad(2) ) )                        &! Quadrupole - yy 
!                 + ( double_cross_product_left( Eyy_quad(1:3), E_mono(1:3), t%quadjy(1:3)  ) )              &
!                 + ( double_cross_product_left( E_mono(1:3), Eyy_quad(1:3), t%quadjy(1:3)  ) )              &
!                 + ( two*double_cross_product_left( Ey_dip(1:3), Ey_dip(1:3), t%quadjy(1:3)  ) )            &
!                 + ( t%quadjxy(1:3)*( three*half*Jir_quad(3)   - eps2*rho_quad(3) ) )                       &! Quadrupole - xy 
!                 + ( double_cross_product_left( Exy_quad(1:3), E_mono(1:3), t%quadjxy(1:3)  ) )             &
!                 + ( double_cross_product_left( E_mono(1:3), Exy_quad(1:3), t%quadjxy(1:3)  ) )             &
!                 + ( double_cross_product_left( Ex_dip(1:3), Ey_dip(1:3), t%quadjxy(1:3)  ) )               &
!                 + ( double_cross_product_left( Ey_dip(1:3), Ex_dip(1:3), t%quadjxy(1:3)  ) ) 
                 
      Jirrxy(1:3)= half*oneoverpi*Jirrxy(1:3)
                 
      Bxy(1:3)   = cross_product(t%monoj(1:3), E_mono(1:3) )                                                &! Monopole
                  +( cross_product(t%dipjx(1:3),Ex_dip(1:3) )   )                                           &! Dipole x
                  +( cross_product(t%dipjy(1:3),Ey_dip(1:3) )   )            !                               &! Dipole y
!                  +( cross_product(t%quadjx(1:3),Exx_quad(1:3) ) )                                          &! Quadrupole xx
!                  +( cross_product(t%quadjy(1:3),Eyy_quad(1:3) ) )                                          &! Quadrupole yy
!                  +( cross_product(t%quadjxy(1:3),Exy_quad(1:3) ) )                                          ! Quadrupole xy
                  
      Bxy(1:3)   = half/vtilde*Bxy(1:3)  
      
      Axy(1:3)   = half*t%monoj(1:3)*( two - logR2e - log(d2)  )                                            &
                 + A1_mono*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )           &! Monopole
                 + half*( t%dipjx(1:3)*( -A1_dip(1) + phi_dip(1) ) )                                        &! Dipole - x 
                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3),A3x_dip(1:3)  ) )        &
                 + A1_dip(1)*( double_cross_product_left(A3_mono(1:3), t%dipjx(1:3),A3_mono(1:3)  ) )       &
                 + half*( t%dipjy(1:3)*( -A1_dip(2) + phi_dip(2) ) )                                        &! Dipole - y 
                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3),A3x_dip(1:3)  ) )        &
                 + A1_dip(2)*( double_cross_product_left(A3_mono(1:3), t%dipjy(1:3),A3_mono(1:3)  ) ) !      &
!                 + half*( t%quadjx(1:3)*( -A1_quad(1) + phi_quad(1) ) )                                     &! Quadrupole - xx 
!                 + A1_mono*( double_cross_product_left( A3xx_quad(1:3), t%quadjx(1:3), A3_mono(1:3)  ) )    &
!                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%quadjx(1:3), A3xx_quad(1:3)  ) )    &
!                 + A1_mono*( two*double_cross_product_left( A3x_dip(1:3), t%quadjx(1:3), A3x_dip(1:3)  ) )  & 
!!                 
!                 + two*A1_dip(1)*( double_cross_product_left( A3x_dip(1:3), t%quadjx(1:3), A3_mono(1:3)  ) )&
!                 + two*A1_dip(1)*( double_cross_product_left( A3_mono(1:3), t%quadjx(1:3), A3x_dip(1:3)  ) )&
!                 + A1_quad(1)*double_cross_product_left( A3_mono(1:3), t%quadjx(1:3), A3_mono(1:3)  )       &
!                 
!                 + half*( t%quadjx(1:3)*( -A1_quad(2) + phi_quad(2) ) )                                     &! Quadrupole - yy
!                 + A1_mono*( double_cross_product_left( A3yy_quad(1:3), t%quadjy(1:3), A3_mono(1:3)  ) )    &
!                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%quadjy(1:3), A3yy_quad(1:3)  ) )    &
!                 + A1_mono*( two*double_cross_product_left( A3y_dip(1:3), t%quadjy(1:3), A3y_dip(1:3)  ) )  &
!                 
!                 + two*A1_dip(2)*( double_cross_product_left( A3y_dip(1:3), t%quadjy(1:3), A3_mono(1:3)  ) )&
!                 + two*A1_dip(2)*( double_cross_product_left( A3_mono(1:3), t%quadjy(1:3), A3y_dip(1:3)  ) )&
!                 + A1_quad(2)*double_cross_product_left( A3_mono(1:3), t%quadjy(1:3), A3_mono(1:3)  )       &
!                 
!                 + half*( t%quadjxy(1:3)*( -A1_quad(3) + phi_quad(3) ) )                                    &! Quadrupole - xy
!                 + A1_mono*( double_cross_product_left( A3xy_quad(1:3), t%quadjxy(1:3), A3_mono(1:3)  ) )   &
!                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%quadjxy(1:3), A3xy_quad(1:3)  ) )   &
!                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%quadjxy(1:3), A3y_dip(1:3)  ) )     &
!                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%quadjxy(1:3), A3x_dip(1:3)  ) )     &
!!                 
!                 + A1_dip(1)*( double_cross_product_left( A3y_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)  ) )  &
!                 + A1_dip(1)*( double_cross_product_left( A3_mono(1:3), t%quadjxy(1:3), A3y_dip(1:3)  ) )  &
!                 + A1_dip(2)*( double_cross_product_left( A3x_dip(1:3), t%quadjxy(1:3), A3_mono(1:3)  ) )  &
!                 + A1_dip(2)*( double_cross_product_left( A3_mono(1:3), t%quadjxy(1:3), A3x_dip(1:3)  ) )  &
!                 + A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%quadjxy(1:3), A3_mono(1:3)  )     
!                 
      Axy(3)     = t%monoj(3)*phi_mono                                                                     &! Monopole
                  +  t%dipjx(3)*phi_dip(1)   + t%dipjy(3)*phi_dip(2)                   !                    &! Dipole
!                  +  t%quadjx(3)*phi_quad(1) + t%quadjy(3)*phi_quad(2)                                     &! Quadrupole
!                  +  t%quadjxy(3)*phi_quad(3)                                                               ! Quadrupole 
                  
      Axy(1:3)   = half/vtilde*Axy(1:3)
      
      
                 
      gradxA(1:3)=  half*( t%monoj(1:3)*( -A1_dip(1) + phi_dip(1) ) )                                       &! Monopole 
                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%monoj(1:3),A3x_dip(1:3)  ) )        &
                 + A1_dip(1)*( double_cross_product_left(A3_mono(1:3), t%monoj(1:3),A3_mono(1:3)  ) )       &
                 + half*( t%dipjx(1:3)*( -A1_quad(1) + phi_quad(1) ) )                                      &! Dipole - x 
                 + A1_mono*( double_cross_product_left( A3xx_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )     &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3xx_quad(1:3)  ) )     &
                 + A1_mono*( two*double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)  ) )   & 
!                 
                 + two*A1_dip(1)*( double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )&
                 + two*A1_dip(1)*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)  ) )&
                 + A1_quad(1)*double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3_mono(1:3)  )       &
                 
                 + half*( t%dipjy(1:3)*( -A1_quad(3) + phi_quad(3) ) )                                    &! Dipole - y
                 + A1_mono*( double_cross_product_left( A3xy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )   &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3xy_quad(1:3)  ) )   &
                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)  ) )     &
                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3x_dip(1:3)  ) )     &
!                 
                 + A1_dip(1)*( double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )  &
                 + A1_dip(1)*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)  ) )  &
                 + A1_dip(2)*( double_cross_product_left( A3x_dip(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )  &
                 + A1_dip(2)*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3x_dip(1:3)  ) )  &
                 + A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%dipjy(1:3), A3_mono(1:3)  )     
                 
      gradxA(3)  = t%monoj(3)*phi_dip(1)                                                                    &! Monopole
                   + t%dipjx(3)*phi_quad(1)   + t%dipjy(3)*phi_quad(3)       !                               &! Dipole
!                   -  t%quadjx(3)*Exx_quad(1) - t%quadjy(3)*Eyy_quad(1)                                     &! Quadrupole
!                   -  t%quadjxy(3)*Exx_quad(2)                                                               ! Quadrupole 
                  
      gradyA(1:3)= half*( t%monoj(1:3)*( -A1_dip(2) + phi_dip(2) ) )                                        &! Monopole 
                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%monoj(1:3),A3y_dip(1:3)  ) )        &
                 + A1_dip(2)*( double_cross_product_left(A3_mono(1:3), t%monoj(1:3),A3_mono(1:3)  ) )       &
                 
                 + half*( t%dipjx(1:3)*( -A1_quad(3) + phi_quad(3) ) )                                    &! Dipole - x
                 + A1_mono*( double_cross_product_left( A3xy_quad(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )   &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3xy_quad(1:3)  ) )   &
                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3y_dip(1:3)  ) )     &
                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%dipjx(1:3), A3x_dip(1:3)  ) )     &
!                 
                 + A1_dip(1)*( double_cross_product_left( A3y_dip(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )  &
                 + A1_dip(1)*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3y_dip(1:3)  ) )  &
                 + A1_dip(2)*( double_cross_product_left( A3x_dip(1:3), t%dipjx(1:3), A3_mono(1:3)  ) )  &
                 + A1_dip(2)*( double_cross_product_left( A3_mono(1:3), t%dipjx(1:3), A3x_dip(1:3)  ) )  &
                 + A1_quad(3)*double_cross_product_left( A3_mono(1:3) , t%dipjx(1:3), A3_mono(1:3)  )    &
                 
                 + half*( t%dipjy(1:3)*( -A1_quad(2) + phi_quad(2) ) )                                     &! Dipole - y
                 + A1_mono*( double_cross_product_left( A3yy_quad(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )    &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3yy_quad(1:3)  ) )    &
                 + A1_mono*( two*double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3y_dip(1:3)  ) )  &
                 
                 + two*A1_dip(2)*( double_cross_product_left( A3y_dip(1:3), t%dipjy(1:3), A3_mono(1:3)  ) )&
                 + two*A1_dip(2)*( double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3y_dip(1:3)  ) )&
                 + A1_quad(2)*double_cross_product_left( A3_mono(1:3), t%dipjy(1:3), A3_mono(1:3)  )                    
      
      gradyA(3)  = t%monoj(3)*phi_dip(2)                                                                     &! Monopole
                   + t%dipjx(3)*phi_quad(2)  + t%dipjy(3)*phi_quad(3)                  !                      &! Dipole
!                   - t%quadjx(3)*Exx_quad(2) - t%quadjy(3)*Eyy_quad(2)                                       &! Quadrupole
!                   - t%quadjxy(3)*Eyy_quad(1)                                                                 ! Quadrupole 
                   
      gradxA(1:3)= half/vtilde*gradxA(1:3)
      gradyA(1:3)= half/vtilde*gradyA(1:3)
      
      
    end subroutine calc_force_darwin_2D3V

!    subroutine calc_force_darwin_2D3V(t, d, d2,eps2, exy, Axy, Jxy, Jirrxy, Bxy, phi)
!      use module_globals, only: vtilde  
!      implicit none
!
!      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
!      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
!      real(kind_physics), intent(out) ::  exy(1:2),Axy(1:3),Jxy(1:3),Jirrxy(1:3),Bxy(1:3),phi
!
!      real(kind_physics) :: dx,dy,rd2,rd4,rd6,rd8,dx2,dy2,dx3,dy3,r2,logR2e,over2,over4,over6,over8
!      real(kind_physics) :: Ax1,Ax2,Ay1,Ay2,Ay3,A0,A1,A2,A3,A4,A5,A6,A7,A8,vdotxM,vdotxDx,vdotxDy
!      real(kind_physics) :: vdotxQxx,vdotxQxy,vdotxQyy, vcrossxM, vcrossxDx,vcrossxDy,vcrossxQx,vcrossxQxy,vcrossxQy
!      real(kind_physics) :: rho0,rho1,rho2,rho3,phix1,phiy1,phix2,phiy2,phixy,E1,E2,E3,E4,Jirr_dipx(1:2),Jirr_dipy(1:2)
!      real(kind_physics) :: Jirr_quadx(1:2),Jirr_quady(1:2),Jirr_quadxy(1:2)
!
!
!
!      dx = d(1)
!      dy = d(2)
!
!      r2 =  dx**2 + dy**2
!      over2 = one/r2
!      over4 = over2*over2
!      over6 = over4*over2
!      over8 = over6*over2
!      rd2 = one/d2
!      rd4 = rd2*rd2
!      rd6 = rd4*rd2
!      rd8 = rd6*rd2
!
!      dx2 = dx *dx
!      dy2 = dy *dy
!      dx3 = dx2*dx
!      dy3 = dy2*dy
!
!      phix1 = - two*rd2*dx
!      phix2 = two*( - rd2 + two*rd4*dx2 )
!      phiy1 = - two*rd2*dy
!      phiy2 = two*( - rd2 + two*rd4*dy2 )
!      phixy = four*rd4*dx*dy
!
!      phi = - t%charge * log(d2)                                          &     !  monopole term
!            + phix1*t%dip(1) + phiy1*t%dip(2)                             &     !  dipole
!            + phix2*t%quad(1) + phiy2*t%quad(2) + phixy*t%xyquad                !  quadrupole
!
!
!      E1     =  two*rd4*( - six + eight*rd2*dx2  )
!      E2     =  two*rd4*( - six + eight*rd2*dy2  )
!      E3     =  two*rd4*( - two + eight*rd2*dx2  )
!      E4     =  two*rd4*( - two + eight*rd2*dy2  )
!
!
!      exy(1) = - phix1*t%charge                                          &     ! monopole
!               - phix2*t%dip(1) - phixy*t%dip(2)                         &     ! dipole
!               + dx*E1*t%quad(1) + dx*E4*t%quad(2) + dy*E3*t%xyquad            ! quadrupole
!
!
!      exy(2) = - phiy1*t%charge                                          &     ! monopole
!               - phixy*t%dip(1) - phiy2*t%dip(2)                         &     ! dipole
!               + dy*E3*t%quad(1) + dy*E2*t%quad(2) + dx*E4*t%xyquad            ! quadrupole
!
!      rho0   = eps2/pi
!      rho1   = four*rho0*rd6*( six*dx2*rd2 - one )
!      rho2   = four*rho0*rd6*( six*dy2*rd2 - one )
!      rho3   = twentyfour*rho0*rd8*dx*dy
!
!!      rho    = rho0*rd4*t%charge                                 &     ! monopole
!!             - four*rho0*rd6*( t%dip(1)*dx + t%dip(2)*dy )       &     ! dipole
!!             + rho1*t%quad(1) + rho2*t%quad(2) + rho3*t%xyquad         ! quadrupole
!
!      Jxy(1) = rho0*rd4*t%monoj(1)                                      &     ! monopole
!             - four*rho0*rd6*( t%dipjx(1)*dx + t%dipjy(1)*dy  )         &     ! dipole
!             + rho1*t%quadjx(1) + rho2*t%quadjy(1) + rho3*t%quadjxy(1)        ! quadrupole
!
!
!      Jxy(2) = rho0*rd4*t%monoj(2)                                      &     ! monopole
!             - four*rho0*rd6*( t%dipjx(2)*dx + t%dipjy(2)*dy  )         &     ! dipole
!             + rho1*t%quadjx(2) + rho2*t%quadjy(2) + rho3*t%quadjxy(2)        ! quadrupole
!
!      Jxy(3) = rho0*rd4*t%monoj(3)                                      &     ! monopole
!             - four*rho0*rd6*( t%dipjx(3)*dx + t%dipjy(3)*dy  )         &     ! dipole
!             + rho1*t%quadjx(3) + rho2*t%quadjy(3) + rho3*t%quadjxy(3)        ! quadrupole
!             
!
!      Bxy(1)     = -rd2*t%monoj(3)*dy                                                                    &   ! - monopole
!                   +two*rd4*dx*dy*t%dipjx(3) + two*rd4*dy2*t%dipjy(3)  -  rd2*t%dipjy(3)                 &   ! - x dipole - y dipole
!                   +two*rd4*dy*t%quadjx(3)   - eight*rd6*dx2*dy*t%quadjx(3)                              &   ! - xx  quadrupole
!                   +six*rd4*dy*t%quadjy(3)   - eight*rd6*dy3*t%quadjy(3)                                 &   ! - yy  quadrupole
!                   +two*rd4*t%quadjxy(3)*dx  - eight*rd6*dx*dy2*t%quadjxy(3)                                 ! - xy  quadrupole
!
!
!      Bxy(2)     =  rd2*t%monoj(3)*dx                                                                    &   ! - monopole
!                   -two*rd4*dx2*t%dipjx(3)   + rd2*t%dipjx(3) - two*rd4*dx*dy*t%dipjy(3)                 &   ! - x dipole - y dipole
!                   -six*rd4*dx*t%quadjx(3)   + eight*rd6*dx3*t%quadjx(3)                                 &   ! - xx  quadrupole
!                   -two*rd4*dx*t%quadjy(3)   - eight*rd6*dx*dy2*t%quadjy(3)                              &   ! - yy  quadrupole
!                   -two*rd4*dy*t%quadjxy(3)  + eight*rd6*dx2*dy*t%quadjxy(3)                                 ! - xy  quadrupole
!
!      Bxy(1:2)   = two *Bxy(1:2)
!
!
!      vcrossxM   = t%monoj(1)*dy - t%monoj(2)*dx
!      vcrossxDx  = t%dipjx(1)*dy - t%dipjx(2)*dx
!      vcrossxDy  = t%dipjy(1)*dy - t%dipjy(2)*dx
!      vcrossxQx  = t%quadjx(1)*dy - t%quadjx(2)*dx
!      vcrossxQy  = t%quadjy(1)*dy - t%quadjy(2)*dx
!      vcrossxQxy = t%quadjxy(1)*dy - t%quadjxy(2)*dx
!
!      Bxy(3)     = two*rd2 *vcrossxM                                                                    &   ! - monopole
!                - two *rd2*t%dipjx(2) - four *rd4*dx*vcrossxDx                                          &   ! - x dipole
!                + two *rd2*t%dipjy(1) - four *rd4*dy*vcrossxDy                                          &   ! - y dipole
!                - four *rd4*( vcrossxQx - two*t%quadjx(2)*dx ) + four*four*rd6 *dx2*vcrossxQx           &   ! - xx  quadrupole
!                - four *rd4*( two*t%quadjy(1)*dy + vcrossxQy ) + four*four*rd6 *dy2*vcrossxQy           &   ! - yy  quadrupole
!                - four *rd4*( t%quadjxy(1)*dx - t%quadjxy(2)*dy ) + four*four *rd6*dx*dy*vcrossxQxy         ! - xy quadrupole
!
!
!      vdotxM   = t%monoj(1)*dx + t%monoj(2)*dy
!
!
!      logR2e = log( one + r2/eps2 )
!
!
!
!      vdotxM   = t%monoj(1)*dx  + t%monoj(2)*dy
!      vdotxDx  = t%dipjx(1)*dx  + t%dipjx(2)*dy
!      vdotxDy  = t%dipjy(1)*dx  + t%dipjy(2)*dy
!      vdotxQxx = t%quadjx(1)*dx + t%quadjx(2)*dy
!      vdotxQyy = t%quadjy(1)*dx + t%quadjy(2)*dy
!      vdotxQxy = t%quadjxy(1)*dx+ t%quadjxy(2)*dy
!
!
!      Jirr_dipx(1)    = eight*dx2*vdotxDx*rd6 - two*vdotxDx*rd4 - four*t%dipjx(1)*dx*rd4
!      Jirr_dipx(2)    = eight*dx*dy*vdotxDx*rd6 - two*t%dipjx(2)*dx*rd4 - two*t%dipjx(1)*dy*rd4
!      Jirr_dipy(1)    = eight*dx*dy*vdotxDy*rd6 - two*t%dipjy(2)*dx*rd4 - two*t%dipjy(1)*dy*rd4
!      Jirr_dipy(2)    = eight*dy2*vdotxDy*rd6 - two*vdotxDy*rd4 - four*t%dipjy(2)*dy*rd4
!
!      Jirr_quadx(1)   = twentyfour*dx*vdotxQxx*rd6 - six*t%quadjx(1)*rd4 - two*twentyfour*dx3*vdotxQxx*rd8 + twentyfour*t%quadjx(1)*dx2*rd6
!      Jirr_quadx(2)   = eight*dy*vdotxQxx*rd6 - two*t%quadjx(2)*rd4 + eight*t%quadjx(2)*dx2*rd6 - two*twentyfour*dx2*dy*vdotxQxx*rd8 + four**2*t%quadjx(1)*dx*dy*rd6
!      Jirr_quady(1)   = eight*dx*vdotxQyy*rd6 - two*t%quadjy(1)*rd4 + eight*t%quadjy(1)*dy2*rd6 - two*twentyfour*dx*dy2*vdotxQyy*rd8 + four**2*t%quadjy(2)*dx*dy*rd6
!      Jirr_quady(2)   = twentyfour*dy*vdotxQyy*rd6 - six*t%quadjy(2)*rd4 - two*twentyfour*dy3*vdotxQyy*rd8 + twentyfour*t%quadjy(2)*dy2*rd6
!      Jirr_quadxy(1)  = eight*dy*vdotxQxy*rd6 - two*t%quadjxy(2)*rd4 + eight*t%quadjxy(2)*dx2*rd6 - two*twentyfour*dx2*dy*vdotxQxy*rd8 + four**2*t%quadjxy(1)*dx*dy*rd6
!      Jirr_quadxy(2)  = eight*dx*vdotxQxy*rd6 - two*t%quadjxy(1)*rd4 + eight*t%quadjxy(1)*dy2*rd6 - two*twentyfour*dx*dy2*vdotxQxy*rd8 + four**2*t%quadjxy(2)*dx*dy*rd6
!
!      Jirrxy(1:2)   = two*( ( t%monoj(1:2)*rd2 - two*d(1:2)*rd4*vdotxM )           &! Monopole
!                        + (  Jirr_dipx(1:2) + Jirr_dipy(1:2)  )                        &! Dipole
!                        + ( Jirr_quadx(1:2) + Jirr_quady(1:2) + Jirr_quadxy(1:2) )  )   ! Quadrupole
!
!      Jirrxy(3)     = zero
!
!
!      A0       = eps2*over2*logR2e  - log(d2)
!      A1       = eps2*over2*logR2e - one
!      A2       = rd2 - over2*logR2e
!      A3       = two*over2*rd2 + rd4 - two*over4*logR2e
!      A4       = eps2*over2*rd2 - rd2 - eps2*over4*logR2e
!      A5       = two*eps2*over2*logR2e - one - eps2*rd2
!      A6       = rd4 - eps2*over2*rd4 + two*over6*eps2*logR2e
!      A7       = rd4 + two*over6*eps2*logR2e
!      A8       = -eps2*over4*rd2 - eps2*over2*rd4 + rd4 + two*over6*eps2*logR2e
!      Ax1      = over2*logR2e - rd2 + four*over2*rd2*dx2 + two*rd4*dx2 - four*over4*dx2*logR2e
!      Ax2      = rd4 - two*eps2*over4*rd2 + two*over6*eps2*logR2e - eps2*over2*rd4
!      Ay1      = over2*logR2e - rd2 + four*over2*rd2*dy2 + two*rd4*dy2 - four*over4*dy2*logR2e
!      !Ay2      = eps2*over2*rd2 - four*eps2*dy2*over4*rd2 - rd2 - over4*eps2*logR2e
!      Ay3      = eps2*over2*rd2 - four*eps2*dy2*over4*rd2 - two*over2*rd4*eps2*dy2 - rd2 - over4*eps2*logR2e
!
!      Axy(1) = half *t%monoj(1)*A0 - one *dx*over2*A1*vdotxM                                               &    !  monopole
!              + one *dx*A4*t%dipjx(1) +  one *dy*A4*t%dipjy(1)                                             &    !  dipole I               V
!              + two *over4*dx2*A5*vdotxDx - one *over2*A1*( vdotxDx + t%dipjx(1)*dx )                      &    !  - x dipole II          V
!              - two *over4*eps2*dx*dy*vdotxDy*A2 + one *over2*dx*A1*( two*over2*dy*vdotxDy - t%dipjy(2) )  &    !  - y dipole II          V
!              + two *dx2*t%quadjx(1)*Ax2 + one *t%quadjx(1)*A4                                             &    !  - xx quadrupole I      V
!              + two *t%quadjy(1)*dy2*Ax2 + one *t%quadjy(1)*A4                                             &    !  - yy quadrupole I      V
!              + two *t%quadjxy(1)*dx*dy*Ax2                                                                &    !  - xy quadrupole I      V
!              + four *over4*t%quadjx(1)*dx2*A1 - four *over4*dx*eps2*A2*( vdotxQxx + t%quadjx(1)*dx )      &    !  - xx quadrupole II     V
!              - eight *over6*dx3*A1*vdotxQxx - two *over2*A1*t%quadjx(1)                                   &
!              + six *over4*dx*A1*vdotxQxx + eight *over6*dx3*eps2*A2*vdotxQxx                              &
!              + two *eps2*over4*dx*Ax1*vdotxQxx                                                            &
!              + two *over4*dx*A1*vdotxQyy + two *over4*dx*eps2*vdotxQyy*Ay1                                &    !  - yy quadrupole II    	V
!              - four *over4*dx*dy*A2*eps2*t%quadjy(2) - eight *over6*dx*dy2*A1*vdotxQyy                    &
!              + four *over4*dx*dy*A1*t%quadjy(2) + eight *over6*eps2*dx*dy2*A2*vdotxQyy                    &
!              + two *over4*dx2*A1*t%quadjxy(2) - two *over4*eps2*dy*A2*vdotxQxy - over2 *A1*t%quadjxy(2)   &    !  - xy quadrupole II		V
!              + two *over4*dy*A1*vdotxQxy +four *over6*eps2*dx2*dy*A2*vdotxQxy                             &
!              + four *over4*eps2*dx2*dy*A3*vdotxQxy - two *over4*eps2*dx2*A2*t%quadjxy(2)                  &
!              - two *over4*eps2*dx*dy*A2*t%quadjxy(1) - eight *over6*dx2*dy*A1*vdotxQxy                    &
!              + two *over4*dx*dy*A1*t%quadjxy(1) + four *over6*eps2*dx2*dy*A2*vdotxQxy
!!!
!      Axy(2) = half *t%monoj(2)*A0   - one *dy*over2*A1*vdotxM                                                &    !   monopole
!              + one *A4*t%dipjx(2)*dx + one *A4*t%dipjy(2)*dy                                                 &    !   dipole I           V
!              - two *eps2*dx*dy*over4*A2*vdotxDx + two *over4*dx*dy*A1*vdotxDx - one *over2*dy*A1*t%dipjx(1)  &    !  - x dipole II
!              + two *over4*dy2*A5*vdotxDy - over2 *A1*( vdotxDy + t%dipjy(2)*dy)                              &    !  - y dipole II
!              + two *dx2*t%quadjx(2)*Ax2 + one *t%quadjx(2)*A4                                                &    !  - xx quadrupole I   V
!              + two *dy2*t%quadjy(2)*Ax2 + one *t%quadjy(2)*A4                                                &    !  - yy quadrupole I   V
!              + two *t%quadjxy(2)*dx*dy*Ax2                                                                   &    !  - xy quadrupole I   V
!              + two *over4*dy*A1*vdotxQxx + two *over4*dy*eps2*Ax1*vdotxQxx                                   &    !  - xx quadrupole II	V
!              - four *eps2*over4*dx*dy*A2*t%quadjx(1) - eight *over6*dx2*dy*A1*vdotxQxx                       &
!              + four *over4*dx*dy*A1*t%quadjx(1) + eight *over6*dx2*dy*eps2*A2*vdotxQxx                       &
!              + four *over4*dy2*A1*t%quadjy(2) - four *over4*dy*eps2*A2*( vdotxQyy + t%quadjy(2)*dy )         &    !  - yy quadrupole II	V
!              - eight *over6*dy3*A1*vdotxQyy - two *over2*A1*t%quadjy(2)                                      &
!              + six *over4*dy*A1*vdotxQyy + eight *over6*eps2*dy3*A2*vdotxQyy                                 &
!              + two *over4*eps2*dy*Ay1*vdotxQyy                                                               &
!              + two *over4*dy2*A1*t%quadjxy(1) - two *over4*eps2*dx*A2*vdotxQxy - over2 *A1*t%quadjxy(1)      &    !  - xy quadrupole II	V
!              + two *over4*dx*A1*vdotxQxy + four *over6*eps2*dx*dy2*A2*vdotxQxy                               &
!              + four *over4*eps2*dx*dy2*A3*vdotxQxy - two *over4*eps2*dx*dy*A2*t%quadjxy(2)                   &
!              - two *over4*eps2*dy2*A2*t%quadjxy(1) - eight *over6*dx*dy2*A1*vdotxQxy                         &
!              + two *over4*dx*dy*A1*t%quadjxy(2) + four *over6*eps2*dx*dy2*A2*vdotxQxy
!
!
!      Axy(3) = -t%monoj(3) *log(d2)                                                        &    !   monopole
!               - two*rd2*t%dipjx(3)*(dx+dy)                                                &    !   dipole
!               + four*rd4*dx2*t%quadjx(3) - two*rd2*t%quadjx(3)                            &    !   quadrupole   - xx
!               + four*rd4*dy2*t%quadjy(3) - two*rd2*t%quadjy(3)                            &    !   quadrupole   - yy
!               + four*rd4*dx*dy*t%quadjxy(3)                                                    !   quadrupole   - xy
!
!
!
!      phi           = half*phi
!      exy           = half*exy
!      Jxy           = half*Jxy
!      Bxy           = half*Bxy
!      Jirrxy        = half*Jirrxy
!      Axy           = half*Axy
!      
!      Axy           =  Axy/vtilde  
!      Bxy           =  Bxy/vtilde
!
!    end subroutine calc_force_darwin_2D3V



    subroutine calc_force_darwin_3D(t, d, d2,eps2, exyz, Axyz, Jxyz, Jirrxyz, Bxyz, phi)
      use module_globals, only: vtilde  
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  ::  d(1:3), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(1:3),Axyz(1:3),Jxyz(1:3), Jirrxyz(1:3),Bxyz(1:3),phi

      real(kind_physics) :: dx,dy,dz,rd,rd2,rd3,rd4,rd5,rd6,rd7,rd8,rd9,dx2,dy2,dx3,dy3,r2,r,rEps
      real(kind_physics) :: vdotxM,vdotxDx,vdotxDy,dz2,dz3,fd1,fd2,fd3,fd4,fd5,fd6,m2rd5,m5rd7
      real(kind_physics) :: v_dot_x_mono,v_dot_x_dip(1:3),v_dot_x_quad(1:3),v_dot_x_quad_offdiag(3),pre1,pre2x,pre2y,pre2z,preQ1,preQ2,preQ3
      real(kind_physics) :: B_monopole(1:3),B_dip_x(1:3),B_dip_y(1:3),B_dip_z(1:3),B_quad_x(1:3),B_quad_y(1:3),B_quad_z(1:3),B_quad_xy(1:3),B_quad_yz(1:3),B_quad_zx(1:3)
      real(kind_physics) :: over,over2,over3,over4,over5,over6,over7,over8,over9,AI_dip,AI_quad(1:3),J_dip,J_quad(1:3),J_quad_offdiag,Jirr_quadI(1:3),Jirr_quad_offdiagI
      real(kind_physics) :: AII(1:3),AII_mono,AII1,AII2,AII3,AII4,AII5,AII6,AII7,AII8,AII_dipx(1:3),AII_dipy(1:3),AII_dipz(1:3),AII_quadx(1:3),AII_quady(1:3),AII_quadz(1:3)
      real(kind_physics) :: AII_quadxy(1:3),AII_quadyz(1:3),AII_quadzx(1:3),AII_quadrupole_mix,AII_dipole_mix_x,AII_dipole_mix_y,AII_dipole_mix_z
      real(kind_physics) :: AI_quad_offdiag,sinh_1,Jirr_dipI,Jirr_dipxII(1:3),Jirr_dipyII(1:3),Jirr_dipzII(1:3),Jirr_quadxII(1:3),Jirr_quadyII(1:3),Jirr_quadzII(1:3)
      real(kind_physics) :: Jirr_quadxyII(1:3),Jirr_quadyzII(1:3),Jirr_quadzxII(1:3),Jirr_mixedII,JirrI(1:3),JirrII(1:3)
!      pi = two*acos(zero)         !<----- pi
!      c = one                     !<----- speed of light

      dx = d(1)
      dy = d(2)
      dz = d(3)

      r2    = dx**2 + dy**2 + dz**2
      r     = sqrt(r2)
      rEps  = sqrt(d2)
      rd2   = one/d2 ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd    = sqrt(rd2)
      rd3   = rd *rd2
      rd5   = rd3*rd2
      rd7   = rd5*rd2
      rd9   = rd7*rd2
      over  = one/r
      over2 = over*over
      over3 = over2*over
      over4 = over3*over
      over5 = over4*over
      over6 = over5*over
      over7 = over6*over
      over8 = over7*over
      over9 = over8*over

      dx2 = dx*dx
      dy2 = dy*dy
      dz2 = dz*dz
      dx3 = dx*dx2
      dy3 = dy*dy2
      dz3 = dz*dz2

      fd1 = three*dx2*rd5 - rd3
      fd2 = three*dy2*rd5 - rd3
      fd3 = three*dz2*rd5 - rd3
      fd4 = three*dx*dy*rd5
      fd5 = three*dy*dz*rd5
      fd6 = three*dx*dz*rd5

      m2rd5 = two*rd5
      m5rd7 = five*rd7
      pre1 = m5rd7*dx*dy*dz
      pre2x = m5rd7*dx2 - rd5
      pre2y = m5rd7*dy2 - rd5
      pre2z = m5rd7*dz2 - rd5
      preQ1 = pre2x*t%quad(1)
      preQ2 = pre2y*t%quad(2)
      preQ3 = pre2z*t%quad(3)

      phi = t%charge*rd                                             &  !  monopole term
            - (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
            + fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3)         &  !  quadrupole
            + fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

      exyz(1) = t%charge*dx*rd3                                     &  ! monopole term
                - (fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3))      &  ! dipole term
                + three * (                                         &  ! quadrupole term
                   dx * (                                           &
                       ( pre2x - m2rd5 )*t%quad(1)                  &
                     + preQ2                                        &
                     + preQ3                                        &
                   )                                                &
                   + dy*pre2x*t%xyquad                              &
                   + dz*pre2x*t%zxquad                              &
                   + pre1*t%yzquad                                  &
                  )

      exyz(2) = t%charge*dy*rd3                                     &  !  monopole term
                - (fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3))      &  ! dipole term
                + three * (                                         &  ! quadrupole term
                   dy * (                                           &
                       ( pre2y - m2rd5 )*t%quad(2)                  &
                     + preQ1                                        &
                     + preQ3                                        &
                   )                                                &
                   + dx*pre2y*t%xyquad                              &
                   + dz*pre2y*t%yzquad                              &
                   + pre1*t%zxquad                                  &
                  )

      exyz(3) = t%charge*dz*rd3                                     &  !  monopole term
                - (fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1))      &  ! dipole term
                + three * (                                         &  ! quadrupole term
                   dz * (                                           &
                     + ( pre2z - m2rd5 )*t%quad(3)                  &
                     + preQ2                                        &
                     + preQ1                                        &
                   )                                                &
                   + dx*pre2z*t%zxquad                              &
                   + dy*pre2z*t%yzquad                              &
                   + pre1*t%xyquad                                  &
                  )

      J_dip               = -five*rd7
      J_quad(1:3)         = (/ thirtyfive*dx2*rd9 - five*rd7, thirtyfive*dy2*rd9 - five*rd7 , thirtyfive*dz2*rd9 - five*rd7 /)
      J_quad_offdiag      = thirtyfive*rd9

      Jxyz(1:3)           = quarter*oneoverpi*three*eps2*( rd5*t%monoj(1:3)                                                          & ! monopole
                           + J_dip*( t%dipjx(1:3)*dx + t%dipjy(1:3)*dy + t%dipjz(1:3)*dz )                                    & ! dipole
                           + t%quadjx(1:3)*J_quad(1) + t%quadjy(1:3)*J_quad(2) + t%quadjz(1:3)*J_quad(3)                              & ! diagonal quadrupole
                           + J_quad_offdiag*( t%quadjxy(1:3)*dx*dy + t%quadjyz(1:3)*dy*dz + t%quadjzx(1:3)*dx*dz )  )   ! off-diagonal quadrupole



      v_dot_x_mono             = t%monoj(1)*dx    + t%monoj(2)*dy    + t%monoj(3)*dz
      v_dot_x_dip(1)           = t%dipjx(1)*dx    + t%dipjx(2)*dy    + t%dipjx(3)*dz
      v_dot_x_dip(2)           = t%dipjy(1)*dx    + t%dipjy(2)*dy    + t%dipjy(3)*dz
      v_dot_x_dip(3)           = t%dipjz(1)*dx    + t%dipjz(2)*dy    + t%dipjz(3)*dz
      v_dot_x_quad(1)          = t%quadjx(1)*dx   + t%quadjx(2)*dy   + t%quadjx(3)*dz
      v_dot_x_quad(2)          = t%quadjy(1)*dx   + t%quadjy(2)*dy   + t%quadjy(3)*dz
      v_dot_x_quad(3)          = t%quadjz(1)*dx   + t%quadjz(2)*dy   + t%quadjz(3)*dz
      v_dot_x_quad_offdiag(1)  = t%quadjxy(1)*dx  + t%quadjxy(2)*dy  + t%quadjxy(3)*dz
      v_dot_x_quad_offdiag(2)  = t%quadjyz(1)*dx  + t%quadjyz(2)*dy  + t%quadjyz(3)*dz
      v_dot_x_quad_offdiag(3)  = t%quadjzx(1)*dx  + t%quadjzx(2)*dy  + t%quadjzx(3)*dz


      Jirr_dipI       = two*rd5 - five*r2*rd7
      Jirr_quadI(1)   = two*rd5 - twenty*dx2*rd7 - five*r2*rd7 + thirtyfive*dx2*r2*rd9
      Jirr_quadI(2)   = two*rd5 - twenty*dy2*rd7 - five*r2*rd7 + thirtyfive*dy2*r2*rd9
      Jirr_quadI(3)   = two*rd5 - twenty*dz2*rd7 - five*r2*rd7 + thirtyfive*dz2*r2*rd9
      Jirr_quad_offdiagI = thirtyfive*r2*rd9 - twenty*rd7


      Jirr_dipxII(1)   = v_dot_x_dip(1)*rd5 + t%dipjx(1)*dx*rd5 - five*dx2*v_dot_x_dip(1)*rd7
      Jirr_mixedII     = t%dipjx(1)*rd5 - five*dx*v_dot_x_dip(1)*rd7
      Jirr_dipxII(2)   = dy*Jirr_mixedII
      Jirr_dipxII(3)   = dz*Jirr_mixedII


      Jirr_mixedII     = t%dipjy(2)*rd5 - five*dy*v_dot_x_dip(2)*rd7
      Jirr_dipyII(1)   = dx*Jirr_mixedII
      Jirr_dipyII(2)   = v_dot_x_dip(2)*rd5 + t%dipjy(2)*dy*rd5 - five*dy2*v_dot_x_dip(2)*rd7
      Jirr_dipyII(3)   = dz*Jirr_mixedII

      Jirr_mixedII     = t%dipjz(3)*rd5 - five*dz*v_dot_x_dip(3)*rd7
      Jirr_dipzII(1)   = dx*Jirr_mixedII
      Jirr_dipzII(2)   = dy*Jirr_mixedII
      Jirr_dipzII(3)   = v_dot_x_dip(3)*rd5 + t%dipjz(3)*dz*rd5 - five*dz2*v_dot_x_dip(3)*rd7

      Jirr_quadxII(1)  = two*t%quadjx(1)*rd5 - fifteen*dx*v_dot_x_quad(1)*rd7 + thirtyfive*dx3*v_dot_x_quad(1)*rd9 - ten*t%quadjx(1)*dx2*rd7
      Jirr_mixedII     = thirtyfive*dx2*v_dot_x_quad(1)*rd9 - five*v_dot_x_quad(1)*rd7 - ten*t%quadjx(1)*dx*rd7
      Jirr_quadxII(2)  = dy*Jirr_mixedII
      Jirr_quadxII(3)  = dz*Jirr_mixedII

      Jirr_mixedII     = thirtyfive*dy2*v_dot_x_quad(2)*rd9 - five*v_dot_x_quad(2)*rd7 - ten*t%quadjy(2)*dy*rd7
      Jirr_quadyII(1)  = dx*Jirr_mixedII
      Jirr_quadyII(2)  = two*t%quadjy(2)*rd5 - fifteen*dy*v_dot_x_quad(2)*rd7 + thirtyfive*dy3*v_dot_x_quad(2)*rd9 - ten*t%quadjy(2)*dy2*rd7
      Jirr_quadyII(3)  = dz*Jirr_mixedII

      Jirr_mixedII     = thirtyfive*dz2*v_dot_x_quad(3)*rd9 - five*v_dot_x_quad(3)*rd7 - ten*t%quadjz(3)*dz*rd7
      Jirr_quadzII(1)  = dx*Jirr_mixedII
      Jirr_quadzII(2)  = dy*Jirr_mixedII
      Jirr_quadzII(3)  = two*t%quadjz(3)*rd5 - fifteen*dz*v_dot_x_quad(3)*rd7 + thirtyfive*dz3*v_dot_x_quad(3)*rd9 - ten*t%quadjz(3)*dz2*rd7

      Jirr_quadxyII(1) = t%quadjxy(2)*rd5 - five*dy*v_dot_x_quad_offdiag(1)*rd7 - five*t%quadjxy(2)*dx2*rd7 + thirtyfive*dx2*dy*v_dot_x_quad_offdiag(1)*rd9 - five*t%quadjxy(1)*dx*dy*rd7
      Jirr_quadxyII(2) = t%quadjxy(1)*rd5 - five*dx*v_dot_x_quad_offdiag(1)*rd7 - five*t%quadjxy(1)*dy2*rd7 + thirtyfive*dx*dy2*v_dot_x_quad_offdiag(1)*rd9 - five*t%quadjxy(2)*dx*dy*rd7
      Jirr_quadxyII(3) = thirtyfive*dx*dy*dz*v_dot_x_quad_offdiag(1)*rd9 - five*t%quadjxy(2)*dx*dz*rd7 - five*t%quadjxy(1)*dy*dz*rd7

      Jirr_quadyzII(1) = thirtyfive*dx*dy*dz*v_dot_x_quad_offdiag(2)*rd9 - five*t%quadjyz(3)*dx*dy*rd7 - five*t%quadjyz(2)*dx*dz*rd7
      Jirr_quadyzII(2) = t%quadjyz(3)*rd5 - five*dz*v_dot_x_quad_offdiag(2)*rd7 - five*t%quadjyz(3)*dy2*rd7 + thirtyfive*dy2*dz*v_dot_x_quad_offdiag(2)*rd9 - five*t%quadjyz(2)*dy*dz*rd7
      Jirr_quadyzII(3) = t%quadjyz(2)*rd5 - five*dy*v_dot_x_quad_offdiag(2)*rd7 - five*t%quadjyz(2)*dz2*rd7 + thirtyfive*dy*dz2*v_dot_x_quad_offdiag(2)*rd9 - five*t%quadjyz(3)*dy*dz*rd7

      Jirr_quadzxII(1) = t%quadjzx(3)*rd5 - five*dz*v_dot_x_quad_offdiag(3)*rd7 - five*t%quadjzx(3)*dx2*rd7 + thirtyfive*dx2*dz*v_dot_x_quad_offdiag(3)*rd9 - five*t%quadjzx(1)*dx*dz*rd7
      Jirr_quadzxII(2) = thirtyfive*dx*dy*dz*v_dot_x_quad_offdiag(3)*rd9 - five*t%quadjzx(3)*dx*dy*rd7 - five*t%quadjzx(1)*dy*dz*rd7
      Jirr_quadzxII(3) = t%quadjzx(1)*rd5 - five*dx*v_dot_x_quad_offdiag(3)*rd7 - five*t%quadjzx(1)*dz2*rd7 + thirtyfive*dx*dz2*v_dot_x_quad_offdiag(3)*rd9 - five*t%quadjzx(3)*dx*dz*rd7


      JirrI(1:3)     =  rd5*( t%monoj(1:3)*r2 )                                                                                        & ! monopole JI
                       + Jirr_dipI*( t%dipjx(1:3)*dx + t%dipjy(1:3)*dy + t%dipjz(1:3)*dz  )                                             & ! dipole JI
                       + t%quadjx(1:3)*Jirr_quadI(1) + t%quadjy(1:3)*Jirr_quadI(2) + t%quadjz(1:3)*Jirr_quadI(3)                       & ! diagonal quadrupole JI
                       + Jirr_quad_offdiagI*( t%quadjxy(1:3)*dx*dy + t%quadjyz(1:3)*dy*dz + t%quadjzx(1:3)*dz*dx )                      ! off-diagonal quadrupole JI

      JirrII(1:3)    =  rd5*v_dot_x_mono*d(1:3)                                                                                        & !  monopole JII
                       + Jirr_dipxII(1:3) + Jirr_dipyII(1:3) + Jirr_dipzII(1:3)                                                        & !  dipole JII
                       + Jirr_quadxII(1:3) + Jirr_quadyII(1:3) + Jirr_quadzII(1:3)                                                     & !  diagonal quadrupole JII
                       + Jirr_quadxyII(1:3) + Jirr_quadyzII(1:3) + Jirr_quadzxII(1:3)                                                    !  off-diagonal quadrupole JII


      Jirrxyz(1:3)  = Jxyz(1:3)/three + quarter*oneoverpi*JirrI(1:3)  - three*quarter*oneoverpi*JirrII(1:3)



!      rho                 = quarter/pi*three*eps2*rd5*t%charge                                                                         & ! monopole
!                           + t%dip(1)*J_dip(1) + t%dip(2)*J_dip(2) + t%dip(3)*J_dip(3)          !                                      & ! dipole
!                           + t%quad(1)*J_quad(1) + t%quad(2)*J_quad(2) + t%quad(3)*J_quad(3)                                          & ! diagonal quadrupole
!                           + t%xyquad*J_quad_offdiag(1) + t%yzquad*J_quad_offdiag(2) + t%zxquad*J_quad_offdiag(3)                       ! off-diagonal quadrupole

      B_monopole(1:3)       = (/ t%monoj(2)*dz - t%monoj(3)*dy , t%monoj(3)*dx - t%monoj(1)*dz , t%monoj(1)*dy - t%monoj(2)*dx /)

      B_dip_x(1:3)          = (/ t%dipjx(2)*dz - t%dipjx(3)*dy , t%dipjx(1)*dz - t%dipjx(3)*dx , t%dipjx(1)*dy - t%dipjx(2)*dx /)
      B_dip_y(1:3)          = (/ t%dipjy(2)*dz - t%dipjy(3)*dy , t%dipjy(1)*dz - t%dipjy(3)*dx , t%dipjy(1)*dy - t%dipjy(2)*dx /)
      B_dip_z(1:3)          = (/ t%dipjz(2)*dz - t%dipjz(3)*dy , t%dipjz(1)*dz - t%dipjz(3)*dx , t%dipjz(1)*dy - t%dipjz(2)*dx /)

      B_quad_x(1:3)         = (/ t%quadjx(2)*dz - t%quadjx(3)*dy , t%quadjx(1)*dz - t%quadjx(3)*dx , t%quadjx(1)*dy - t%quadjx(2)*dx /)
      B_quad_y(1:3)         = (/ t%quadjy(2)*dz - t%quadjy(3)*dy , t%quadjy(1)*dz - t%quadjy(3)*dx , t%quadjy(1)*dy - t%quadjy(2)*dx /)
      B_quad_z(1:3)         = (/ t%quadjz(2)*dz - t%quadjz(3)*dy , t%quadjz(1)*dz - t%quadjz(3)*dx , t%quadjz(1)*dy - t%quadjz(2)*dx /)

      B_quad_xy(1:3)        = (/ t%quadjxy(2)*dz - t%quadjxy(3)*dy , t%quadjxy(1)*dz - t%quadjxy(3)*dx , t%quadjxy(1)*dy - t%quadjxy(2)*dx /)
      B_quad_yz(1:3)        = (/ t%quadjyz(2)*dz - t%quadjyz(3)*dy , t%quadjyz(1)*dz - t%quadjyz(3)*dx , t%quadjyz(1)*dy - t%quadjyz(2)*dx /)
      B_quad_zx(1:3)        = (/ t%quadjzx(2)*dz - t%quadjzx(3)*dy , t%quadjzx(1)*dz - t%quadjzx(3)*dx , t%quadjzx(1)*dy - t%quadjzx(2)*dx /)


      Bxyz(1:3)  = rd3 *B_monopole(1:3)                            !                                                ! monopole

      Bxyz(1)    = Bxyz(1) - three *rd5*B_dip_x(1)*dx                                                              &! dipole
                  -rd3 *t%dipjy(3) - three *rd5*B_dip_y(1)*dy                                                     &
                  +rd3 *t%dipjz(2) - three *rd5*B_dip_z(1)*dz                                                     &
                  +three *rd5*B_quad_x(1)*( five*dx2*rd2 - one )                                                   &! xx - quadrupole
                  +three *rd5*B_quad_y(1)*( five*dy2*rd2 - one ) + six *rd5*t%quadjy(3)*dy                        &! yy
                  +three *rd5*B_quad_z(1)*( five*dz2*rd2 - one ) - six *rd5*t%quadjz(2)*dz                        &! zz
                  +three *rd5*dx*( five*B_quad_xy(1)*dy*rd2 + t%quadjxy(3) )                                       &! xy
                  +three *rd5*( five*B_quad_yz(1)*dy*dz*rd2 + t%quadjyz(3)*dz - t%quadjyz(2)*dy )                  &! yz
                  +three *rd5*dx*( five*B_quad_zx(1)*dz*rd2 - t%quadjzx(2) )                                        ! zx


      Bxyz(2)    = Bxyz(2) + t%dipjx(3)*rd3  + three *rd5*B_dip_x(2)*dx                                           &! dipole
                   + three *rd5*B_dip_y(2)*dy                                                                      &
                  -rd3 *t%dipjz(1) + three *rd5*B_dip_z(2)*dz                                                     &
                  -three *rd5*B_quad_x(2)*( five*dx2*rd2 - one ) - six *rd5*t%quadjx(3)*dx                        &! xx - quadrupole
                  -three *rd5*B_quad_y(2)*( five*dy2*rd2 - one )                                                   &! yy
                  -three *rd5*B_quad_z(2)*( five*dz2*rd2 - one ) + six *rd5*t%quadjz(1)*dz                        &! zz
                  -three *rd5*dy*( five*B_quad_xy(2)*dx*rd2 + t%quadjxy(3) )                                       &! xy
                  -three *rd5*dy*( five*B_quad_yz(2)*dz*rd2 - t%quadjyz(1) )                                       &! yz
                  -three *rd5*( five*B_quad_zx(2)*dx*dz*rd2 + t%quadjzx(3)*dz - t%quadjzx(1)*dx )                   ! xz


      Bxyz(3)    = Bxyz(3) - rd3 *t%dipjx(2) - three *rd5*B_dip_x(3)*dx                                           &! dipole
                  +rd3 *t%dipjy(1) - three *rd5*B_dip_y(3)*dy                                                     &
                  - three *rd5*B_dip_z(3)*dz                                                                       &
                  +three *rd5*B_quad_x(3)*( five*dx2*rd2 - one ) + six *rd5*t%quadjx(2)*dx                        &! xx - quadrupole
                  +three *rd5*B_quad_y(3)*( five*dy2*rd2 - one ) - six *rd5*t%quadjy(1)*dy                        &! yy
                  +three *rd5*B_quad_z(3)*( five*dz2*rd2 - one )                                                   &! zz
                  +three *rd5*( five*B_quad_xy(3)*dy*dx*rd2 + t%quadjxy(2)*dy - t%quadjxy(1)*dx  )                 &! xy
                  +three *rd5*( five*B_quad_yz(3)*dy*dz*rd2 - t%quadjyz(1)*dz )                                    &! yz
                  +three *rd5*dz*( five*B_quad_zx(3)*dx*rd2 + t%quadjzx(2) )                                        ! zx




      sinh_1                   = asinh( sqrt( r2/eps2 ) )

      AI_dip                   = rd3 + half*over2*rd - rEps*over4 - half*eps2*rd*over4 + half*three*eps2*over5*sinh_1

      AI_quad(1)               = rd3 - three*dx2*rd5 + half*over2*rd - rEps*over4 + half*three*eps2*over5*sinh_1           &
                                - two*dx2*over4*rd - half*dx2*over2*rd3 + four*dx2*rEps*over6 - half*eps2*rd*over4 + half*eps2*dx2*rd3*over4 &
                                - half*fifteen*eps2*dx2*over7*sinh_1 + half*seven*eps2*dx2*rd*over6

      AI_quad(2)               = rd3 - three*dy2*rd5 + half*over2*rd - rEps*over4 + half*three*eps2*over5*sinh_1           &
                                - two*dy2*over4*rd - half*dy2*over2*rd3 + four*dy2*rEps*over6 - half*eps2*rd*over4 + half*eps2*dy2*rd3*over4 &
                                - half*fifteen*eps2*dy2*over7*sinh_1 + half*seven*eps2*dy2*rd*over6

      AI_quad(3)               = rd3 - three*dz2*rd5 + half*over2*rd - rEps*over4 + half*three*eps2*over5*sinh_1           &
                                - two*dz2*over4*rd - half*dz2*over2*rd3 + four*dz2*rEps*over6 - half*eps2*rd*over4 + half*eps2*dz2*rd3*over4 &
                                - half*fifteen*eps2*dz2*over7*sinh_1 + half*seven*eps2*dz2*rd*over6

      AI_quad_offdiag          =  three*rd5 + two*over4*rd  + half*over2*rd3 - four*rEps*over6                       &
                                - half*eps2*over4*rd3 + half*fifteen*eps2*sinh_1*over7 - seven*half*eps2*rd*over6


      AII1                     =  three*over*sinh_1 - (r2/eps2 + three)*rd
      AII2                     =  (r2/eps2 + three)*rd3  - two*rd/eps2 -three*over3*sinh_1 + three*over2*rd
      AII3                     =  two*rd/eps2 + three*over3*sinh_1 - four*dx2/eps2*rd3 - ( two*eps2 + d2 )*rd3/eps2 - nine*dx2*sinh_1*over5 &
                                  - three*over2*rd + dx2*( six*eps2 + three*d2 )*rd5/eps2 + nine*dx2*rd*over4 + three*dx2*rd3*over2

      AII4                     =  two*rd/eps2 + three*over3*sinh_1 - (r2/eps2 + three)*rd3 - four*dx2*rd3/eps2 &
                                  - three*rd*over2 - nine*dx2*sinh_1*over5 + three*dx2*(r2/eps2 + three)*rd5   &
                                  + nine*dx2*rd*over4 + three*dx2*rd3*over2

      AII5                     =  two*rd/eps2 + three*over3*sinh_1 - (r2/eps2 + three)*rd3 - four*dy2*rd3/eps2 &
                                  - three*rd*over2 - nine*dy2*sinh_1*over5 + three*dy2*(r2/eps2 + three)*rd5   &
                                  + nine*dy2*rd*over4 + three*dy2*rd3*over2

      AII6                     =  two*rd/eps2 + three*over3*sinh_1 - (r2/eps2 + three)*rd3 - four*dz2*rd3/eps2 &
                                  - three*rd*over2 - nine*dz2*sinh_1*over5 + three*dz2*(r2/eps2 + three)*rd5   &
                                  + nine*dz2*rd*over4 + three*dz2*rd3*over2

      AII7                     =  three*(r2/eps2 + three)*rd5 - nine*over5*sinh_1 - four*rd3/eps2         &
                                  + nine*rd*over4 + three*rd3*over2

      AII8                     =  r**3 + three*eps2*r - three*eps2*sinh_1*rEps




      AII_mono                 = -over4*AII1*v_dot_x_mono

      AII_dipx(1)              =  over4*( four*dx2*AII1*v_dot_x_dip(1)*over2 -  t%dipjx(1)*dx*AII1 - dx2*AII2*v_dot_x_dip(1) - AII1*v_dot_x_dip(1) )

      AII_dipole_mix_x           = ( four*dx*AII1*v_dot_x_dip(1)*over2 - dx*AII2*v_dot_x_dip(1) - AII1*t%dipjx(1) )
      AII_dipole_mix_y           = ( four*dy*AII1*v_dot_x_dip(2)*over2 - dy*AII2*v_dot_x_dip(2) - AII1*t%dipjy(2) )
      AII_dipole_mix_z           = ( four*dz*AII1*v_dot_x_dip(3)*over2 - dz*AII2*v_dot_x_dip(3) - AII1*t%dipjz(3) )

      AII_dipy(1)              =  over4*dx*AII_dipole_mix_y

      AII_dipz(1)              =  over4*dx*AII_dipole_mix_z



      AII_dipx(2)              = over4*dy*AII_dipole_mix_x

      AII_dipy(2)              = over4*( four*dy2*AII1*v_dot_x_dip(2)*over2 -  t%dipjy(2)*dy*AII1 - dy2*AII2*v_dot_x_dip(2) - AII1*v_dot_x_dip(2) )

      AII_dipz(2)              = over4*dy*AII_dipole_mix_z




      AII_dipx(3)              = over4*dz*AII_dipole_mix_x

      AII_dipy(3)              = over4*dz*AII_dipole_mix_y

      AII_dipz(3)              = over4*( four*dz2*AII1*v_dot_x_dip(3)*over2 -  t%dipjz(3)*dz*AII1 - dz2*AII2*v_dot_x_dip(3) - AII1*v_dot_x_dip(3) )

      AII_quadx(1)             = over4*( eight*t%quadjx(1)*dx2*AII1*over2 - two*dx*v_dot_x_quad(1)*AII2 - two*t%quadjx(1)*AII1 + dx*v_dot_x_quad(1)*AII4   &
                                 + eight*dx3*v_dot_x_quad(1)*AII2*over2 + twelve*dx*AII1*v_dot_x_quad(1)*over2 - twentyfour*dx3*AII1*v_dot_x_quad(1)*over4 - two*t%quadjx(1)*dx2*AII2  )

      AII_quadrupole_mix       = v_dot_x_quad(1)*AII4 + four*AII1*v_dot_x_quad(1)*over2 - two*t%quadjx(1)*dx*AII2 + eight*t%quadjx(1)*dx*AII1*over2 &
                                 + eight*dx2*v_dot_x_quad(1)*AII2*over2 - twentyfour*dx2*AII1*over4*v_dot_x_quad(1)


      AII_quadx(2)             = over4*dy*AII_quadrupole_mix

      AII_quadx(3)             = over4*dz*AII_quadrupole_mix


      AII_quadrupole_mix       = v_dot_x_quad(2)*AII5 + four*AII1*v_dot_x_quad(2)*over2 - two*t%quadjy(2)*dy*AII2 + eight*t%quadjy(2)*dy*AII1*over2 &
                                 + eight*dy2*v_dot_x_quad(2)*AII2*over2 - twentyfour*dy2*AII1*over4*v_dot_x_quad(2)


      AII_quady(1)             = over4*dx*AII_quadrupole_mix

      AII_quady(2)             = over4*( eight*t%quadjy(2)*dy2*AII1*over2 - two*dy*v_dot_x_quad(2)*AII2 - two*t%quadjy(2)*AII1 + dy*v_dot_x_quad(2)*AII5 &
                                 + eight*dy3*v_dot_x_quad(2)*AII2*over2  + twelve*dy*AII1*v_dot_x_quad(2)*over2 - twentyfour*dy3*AII1*v_dot_x_quad(2)*over4 - two*t%quadjy(2)*dy2*AII2   )

      AII_quady(3)             = over4*dz*AII_quadrupole_mix


      AII_quadrupole_mix       = v_dot_x_quad(3)*AII6 - two*t%quadjz(3)*dz*AII2 + four*AII1*v_dot_x_quad(3)*over2 + eight*t%quadjz(3)*dz*AII1*over2 &
                                 + eight*dz2*v_dot_x_quad(3)*AII2*over2 - twentyfour*dz2*AII1*v_dot_x_quad(3)*over4

      AII_quadz(1)             = over4*dx*AII_quadrupole_mix

      AII_quadz(2)             = over4*dy*AII_quadrupole_mix

      AII_quadz(3)             = over4*( eight*t%quadjz(3)*dz2*AII1*over2 - two*dz*v_dot_x_quad(3)*AII2 - two*t%quadjz(3)*AII1 + dz*v_dot_x_quad(3)*AII6 &
                                 + eight*dz3*v_dot_x_quad(3)*AII2*over2  + twelve*dz*AII1*v_dot_x_quad(3)*over2 - twentyfour*dz3*AII1*v_dot_x_quad(3)*over4 - two*t%quadjz(3)*dz2*AII2   )



      AII_quadxy(1)            = over4*( four*t%quadjxy(2)*dx2*AII1*over2 - v_dot_x_quad_offdiag(1)*dy*AII2 - t%quadjxy(1)*dx*dy*AII2 - t%quadjxy(2)*AII1  &
                                 + four*dx2*dy*v_dot_x_quad_offdiag(1)*AII2*over2 + four*dy*AII1*v_dot_x_quad_offdiag(1)*over2 + v_dot_x_quad_offdiag(1)*dx2*dy*AII7 &
                                 - t%quadjxy(2)*dx2*AII2 + four*t%quadjxy(1)*dx*dy*AII1*over2 + four*dx2*dy*v_dot_x_quad_offdiag(1)*AII2*over2 - twentyfour*dx2*dy*AII1*v_dot_x_quad_offdiag(1)*over4  )

      AII_quadxy(2)            = over4*( four*t%quadjxy(1)*dy2*AII1*over2 - v_dot_x_quad_offdiag(1)*dx*AII2 - t%quadjxy(1)*dy2*AII2 - t%quadjxy(1)*AII1  &
                                 + four*dx*dy2*v_dot_x_quad_offdiag(1)*AII2*over2 + four*dx*AII1*v_dot_x_quad_offdiag(1)*over2 + v_dot_x_quad_offdiag(1)*dx*dy2*AII7 &
                                 - t%quadjxy(2)*dx*dy*AII2 + four*t%quadjxy(2)*dx*dy*AII1*over2 + four*dx*dy2*v_dot_x_quad_offdiag(1)*AII2*over2 - twentyfour*dx*dy2*AII1*v_dot_x_quad_offdiag(1)*over4  )

      AII_quadxy(3)            = dz*over4*( v_dot_x_quad_offdiag(1)*dx*dy*AII7 + ( t%quadjxy(1)*dy + t%quadjxy(2)*dx )* ( - AII2 + four*AII1*over2)                  &
                                 + eight*dx*dy*v_dot_x_quad_offdiag(1)*AII2*over2 - twentyfour*dx*dy*AII1*v_dot_x_quad_offdiag(1)*over4   )


      AII_quadyz(1)            = dx*over4*( dy*dz*v_dot_x_quad_offdiag(2)*AII7 + ( t%quadjyz(2)*dz + t%quadjyz(3)*dy )*( - AII2 + four*AII1*over2)                   &
                                 + eight*dy*dz*v_dot_x_quad_offdiag(2)*AII2*over2 - twentyfour*dy*dz*AII1*v_dot_x_quad_offdiag(2)*over4   )


      AII_quadyz(2)            = over4*( four*t%quadjyz(3)*dy2*AII1*over2 - v_dot_x_quad_offdiag(2)*dz*AII2 - t%quadjyz(3)*dy2*AII2 - t%quadjyz(2)*dy*dz*AII2  &
                                 - t%quadjyz(3)*AII1 + four*dy2*dz*v_dot_x_quad_offdiag(2)*AII2*over2 + four*dz*AII1*v_dot_x_quad_offdiag(2)*over2 + v_dot_x_quad_offdiag(2)*dy2*dz*AII7 &
                                 + four*t%quadjyz(2)*dy*dz*AII1*over2 + four*v_dot_x_quad_offdiag(2)*dy2*dz*AII2*over2 - twentyfour*dy2*dz*v_dot_x_quad_offdiag(2)*AII1*over4  )

      AII_quadyz(3)            = over4*( four*t%quadjyz(2)*dz2*AII1*over2 - v_dot_x_quad_offdiag(2)*dy*AII2 - t%quadjyz(3)*dy*dz*AII2 - t%quadjyz(2)*dz2*AII2  &
                                 - t%quadjyz(2)*AII1 + four*dz2*dy*v_dot_x_quad_offdiag(2)*AII2*over2 + four*dy*AII1*v_dot_x_quad_offdiag(2)*over2 + v_dot_x_quad_offdiag(2)*dy*dz2*AII7 &
                                 + four*t%quadjyz(3)*dy*dz*AII1*over2 + four*v_dot_x_quad_offdiag(2)*dy*dz2*AII2*over2 - twentyfour*dy*dz2*v_dot_x_quad_offdiag(2)*AII1*over4  )


      AII_quadzx(1)            = over4*( four*t%quadjzx(3)*dx2*AII1*over2 - v_dot_x_quad_offdiag(3)*dz*AII2 - t%quadjzx(1)*dx*dz*AII2 - t%quadjzx(3)*AII1  &
                                 + four*v_dot_x_quad_offdiag(3)*dx2*dz*AII2*over2 + four*dz*v_dot_x_quad_offdiag(3)*AII1*over2  + v_dot_x_quad_offdiag(3)*dx2*dz*AII7 &
                                 - t%quadjzx(3)*dx2*AII2 + four*t%quadjzx(1)*dx*dz*AII1*over2 + four*v_dot_x_quad_offdiag(3)*dx2*dz*AII2*over2 - twentyfour*dx2*dz*v_dot_x_quad_offdiag(3)*AII1*over4  )


      AII_quadzx(2)            = dy*over4*( dx*dz*v_dot_x_quad_offdiag(3)*AII7 + ( t%quadjzx(1)*dz + t%quadjyz(3)*dx )*( - AII2 + four*AII1*over2)                   &
                                 + eight*dx*dz*v_dot_x_quad_offdiag(3)*AII2*over2 - twentyfour*dx*dz*AII1*v_dot_x_quad_offdiag(3)*over4   )

      AII_quadzx(3)            = over4*( four*t%quadjzx(1)*dz2*AII1*over2 - v_dot_x_quad_offdiag(3)*dx*AII2 - t%quadjzx(3)*dx*dz*AII2 - t%quadjzx(1)*dz2*AII2  &
                                 - t%quadjzx(1)*AII1 + four*dz2*dx*v_dot_x_quad_offdiag(3)*AII2*over2 + four*dx*AII1*v_dot_x_quad_offdiag(3)*over2 + v_dot_x_quad_offdiag(3)*dz2*dx*AII7 &
                                 + four*t%quadjzx(3)*dx*dz*AII1*over2 + four*v_dot_x_quad_offdiag(3)*dx*dz2*AII2*over2 - twentyfour*dz2*dx*v_dot_x_quad_offdiag(3)*AII1*over4  )



      AII(1:3)  = AII_mono*d(1:3) + AII_dipx(1:3) + AII_dipy(1:3) + AII_dipz(1:3) &
                  + AII_quadx(1:3) + AII_quady(1:3) + AII_quadz(1:3)              &
                  + AII_quadxy(1:3) + AII_quadyz(1:3) + AII_quadzx(1:3)

      Axyz(1:3)  = t%monoj(1:3)*( rd - half*rEps*over2 + half*eps2*over3*sinh_1 )                                                      & ! monopole
                 - AI_dip*( t%dipjx(1:3)*dx + t%dipjy(1:3)*dy + t%dipjz(1:3)*dz )                                                      & ! dipole I
                 - ( t%quadjx(1:3)*AI_quad(1) + t%quadjy(1:3)*AI_quad(2) + t%quadjz(1:3)*AI_quad(3)   )                                & ! diagonal quadrupole I
                 + AI_quad_offdiag*( t%quadjxy(1:3)*dx*dy + t%quadjyz(1:3)*dz*dy + t%quadjzx(1:3)*dz*dx   )                              ! off-diagonal quadrupole I


      Axyz(1:3)  = Axyz(1:3)  + half*eps2 *AII(1:3)


      phi            = half*phi
      exyz           = half*exyz
      Jxyz           = half*Jxyz
      Bxyz           = half*Bxyz
      Jirrxyz        = half*Jirrxyz
      Axyz           = half*Axyz
      
      Axyz           =  Axyz/vtilde  
      Bxyz           =  Bxyz/vtilde

    end subroutine calc_force_darwin_3D


    !>
    !> CALC_FORCE_LJ
    !>
    !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
    !> shifted by the lattice vector vbox
    !> results are returned exyz, phi
    !>
    subroutine calc_force_LJ(t, d, r2, aii2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), r2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi
      real(kind_physics), intent(in) :: aii2
      real(kind_physics) :: flj, epsc2, aii2_r2,aii4_r4, r, fljrd

      ! epsc should be > a_ii to get evenly spaced ions
      epsc2 = 0.8_kind_physics * aii2

      ! Force is repulsive up to and just beyond aii
      if (r2 > epsc2) then
          aii2_r2 = aii2/r2
      else
          aii2_r2 = aii2/epsc2
      endif

      aii4_r4 = aii2_r2*aii2_r2

      flj = two*(aii4_r4*aii4_r4) - aii4_r4

      ! potential
      phi = zero

      !  forces
      r     = sqrt(r2)
      fljrd = flj/r

      exyz = d*fljrd
    end subroutine calc_force_LJ


    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_direct(t, d, dist2, exyz, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exyz(3), phi

      real(kind_physics) :: rd,r,rd3charge

      r         = sqrt(dist2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd        = one/r
      rd3charge = t%charge*rd*rd*rd

      phi  = t%charge*rd
      exyz = rd3charge*d
    end subroutine calc_force_coulomb_3D_direct


    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox    !>

    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D_direct(t, d, d2, exy, phi)
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in) :: d(2), d2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(2), phi

      real(kind_physics) :: rd2charge

      phi       = - half*t%charge*log(d2)
      rd2charge = t%charge/d2
      exy       = rd2charge*d

    end subroutine calc_force_coulomb_2D_direct



    subroutine calc_force_darwin_2D_direct(t, d, d2,eps2, exy, Axy, Jxy, Jirrxy, Bz, phi)
      use module_globals, only: vtilde  
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) :: exy(1:2),Axy(1:2),Jxy(1:2), Jirrxy(1:2),Bz,phi

      real(kind_physics)              :: vdotxM,over4,dx,dy,rd2,rd4,r2,logR2e,over2,A1,A2,dx2,dy2

      dx = d(1)
      dy = d(2)
      dx2= dx**2
      dy2= dy**2

      r2 =  dx**2 + dy**2
      over2 = one/r2
      over4 = over2**2
      rd2 = one/d2
      rd4 = rd2*rd2

      phi = - t%charge * log(d2)

      exy = two*t%charge * d* rd2


!      Jxy = eps2/pi*rd4*t%monoj(1:2)
      Jxy = eps2/pi*rd4*t%charge

      Bz  = two*rd2 *( t%monoj(1)*dy - t%monoj(2)*dx )

      vdotxM   = t%monoj(1)*dx + t%monoj(2)*dy

      Jirrxy(1:2) = two*( t%monoj(1:2)*rd2 - two*d(1:2)*rd4*vdotxM )


      logR2e = log( one + r2/eps2 )
      A1     = eps2*over2*logR2e - one
      A2     = rd2 - over2*logR2e

      Axy(1:2) = half *t%monoj(1:2)*( eps2*over2*logR2e  - log(d2) )                                             &
              - A1*vdotxM*d(1:2)*over2
      
      phi           = half*phi
      exy           = half*exy
      Jxy           = half*Jxy
      Bz            = half*Bz/vtilde
      Jirrxy        = half*Jirrxy
      Axy           = half*Axy/vtilde
      

    end subroutine calc_force_darwin_2D_direct



    subroutine calc_force_darwin_2D3V_direct(t, d, d2,eps2, exy, Axy, Jxy, Jirrxy, Bxy, phi, gradxA, gradyA)
      use module_globals, only: vtilde 
      use module_tool   , only: cross_product,double_cross_product_left
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(2), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) ::  exy(1:2),Axy(1:3),Jxy(1:3),Jirrxy(1:3),Bxy(1:3),phi,gradxA(1:3),gradyA(1:3)

      real(kind_physics) :: dx,dy,rd2,rd4,rd6,rd8,dx2,dy2,dx3,dy3,r2,logR2e,over2,over4,over6,over8, &
                            phi_mono,rho_mono, E_mono(1:3),Ex_dip(1:3),Ey_dip(1:3),logTmp,A1_mono,   &
                            A1_dip(1:2),A1_quad(1:3),A2_mono,A3_mono(1:3),A3x_dip(1:3),A3y_dip(1:3), &   
                            dist(1:3),Jir_mono,over,phi_dip(1:2),phi_quad(1:3)

      dx = d(1)
      dy = d(2)
      dist(1:2) = d(1:2)
      dist(3)   = zero
      r2 =  dx**2 + dy**2
      over  = sqrt(one/r2)
      over2 = one/r2
      over4 = over2*over2
      over6 = over4*over2
      over8 = over6*over2
      rd2 = one/d2
      rd4 = rd2*rd2
      rd6 = rd4*rd2
      rd8 = rd6*rd2

      dx2 = dx *dx
      dy2 = dy *dy
      dx3 = dx2*dx
      dy3 = dy2*dy
      
      logTmp     = log(r2/eps2 + one)
      logR2e     = eps2*over2*logTmp 
      
      !!! Rho/J must be multiply by eps2/pi 
      rho_mono   =  rd4
      
      !!! Jirr must be multiply by 2 
      Jir_mono      = rd2
      
      phi_mono   = -log(d2)
      phi_dip    = -two*d(1:2)*rd2
      phi_quad(1)=  four*dx2*rd4 - two*rd2
      phi_quad(2)=  four*dy2*rd4 - two*rd2
      phi_quad(3)=  four*dy*dx*rd4 

      E_mono(1:2)= -phi_dip
      E_mono(3)  =  zero
      A1_mono    =  logR2e - one
      A1_dip(1:2)=  two*eps2*d(1:2)*rd2*over2 - two*over2*logR2e*d(1:2)
      
      
      A3_mono(1:3)  = over*dist(1:3)
      A3x_dip(1)    = over - dx2*over*over2
      A3x_dip(2)    =      - dx*dy*over*over2
      A3x_dip(3)    = zero
      A3y_dip(1)    =      - dx*dy*over*over2 
      A3y_dip(2)    = over - dy2*over*over2
      A3y_dip(3)    = zero
      
      phi        = ( t%charge*phi_mono      )                                                              ! Quadrupole
      phi        = half*phi
                   
      exy(1:2)   = ( t%charge*E_mono(1:2)   )            
      exy(1:2)   = half*exy(1:2)
      
      Jxy(1:3)   = ( t%monoj(1:3)*rho_mono )                
      Jxy(1:3)   = half*eps2/pi*Jxy(1:3)
                  
                  
      Jirrxy(1:3)= ( t%monoj(1:3)*(three*half*Jir_mono - eps2*rho_mono)                                     &
                   + double_cross_product_left( E_mono(1:3), E_mono(1:3), t%monoj(1:3)  ) )!                 
                 
      Jirrxy(1:3)= half*oneoverpi*Jirrxy(1:3)
                 
      Bxy(1:3)   = cross_product(t%monoj(1:3), E_mono(1:3) )   
                  
      Bxy(1:3)   = half/vtilde*Bxy(1:3)  
      
      Axy(1:3)   = half*t%monoj(1:3)*( two - logR2e - log(d2)  )                                            &
                 + A1_mono*double_cross_product_left( A3_mono(1:3), t%monoj(1:3), A3_mono(1:3)  )  
!                 
      Axy(3)     = t%monoj(3)*phi_mono                                                              ! Quadrupole 
                  
      Axy(1:3)   = half/vtilde*Axy(1:3)
      
      
                 
      gradxA(1:3)=  half*( t%monoj(1:3)*( -A1_dip(1) + phi_dip(1) ) )                                       &! Monopole 
                 + A1_mono*( double_cross_product_left( A3x_dip(1:3), t%monoj(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%monoj(1:3),A3x_dip(1:3)  ) )        &
                 + A1_dip(1)*( double_cross_product_left(A3_mono(1:3), t%monoj(1:3),A3_mono(1:3)  ) )       &
                 + half*( t%dipjx(1:3)*( -A1_quad(1) + phi_quad(1) ) )                                      
                 
      gradxA(3)  = t%monoj(3)*phi_dip(1)                                                                     
                  
      gradyA(1:3)= half*( t%monoj(1:3)*( -A1_dip(2) + phi_dip(2) ) )                                        &! Monopole 
                 + A1_mono*( double_cross_product_left( A3y_dip(1:3), t%monoj(1:3), A3_mono(1:3)  ) )       &
                 + A1_mono*( double_cross_product_left( A3_mono(1:3), t%monoj(1:3),A3y_dip(1:3)  ) )        &
                 + A1_dip(2)*( double_cross_product_left(A3_mono(1:3), t%monoj(1:3),A3_mono(1:3)  ) )                  
      
      gradyA(3)  = t%monoj(3)*phi_dip(2)                                                                     
      
      gradxA(1:3)= half/vtilde*gradxA(1:3)
      gradyA(1:3)= half/vtilde*gradyA(1:3)
      
      

    end subroutine calc_force_darwin_2D3V_direct

    subroutine calc_force_darwin_3D_direct(t, d, d2,eps2, exyz, Axyz, Jxyz, Jirrxyz, Bxyz, phi)
      use module_globals, only: vtilde  
      implicit none

      type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
      real(kind_physics), intent(in)  :: d(1:3), d2, eps2 !< separation vector and magnitude**2 precomputed in walk_single_particle
      real(kind_physics), intent(out) :: exyz(1:3),Axyz(1:3),Jxyz(1:3),Jirrxyz(1:3),Bxyz(1:3),phi
      real(kind_physics)              :: vcrossxM(1:3),v_dot_x_mono,dx,dy,dz,rEps,rd,rd3
      real(kind_physics)              :: rd5,rd3charge,over,over2,over3,over4,r2


      dx            = d(1)
      dy            = d(2)
      dz            = d(3)
      r2            = dx**2 + dy**2 + dz**2
      rEps          = sqrt(d2) ! eps2 is added in calling routine to have plummer intead of coulomb here
      rd            = one/rEps
      rd3           = rd*rd*rd
      rd5           = rd3*rd*rd

      over          = one/sqrt( r2 )
      over2         = over**2
      over3         = over2*over
      over4         = over3*over


      rd3charge     = t%charge*rd3
      v_dot_x_mono  = t%monoj(1)*dx    + t%monoj(2)*dy    + t%monoj(3)*dz

      vcrossxM(1)   = t%monoj(2)*dz - t%monoj(3)*dy
      vcrossxM(2)   = t%monoj(3)*dx - t%monoj(1)*dz
      vcrossxM(3)   = t%monoj(1)*dy - t%monoj(2)*dx

      phi           = t%charge*rd
      exyz          = rd3charge*d

      Jxyz(1:3)     = three*(quarter*oneoverpi)*eps2*t%monoj(1:3)*rd5
      Jirrxyz(1:3)  = Jxyz(1:3)/three  + quarter*oneoverpi*rd5*t%monoj(1:3)*r2 - three*quarter*oneoverpi*rd5*v_dot_x_mono*d(1:3)
!      rho           = three*(quarter/pi)*eps2*t%charge*rd5
      Bxyz(1:3)     = one *rd3*vcrossxM(1:3)
      Axyz(1:3)     = one *t%monoj(1:3)*( rd - half*rEps*over2 + half*eps2*over3*asinh( sqrt( r2/eps2 ) )  ) &
                    - half *over4*eps2*d(1:3)*( three*over*asinh( sqrt( r2/eps2 ) ) - ( r2/eps2 + three )*rd )*v_dot_x_mono

      phi            = half*phi
      exyz           = half*exyz
      Jxyz           = half*Jxyz
      Bxyz           = half*Bxyz/vtilde
      Jirrxyz        = half*Jirrxyz
      Axyz           = half*Axyz/vtilde

      
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
      real(kind_physics) :: rd,r,rd3
      real(kind_physics), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kind_physics) :: ome, rol, lambda, q, fprefac

      q = t%charge

      ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
      if (particle%data%q * q < 0.) then
        ! e-i or i-e interaction
        lambda = 1.00027227_8 * kelbg_invsqrttemp
      else
        if ( q > 0. ) then
          ! i-i interaction
          lambda = 0.03300355_8 * kelbg_invsqrttemp
        else
          ! e-e interaction
          lambda = 1.41421356_8 * kelbg_invsqrttemp
        endif
      endif

      r   = sqrt(dist2)
      rd  = one / r
      rd3 = rd*rd*rd
      rol = r  / lambda        !< "r over lambda"
      ome = 1  - exp(-rol*rol) !< "one minus exp(stuff)"

      ! potential
      phi = q * rd  * (ome + sqrtpi*rol*(1-erf(rol)))
      !  forces
      fprefac = q * rd3 * ome
      exyz    = fprefac * d
    end subroutine calc_force_kelbg_3D_direct
end module



