! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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
  use module_pepc_types
  use module_interaction_specific_types
  #ifndef NO_SPATIAL_INTERACTION_CUTOFF
  use module_mirror_boxes, only: spatial_interaction_cutoff
  #endif
  implicit none
  save
  private

    integer, parameter :: kfp       =  8 ! numeric precision (kind value)
    ! shortcut notations
    real(kfp), parameter :: zero    =  0._kfp
    real(kfp), parameter :: one     =  1._kfp
    real(kfp), parameter :: two     =  2._kfp
    real(kfp), parameter :: three   =  3._kfp
    real(kfp), parameter :: four    =  4._kfp
    real(kfp), parameter :: five    =  5._kfp
    real(kfp), parameter :: eight   =  8._kfp
    real(kfp), parameter :: nine    =  9._kfp
    real(kfp), parameter :: half    =  0.5_kfp

    public calc_force_coulomb_3D
    public calc_force_coulomb_3D_direct
    public calc_force_coulomb_2D
    public calc_force_coulomb_2D_direct
    public calc_force_LJ
    public calc_force_kelbg_3D_direct


  contains

    #if defined(__TOS_BGQ__) && ( defined(__IBMC__) || defined(__IBMCPP__) )
      !> helper macros for QPX, see
      !> http://pic.dhe.ibm.com/infocenter/compbg/v121v141/topic/com.ibm.xlf141.bg.doc/language_ref/vmxintrinsics.html
      !> for the basic instructions
      !> TODO: in some places, even reordering of the expressions lead to performance improvement around 2%... - someone might want to try this out
      ! a * b * c
      #define VEC_MUL3(a,b,c) VEC_MUL(a, VEC_MUL(b,c))
      ! a * b + c * d
      #define VEC_aTbMcTd(a,b,c,d) VEC_MSUB(a,b,VEC_MUL(c,d))
      ! a * b + c * d + e * f
      #define VEC_aTbPcTdPeTf(a,b,c,d,e,f) VEC_MADD(a,b,VEC_MADD(c,d,VEC_MUL(e,f)))
      ! a * b * c -d
      #define VEC_MUL3SUB(a,b,c,d) VEC_MSUB(a, VEC_MUL(b,c), d)
      ! performs variable(index) += summand
      #define VEC_ACCUMULATE(index, variable, summand) VEC_STA(VEC_ADD(summand,VEC_LDA(index,variable)),index,variable)
      ! performs variable(index) += a * bsummand
      #define VEC_MACCUMULATE(index, variable, a, b) VEC_STA(VEC_MADD(a,b,VEC_LDA(index,variable)),index,variable)
    #endif

    !>
    !> Calculates 3D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D(delta, dist2, particle_pack, t, eps2)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

  #if defined(__TOS_BGQ__) && ( defined(__IBMC__) || defined(__IBMCPP__) )
      vector(real(kfp)) :: rd,dx,dy,dz,dx2,dy2,dz2,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
      vector(real(kfp)) :: veps2, vone, vthree, vfive, vhalf, vdist2, vtcharge, vtdip(3), vtquad(3), vtxyquad, vtzxquad, vtyzquad, &
                           vfive_rd7, dx_rd5, dy_rd5, dz_rd5, vtcharge_rd3, vfive_rd7_dx, vfive_rd7_dy, vfive_rd7_dz,vfive_rd7_dxdydz
      integer(kind_particle) :: ip, np

      
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        VECTOR(REAL(kfp)) :: vcut(3), include_mask
        
        do ip=1,3
          vcut(ip) = VEC_SPLATS(spatial_interaction_cutoff(ip))
        end do
      #endif  
      
      veps2    = VEC_SPLATS(eps2)
      vone     = VEC_SPLATS(one)
      vthree   = VEC_SPLATS(three)
      vfive    = VEC_SPLATS(five)
      vhalf    = VEC_SPLATS(half)

      vtcharge = VEC_SPLATS(t%charge)
      do ip=1,3
        vtdip(ip)  = VEC_SPLATS(t%dip(ip))
        vtquad(ip) = VEC_SPLATS(t%quad(ip))
      end do
      vtxyquad = VEC_SPLATS(t%xyquad)
      vtyzquad = VEC_SPLATS(t%yzquad)
      vtzxquad = VEC_SPLATS(t%zxquad)

      np = 8*size(dist2, kind = kind_particle)
      
      !ibm* assert(nodeps)
      do ip=0,np-8,32
        vdist2 = VEC_LDA(ip, dist2)

        rd  = VEC_SWDIV_NOCHK(vone, VEC_SWSQRT_NOCHK( VEC_ADD(vdist2, veps2) ))
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
          include_mask = VEC_AND( VEC_AND( 
              VEC_CMPLT(VEC_LDA(      ip, delta), vcut(1))  , &
              VEC_CMPLT(VEC_LDA(   np+ip, delta), vcut(2)) ), &
              VEC_CMPLT(VEC_LDA(np+np+ip, delta), vcut(3)) )
          rd = VEC_SEL(vzero, rd, include_mask)
        #endif
        rd2 = VEC_MUL(rd,  rd )
        rd3 = VEC_MUL(rd,  rd2)
        rd5 = VEC_MUL(rd3, rd2)
        rd7 = VEC_MUL(rd5, rd2)

        dx = VEC_LDA(      ip, delta)
        dy = VEC_LDA(   np+ip, delta)
        dz = VEC_LDA(np+np+ip, delta)
        dx2 = VEC_MUL(dx, dx)
        dy2 = VEC_MUL(dy, dy)
        dz2 = VEC_MUL(dz, dz)

        fd1 = VEC_SUB(VEC_MUL(VEC_MUL(vthree, dx2), rd5), rd3)
        fd2 = VEC_SUB(VEC_MUL(VEC_MUL(vthree, dy2), rd5), rd3)
        fd3 = VEC_SUB(VEC_MUL(VEC_MUL(vthree, dz2), rd5), rd3)
        fd4 = VEC_MUL(VEC_MUL(VEC_MUL(vthree, dx), dy), rd5)
        fd5 = VEC_MUL(VEC_MUL(VEC_MUL(vthree, dy), dz), rd5)
        fd6 = VEC_MUL(VEC_MUL(VEC_MUL(vthree, dx), dz), rd5)
        
        vfive_rd7    = VEC_MUL(vfive, rd7)
        vfive_rd7_dx = VEC_MUL(vfive_rd7, dx)
        vfive_rd7_dy = VEC_MUL(vfive_rd7, dy)
        vfive_rd7_dz = VEC_MUL(vfive_rd7, dz)
        vfive_rd7_dxdydz = VEC_MUL3(vfive_rd7_dx, dy, dz)
        dx_rd5       = VEC_MUL(dx, rd5)
        dy_rd5       = VEC_MUL(dy, rd5)
        dz_rd5       = VEC_MUL(dz, rd5)
        vtcharge_rd3 = VEC_MUL(vtcharge, rd3)

        call VEC_ACCUMULATE(ip, particle_pack%pot,
                VEC_MADD(vtcharge,
                         rd,
                         VEC_MADD(rd3,
                                  VEC_aTbPcTdPeTf(dx,vtdip(1),dy,vtdip(2),dz,vtdip(3)),
                                  VEC_MADD(vhalf,
                                           VEC_aTbPcTdPeTf(fd1,vtquad(1),fd2,vtquad(2),fd3,vtquad(3)),
                                           VEC_aTbPcTdPeTf(fd4,vtxyquad,fd5,vtyzquad,fd6,vtzxquad)
                                           )
                                  )
                         )
                        )

        call VEC_ACCUMULATE(ip,particle_pack%ex,
                  VEC_MADD( vthree,
                            VEC_MADD( vhalf,
                                      VEC_aTbPcTdPeTf(VEC_aTbMcTd(vfive_rd7_dx,dx2,dx_rd5,vthree),vtquad(1),
                                                      VEC_MSUB(   vfive_rd7_dx,dy2,dx_rd5)       ,vtquad(2),
                                                      VEC_MSUB(   vfive_rd7_dx,dz2,dx_rd5)       ,vtquad(3)
                                                      ),
                                      VEC_aTbPcTdPeTf(VEC_MSUB(vfive_rd7_dy,dx2,dy_rd5), vtxyquad,
                                                      VEC_MSUB(vfive_rd7_dz,dx2,dz_rd5), vtzxquad,
                                                      vfive_rd7_dxdydz,                  vtyzquad
                                                      ) 
                                      ),
                            VEC_MADD( vtcharge_rd3,
                                      dx,
                                      VEC_aTbPcTdPeTf(fd1,vtdip(1),
                                                      fd4,vtdip(2),
                                                      fd6,vtdip(3))
                                     )
                           )
                          )

        call VEC_ACCUMULATE(ip,particle_pack%ey,
                  VEC_MADD( vthree,
                            VEC_MADD( vhalf,
                                      VEC_aTbPcTdPeTf(VEC_aTbMcTd(vfive_rd7_dy,dy2,dy_rd5,vthree),vtquad(2),
                                                      VEC_MSUB(   vfive_rd7_dy,dx2,dy_rd5)       ,vtquad(1),
                                                      VEC_MSUB(   vfive_rd7_dy,dz2,dy_rd5)       ,vtquad(3)
                                                      ),
                                      VEC_aTbPcTdPeTf(VEC_MSUB(vfive_rd7_dx,dy2,dx_rd5 ),vtxyquad, 
                                                      VEC_MSUB(vfive_rd7_dz,dy2,dz_rd5 ),vtyzquad,
                                                      vfive_rd7_dxdydz,                  vtzxquad
                                                      )
                                      ),
                            VEC_MADD(vtcharge_rd3,
                                     dy,
                                     VEC_aTbPcTdPeTf(fd2,vtdip(2),
                                                     fd4,vtdip(1),
                                                     fd5,vtdip(3))
                                     )
                          )
                         )

        call VEC_ACCUMULATE(ip,particle_pack%ez,
                  VEC_MADD( vthree,
                            VEC_MADD( vhalf,
                                      VEC_aTbPcTdPeTf(VEC_aTbMcTd(vfive_rd7_dz,dz2,dz_rd5,vthree),vtquad(3),
                                                      VEC_MSUB(   vfive_rd7_dz,dy2,dz_rd5)       ,vtquad(2),
                                                      VEC_MSUB(   vfive_rd7_dz,dx2,dz_rd5)       ,vtquad(1)
                                                      ),
                                      VEC_aTbPcTdPeTf(VEC_MSUB(vfive_rd7_dx,dz2,dx_rd5 ),vtzxquad, 
                                                      VEC_MSUB(vfive_rd7_dy,dz2,dy_rd5 ),vtyzquad,
                                                      vfive_rd7_dxdydz,                  vtxyquad
                                                      )
                                      ),
                            VEC_MADD(vtcharge_rd3,
                                     dz,
                                     VEC_aTbPcTdPeTf(fd3,vtdip(3),
                                                     fd5,vtdip(2),
                                                     fd6,vtdip(1))
                                     )
                          )
                         )
      end do
  #else
      real(kfp) :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      do ip = 1, np
        dx = delta(ip, 1)
        dy = delta(ip, 2)
        dz = delta(ip, 3)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (dx >= spatial_interaction_cutoff(1) .or. dy >= spatial_interaction_cutoff(2) .or. dz >= spatial_interaction_cutoff(3)) &
          cycle
        #endif

        r  = sqrt(dist2(ip) + eps2)
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

        particle_pack%pot(ip) = particle_pack%pot(ip) &
              + t%charge*rd                                           &  !  monopole term
              + (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3         &  !  dipole
              + half*(fd1*t%quad(1) + fd2*t%quad(2) + fd3*t%quad(3))  &  !  quadrupole
              +       fd4*t%xyquad  + fd5*t%yzquad  + fd6*t%zxquad

        particle_pack%ex(ip) = particle_pack%ex(ip) &
                  + t%charge*dx*rd3                                  &  ! monopole term
                  + fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3)       &  ! dipole term
                  + three * (                                        &  ! quadrupole term
                    half * (                                         &
                        ( five*dx3   *rd7 - three*dx*rd5 )*t%quad(1) &
                      + ( five*dx*dy2*rd7 -       dx*rd5 )*t%quad(2) &
                      + ( five*dx*dz2*rd7 -       dx*rd5 )*t%quad(3) &
                    )                                                &
                    + ( five*dy*dx2  *rd7 - dy*rd5 )*t%xyquad        &
                    + ( five*dz*dx2  *rd7 - dz*rd5 )*t%zxquad        &
                    + ( five*dx*dy*dz*rd7          )*t%yzquad        &
                    )

        particle_pack%ey(ip) = particle_pack%ey(ip) &
                  + t%charge*dy*rd3                                  &
                  + fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3)       &
                  + three * (                                        &
                    half * (                                         &
                        ( five*dy3*rd7    - three*dy*rd5 )*t%quad(2) &
                      + ( five*dy*dx2*rd7 -       dy*rd5 )*t%quad(1) &
                      + ( five*dy*dz2*rd7 -       dy*rd5 )*t%quad(3) &
                    )                                                &
                    + ( five*dx*dy2  *rd7 - dx*rd5 )*t%xyquad        &
                    + ( five*dz*dy2  *rd7 - dz*rd5 )*t%yzquad        &
                    + ( five*dx*dy*dz*rd7          )*t%zxquad        &
                    )

        particle_pack%ez(ip) = particle_pack%ez(ip) &
                  + t%charge*dz*rd3                                  &
                  + fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1)       &
                  + three * (                                        &
                    half * (                                         &
                      + ( five*dz3   *rd7 - three*dz*rd5 )*t%quad(3) &
                      + ( five*dz*dy2*rd7 -       dz*rd5 )*t%quad(2) &
                      + ( five*dz*dx2*rd7 -       dz*rd5 )*t%quad(1) &
                    )                                                &
                    + ( five*dx*dz2  *rd7 - dx*rd5 )*t%zxquad        &
                    + ( five*dy*dz2  *rd7 - dy*rd5 )*t%yzquad        &
                    + ( five*dx*dy*dz*rd7          )*t%xyquad        &
                    )
      end do
  #endif
    end subroutine calc_force_coulomb_3D


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D(delta, dist2, particle_pack, t, eps2)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      do ip = 1, np
        dx = delta(ip,1)
        dy = delta(ip,2)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (dx >= spatial_interaction_cutoff(1) .or. dy >= spatial_interaction_cutoff(2)) cycle
        #endif

        dx2 = dx *dx
        dy2 = dy *dy
        dx3 = dx2*dx
        dy3 = dy2*dy

        d2 = dx2 + dy2 + eps2

        rd2 = one/d2
        rd4 = rd2*rd2
        rd6 = rd4*rd2

        particle_pack%pot(ip) = particle_pack%pot(ip) &
            - half*t%charge*log(d2)               & !  monopole term
            + (dx*t%dip(1) + dy*t%dip(2))*rd2     & !  dipole
            + half*(t%quad(1)*(dx2*rd4 - rd2)     & ! quadrupole
            +       t%quad(2)*(dy2*rd4 - rd2))    &
            + t%xyquad*dx*dy*rd4

        particle_pack%ex(ip) = particle_pack%ex(ip) &
            + t%charge*dx*rd2                              & ! monopole
            + t%dip(1)*(two*dx2  *rd4 - rd2)               & ! dipole
            + t%dip(2)* two*dx*dy*rd4                      &
            + t%quad(1)*(four *dx3    *rd6 - three*dx*rd4) & ! quadrupole
            + t%quad(2)*(four *dx *dy2*rd6 -       dx*rd4) &
            + t%xyquad *(eight*dx2*dy *rd6 -   two*dy*rd4)

        particle_pack%ey(ip) = particle_pack%ey(ip) &
            + t%charge*dy*rd2                              & ! monopole
            + t%dip(2)*(two*dy2  *rd4 - rd2)               & ! dipole
            + t%dip(1)* two*dx*dy*rd4                      &
            + t%quad(2)*(four *dy3    *rd6 - three*dy*rd4) & ! quadrupole
            + t%quad(1)*(four *dy *dx2*rd6 -       dy*rd4) &
            + t%xyquad *(eight*dy2*dx *rd6 -   two*dx*rd4)
      end do
    end subroutine calc_force_coulomb_2D


    !>
    !> CALC_FORCE_LJ
    !>
    !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
    !> shifted by the lattice vector vbox
    !> results are returned exyz, phi
    !>
    subroutine calc_force_LJ(delta, dist2, particle_pack, t, aii2)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: aii2

      real(kfp) :: flj, epsc2, aii2_r2,aii4_r4, r, fljrd
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      ! epsc should be > a_ii to get evenly spaced ions
      epsc2 = 0.8_kfp * aii2

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        ! Force is repulsive up to and just beyond aii
        !aii2_r2 = aii2 / max(dist2(ip), epsc2)

        !aii4_r4 = aii2_r2*aii2_r2

        !flj = two*(aii4_r4*aii4_r4) - aii4_r4

        !  forces
        r     = sqrt(dist2(ip))
        !fljrd = flj/r
        fljrd = 1.0_kfp/r

        particle_pack%ex(ip) = particle_pack%ex(ip) + delta(ip,1) * fljrd
        particle_pack%ey(ip) = particle_pack%ey(ip) + delta(ip,2) * fljrd
        particle_pack%ez(ip) = particle_pack%ez(ip) + delta(ip,3) * fljrd
      end do
    end subroutine calc_force_LJ


    !>
    !> Calculates 3D Coulomb interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_coulomb_3D_direct(delta, dist2, particle_pack, t, eps2)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

  #if defined(__TOS_BGQ__) && ( defined(__IBMC__) || defined(__IBMCPP__) )
      VECTOR(REAL(kfp)) :: rd,r,rd3tcharge,rdtcharge, vzero, vdist2, vone, vtcharge, veps2, include_mask
      integer(kind_particle) :: ip, np
      
      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        VECTOR(REAL(kfp)) :: vcut(3)
        
        do ip=1,3
          vcut(ip) = VEC_SPLATS(spatial_interaction_cutoff(ip))
        end do
      #endif      

      vzero    = VEC_SPLATS(0._kfp)
      vone     = VEC_SPLATS(one)
      vtcharge = VEC_SPLATS(t%charge)
      veps2    = VEC_SPLATS(eps2)
      
      np = 8*size(dist2, kind = kind_particle)

      !ibm* assert(nodeps)
      do ip=0,np-8,32
          vdist2 = VEC_LDA(ip, dist2)

          r  = VEC_SWSQRT_NOCHK( VEC_ADD(vdist2, veps2) )
          include_mask = VEC_CMPGT(vdist2, vzero)
          #ifndef NO_SPATIAL_INTERACTION_CUTOFF
            include_mask = VEC_AND( VEC_AND( VEC_AND( 
              VEC_CMPLT(VEC_LDA(      ip, delta), vcut(1))  , &
              VEC_CMPLT(VEC_LDA(   np+ip, delta), vcut(2)) ), &
              VEC_CMPLT(VEC_LDA(np+np+ip, delta), vcut(3)) ), &
              include_mask )
          #endif
          rd = VEC_SEL(vzero, VEC_SWDIV_NOCHK(vone, r), include_mask)
          rdtcharge  = VEC_MUL(vtcharge, rd)
          call VEC_ACCUMULATE(ip, particle_pack%pot, rdtcharge)
          rd3tcharge = VEC_MUL3(rdtcharge,rd,rd)
          call VEC_MACCUMULATE(ip, particle_pack%ex, rd3tcharge, VEC_LDA(      ip, delta))
          call VEC_MACCUMULATE(ip, particle_pack%ey, rd3tcharge, VEC_LDA(   np+ip, delta))
          call VEC_MACCUMULATE(ip, particle_pack%ez, rd3tcharge, VEC_LDA(np+np+ip, delta))
      end do
  #else
      np = size(dist2, kind = kind_particle)

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        if (dist2(ip) > 0.0_kfp) then
          r         = sqrt(dist2(ip) + eps2)
          rd        = one/r
          rd3charge = t%charge*rd*rd*rd

          particle_pack%pot(ip) = particle_pack%pot(ip) + t%charge*rd
          particle_pack%ex(ip) = particle_pack%ex(ip) + rd3charge*delta(ip,1)
          particle_pack%ey(ip) = particle_pack%ey(ip) + rd3charge*delta(ip,2)
          particle_pack%ez(ip) = particle_pack%ez(ip) + rd3charge*delta(ip,3)
        end if
      end do
  #endif
    end subroutine calc_force_coulomb_3D_direct


    !>
    !> Calculates 2D Coulomb interaction of particle p with tree node inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exy, phi
    !> Unregularized force law is:
    !>   Phi = -2q log R
    !>   Ex = -dPhi/dx = 2 q x/R^2 etc
    !>
    subroutine calc_force_coulomb_2D_direct(delta, dist2, particle_pack, t, eps2)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: eps2

      real(kfp) :: rd2charge, dx, dy, d2
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      do ip = 1, np
        dx = delta(ip,1)
        dy = delta(ip,2)

        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if (dx >= spatial_interaction_cutoff(1) .or. dy >= spatial_interaction_cutoff(2)) cycle
        #endif

        d2 = dx * dx + dy * dy
        if (d2 > 0.0_kfp) then
          d2 = d2 + eps2
          rd2charge = t%charge/d2

          particle_pack%pot(ip) = particle_pack%pot(ip) - half*t%charge*log(d2)
          particle_pack%ex(ip) = particle_pack%ex(ip) + rd2charge * delta(ip,1)
          particle_pack%ey(ip) = particle_pack%ey(ip) + rd2charge * delta(ip,2)
        end if
      end do
    end subroutine calc_force_coulomb_2D_direct


    !>
    !> Calculates 3D Kelbg interaction of particle p with particle inode
    !> that is shifted by the lattice vector vbox
    !> results are returned in exyz, phi
    !>
    subroutine calc_force_kelbg_3D_direct(delta, dist2, particle_pack, t, kelbg_invsqrttemp)
      implicit none

      real(kfp), intent(in) :: delta(:,:)
      real(kfp), intent(in) :: dist2(:)
      type(t_particle_pack), intent(inout) :: particle_pack
      type(t_tree_node_interaction_data), intent(in) :: t
      real(kfp), intent(in) :: kelbg_invsqrttemp

      real(kfp) :: rd,r,rd3
      real(kfp), parameter :: sqrtpi = sqrt(acos(-1.0_8))
      real(kfp) :: ome, rol, lambda, q, fprefac
      integer(kind_particle) :: ip, np

      np = size(dist2, kind = kind_particle)

      q = t%charge

      do ip = 1, np
        #ifndef NO_SPATIAL_INTERACTION_CUTOFF
        if ( &
          delta(ip,1) >= spatial_interaction_cutoff(1) .or. &
          delta(ip,2) >= spatial_interaction_cutoff(2) .or. &
          delta(ip,3) >= spatial_interaction_cutoff(3) &
        ) cycle
        #endif

        if (dist2(ip) > 0.0_kfp) then
          ! TODO: lambda must be adjusted depending on mass and temperature of interacting partners - currently it is fixed for electron-proton interactions
          if (particle_pack%q(ip) * q < 0.) then
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

          r   = sqrt(dist2(ip))
          rd  = one / r
          rd3 = rd*rd*rd
          rol = r  / lambda        !< "r over lambda"
          ome = 1  - exp(-rol*rol) !< "one minus exp(stuff)"

          ! potential
          particle_pack%pot(ip) = particle_pack%pot(ip) + q * rd  * (ome + sqrtpi*rol*(1-erf(rol)))
          !  forces
          fprefac = q * rd3 * ome
          particle_pack%ex(ip) = particle_pack%ex(ip) + fprefac * delta(ip,1)
          particle_pack%ey(ip) = particle_pack%ey(ip) + fprefac * delta(ip,2)
          particle_pack%ez(ip) = particle_pack%ez(ip) + fprefac * delta(ip,3)
        end if
      end do
    end subroutine calc_force_kelbg_3D_direct
end module
