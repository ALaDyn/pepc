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
!> Contains mapper functions for space-filling curves
!>
module module_spacefilling
      use module_pepc_kinds
      implicit none

      integer(kind_key), parameter, public :: KEY_INVALID = 0_kind_key

      interface is_ancestor_of
        module procedure is_ancestor_of, is_ancestor_of_with_level
      end interface is_ancestor_of

      interface coord_to_intcoord
        module procedure coord_to_intcoord, coord_to_intcoord_with_refinement_length
      end interface coord_to_intcoord

      interface is_ancestor_of_particle
        module procedure is_ancestor_of_particle, is_ancestor_of_particle_with_level
      end interface is_ancestor_of_particle

      contains

        !>
        !> calculates level from key by finding position of placeholder bit
        !>
        elemental function level_from_key(key)
          use treevars, only: idim
          implicit none
          integer(kind_key), intent(in) :: key
          integer(kind_level) :: level_from_key

          ! using log_{2**idim}(key):
          ! level_from_key = int( log(1._8*key) / log((2._8)**idim))
          ! counting leading zeros (faster):
          level_from_key = int((bit_size(key) - leadz(key) - 1_kind_key) / idim, kind_level)
        end function


        !>
        !> returns the key of a key's parent
        !>
        DEBUG_ELEMENTAL function parent_key_from_key(key)
          use module_debug
          implicit none
          integer(kind_key), intent(in) :: key
          integer(kind_key) :: parent_key_from_key

          DEBUG_ASSERT(level_from_key(key) > 0)

          parent_key_from_key = shift_key_by_level(key, -1_kind_level)
        end function parent_key_from_key


        !>
        !> returns the child number of `key` with respect to its parent
        !>
        DEBUG_ELEMENTAL function child_number_from_key(key)
          use treevars, only: idim
          use module_debug
          implicit none

          integer(kind_key), intent(in) :: key
          integer(kind_byte) :: child_number_from_key

          DEBUG_ASSERT(level_from_key(key) > 0)

          child_number_from_key = int(ibits(key, 0, idim), kind_byte)
        end function child_number_from_key


        !>
        !> returns `key` shifted up (negative argument) or down by a number
        !> of levels
        !>
        DEBUG_PURE function shift_key_by_level(key, lvl)
          use treevars
          use module_debug
          implicit none

          integer(kind_key), intent(in) :: key
          integer(kind_level), intent(in) :: lvl
          integer(kind_key) :: shift_key_by_level

          DEBUG_ASSERT(0 <= level_from_key(key) + lvl)
          DEBUG_ASSERT(level_from_key(key) + lvl <= maxlevel)

          shift_key_by_level = ishft(key, idim * lvl)
        end function shift_key_by_level


        !>
        !> returns the key for child `n` of a node with key `key`
        !>
        DEBUG_PURE function child_key_from_parent_key(key, n)
          implicit none

          integer(kind_key), intent(in) :: key
          integer, intent(in) :: n
          integer(kind_key) :: child_key_from_parent_key

          child_key_from_parent_key = shift_key_by_level(key, 1_kind_level) + n
        end function child_key_from_parent_key


        !>
        !> checks whether `ka` is an ancestor of `kc`
        !>
        DEBUG_PURE function is_ancestor_of(ka, kc)
          implicit none
          logical :: is_ancestor_of
          integer(kind_key), intent(in) :: ka, kc

          integer(kind_level) :: la, lc

          la = level_from_key(ka)
          lc = level_from_key(kc)
          is_ancestor_of = is_ancestor_of_with_level(ka, la, kc, lc)
        end function


        !>
        !> checks whether `ka` is an ancestor of `kc`, respective levels are `la` and `lc`
        !>
        DEBUG_PURE function is_ancestor_of_with_level(ka, la, kc, lc)
          use module_debug
          implicit none
          logical :: is_ancestor_of_with_level
          integer(kind_key), intent(in) :: ka, kc
          integer(kind_level), intent(in) :: la, lc

          DEBUG_ASSERT(level_from_key(ka) == la)
          DEBUG_ASSERT(level_from_key(kc) == lc)

          is_ancestor_of_with_level = lc >= la
          if (.not. is_ancestor_of_with_level) return
          is_ancestor_of_with_level = ka == shift_key_by_level(kc, la - lc)
        end function


        !>
        !> checks whether `ka` is an ancestor of `kp` (which must be at highest tree level `maxlevel`, i.e. a particle key)
        !>
        DEBUG_PURE function is_ancestor_of_particle(ka, kp)
          use treevars, only: maxlevel
          use module_debug
          implicit none
          logical :: is_ancestor_of_particle
          integer(kind_key), intent(in) :: ka, kp

          integer(kind_level) :: la

          DEBUG_ASSERT(level_from_key(kp) == maxlevel)

          la = level_from_key(ka)
          is_ancestor_of_particle = is_ancestor_of_with_level(ka, la, kp, maxlevel)
        end function


        !>
        !> checks whether `ka` at level `la` is an ancestor of `kp`
        !> (which must be at highest tree level `maxlevel`, i.e. a particle key)
        !>
        DEBUG_PURE function is_ancestor_of_particle_with_level(ka, la, kp)
          use treevars, only: maxlevel
          use module_debug
          implicit none
          logical :: is_ancestor_of_particle_with_level
          integer(kind_key), intent(in) :: ka
          integer(kind_level), intent(in) :: la
          integer(kind_key), intent(in) :: kp

          DEBUG_ASSERT(level_from_key(ka) == la)
          DEBUG_ASSERT(level_from_key(kp) == maxlevel)

          is_ancestor_of_particle_with_level = is_ancestor_of_with_level(ka, la, kp, maxlevel)
        end function


        function coord_to_intcoord(b, x) result(intcoord)
          use module_box, only: t_box
          use treevars, only: idim, maxlevel
          implicit none

          integer(kind_key) :: intcoord(idim)

          type(t_box), intent(in) :: b
          real(kind_physics), intent(in) :: x(3)

          real(kind_physics) :: s(3)

          s = b%boxsize / 2_kind_key**maxlevel       ! refinement length
          intcoord = coord_to_intcoord_with_refinement_length(b, s, x)
        end function coord_to_intcoord


        function coord_to_intcoord_with_refinement_length(b, s, x) result(intcoord)
          use module_box, only: t_box
          use treevars, only: idim
          implicit none

          integer(kind_key) :: intcoord(idim)

          type(t_box), intent(in) :: b
          real(kind_physics), intent(in) :: s(3)
          real(kind_physics), intent(in) :: x(3)

          intcoord(:) = int(( x(1:idim) - b%boxmin(1:idim) ) / s(1:idim), kind = kind_key) ! partial keys
        end function coord_to_intcoord_with_refinement_length


        !>
        !> calculates keys from local particles (faster than per-particle call to coord_to_key())
        !>
        subroutine compute_particle_keys(b, particles)
          use treevars, only: idim, maxlevel
          use module_pepc_types, only: t_particle
          use module_box, only: t_box
          use module_debug
          implicit none

          type(t_box), intent(in) :: b
          type(t_particle), intent(inout) :: particles(:)

          real(kind_physics) :: s(3)
          integer(kind_particle) :: j, nl

          nl = ubound(particles, 1)
          s = b%boxsize / 2_kind_key**maxlevel       ! refinement length

          ! construct particle keys
          select case (idim)
            case (1) ! fallback to Z-curve
              do j = 1, nl
                particles(j)%key = intcoord_to_key_morton1D(coord_to_intcoord(b, s, particles(j)%x))
              end do
            case (2) ! 2D hilbert curve
              do j = 1, nl
                particles(j)%key = intcoord_to_key_hilbert2D(coord_to_intcoord(b, s, particles(j)%x))
              end do
            case (3) ! 3D hilbert curve
              do j = 1, nl
                particles(j)%key = intcoord_to_key_hilbert3D(coord_to_intcoord(b, s, particles(j)%x))
              end do
            case default
              DEBUG_ERROR(*, "Key generation implemented for 1D, 2D and 3D")
          end select
        end subroutine compute_particle_keys


        !>
        !> calculates key from particle coordinate as vector on top level
        !>
        function coord_to_key(b, x)
          use treevars, only : idim
          use module_box
          use module_debug
          implicit none

          integer(kind_key) :: coord_to_key
          type(t_box), intent(in) :: b
          real(kind_physics), intent(in) :: x(3)

          ! construct particle keys
          select case (idim)
            case (1) ! fallback to Z-curve
              coord_to_key = intcoord_to_key_morton1D(coord_to_intcoord(b, x))
            case (2) ! 2D hilbert curve
              coord_to_key = intcoord_to_key_hilbert2D(coord_to_intcoord(b, x))
            case (3) ! 3D hilbert curve
              coord_to_key = intcoord_to_key_hilbert3D(coord_to_intcoord(b, x))
            case default
              coord_to_key = -1
              DEBUG_ERROR(*, "Key generation implemented for 1D, 2D and 3D")
          end select
        end function coord_to_key


        !>
        !> (Morton-)Z-curve
        !> construct keys by interleaving coord bits and add placeholder bit
        !> note use of 64-bit constants to ensure correct arithmetic
        !>
        function intcoord_to_key_morton1D(ic)
          use treevars
          use module_debug
          implicit none
          integer(kind_key), intent(in) :: ic(1)
          integer(kind_key) :: intcoord_to_key_morton1D
          integer(kind_level) :: i

          DEBUG_ASSERT(idim == 1)

          ! set placeholder bit
          intcoord_to_key_morton1D = 1_kind_key

          ! key generation
          do i = maxlevel - 1_kind_level, 0_kind_level, -1_kind_level
            intcoord_to_key_morton1D = ior(ishft(intcoord_to_key_morton1D, 1), ibits(ic(1), i, 1_kind_key))
          end do
        end function intcoord_to_key_morton1D


        !>
        !> Hilbert-curve,
        !>
        !> algorithm from
        !>
        !> Kamata, S.-I.; Eason, R.O.; Bandou, Y.; ,
        !> "A new algorithm for N-dimensional Hilbert scanning",
        !> Image Processing, IEEE Transactions on , vol.8, no.7, pp.964-973, Jul 1999
        !> doi: 10.1109/83.772242
        !> http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=772242&isnumber=16772
        !>
        function intcoord_to_key_hilbert2D(ic)
          use treevars
          use module_debug
          implicit none
          integer(kind_key), intent(in) :: ic(2)
          integer(kind_key) :: intcoord_to_key_hilbert2D
          integer(kind_level) :: i
          integer :: j

          integer, parameter :: CI2(0:3)    = [0,3,1,2] ! 2D - inverse hilbert cell
          integer, parameter :: G2(0:3,0:1) = reshape([3,0,0,3,0,0,0,3],shape(G2))
          integer(kind_key) :: horder ! order of the hilbert cell C, has to be kind_key to assuage xlf s strict interpretation of ior-parameter
          integer :: exchange, reverse
          integer(kind_key) :: itemp(2), change
          integer(kind_key) :: cval

          DEBUG_ASSERT(idim == 2)

          ! copy, because construction alters original values
          itemp=ic

          ! set placeholder bit
          intcoord_to_key_hilbert2D = 1

          ! key generation
          do i=maxlevel-1_kind_level,0_kind_level,-1_kind_level

            cval = 0_kind_key

            do j=2,1,-1
              cval = ior(ishft(cval, 1), ibits(itemp(j), i, 1_kind_key))
              itemp(j) = ibclr(itemp(j), i)
            end do

            ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
            horder = CI2(cval)
            exchange = G2(horder,0)
            reverse = G2(horder,1)

            ! appending H-order to hkey
            intcoord_to_key_hilbert2D = ior(ishft(intcoord_to_key_hilbert2D, 2), horder)

            ! transform partial curve with the gene rules for the next level step
            ! exchange
            select case (exchange)
              case (3)
                change = itemp(1)
                itemp(1) = itemp(2)
                itemp(2) = change
            end select

            ! reverse
            do j=1,2
              if (btest(reverse, j - 1)) itemp(j) = iand(not(itemp(j)), 2_kind_key**(i) - 1)
            end do
          end do

        end function intcoord_to_key_hilbert2D


        function intcoord_to_key_hilbert3D(ic)
          use treevars
          use module_debug
          implicit none
          integer(kind_key), intent(in) :: ic(3)
          integer(kind_key) :: intcoord_to_key_hilbert3D
          integer(kind_level) :: i
          integer :: j

          integer, parameter :: CI3(0:7)    = [0,1,3,2,7,6,4,5] ! 3D - inverse hilbert cell
          integer, parameter :: G3(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G3))     ! 3D - hilbert gene
          integer(kind_key) :: horder ! order of the hilbert cell C, has to be kind_key to assuage xlf s strict interpretation of ior-parameter
          integer :: exchange, reverse
          integer(kind_key) :: itemp(3), change
          integer(kind_key) :: cval

          DEBUG_ASSERT(idim == 3)

          ! copy, because construction alters original values
          itemp=ic

          ! set placeholder bit
          intcoord_to_key_hilbert3D = 1

          ! key generation
          do i=maxlevel-1_kind_level,0_kind_level,-1_kind_level

            cval = 0_kind_key

            do j=3,1,-1
              cval = ior(ishft(cval, 1), ibits(itemp(j), i, 1_kind_key))
              itemp(j) = ibclr(itemp(j), i)
            end do

            ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
            horder = CI3(cval)
            exchange = G3(horder,0)
            reverse = G3(horder,1)

            ! appending H-order to hkey
            intcoord_to_key_hilbert3D = ior(ishft(intcoord_to_key_hilbert3D, 3), horder)

            ! transform partial curve with the gene rules for the next level step
            ! exchange
            select case (exchange)
              case (5)
                change = itemp(1)
                itemp(1) = itemp(3)
                itemp(3) = change
              case (6)
                change = itemp(2)
                itemp(2) = itemp(3)
                itemp(3) = change
            end select

            ! reverse
            do j=1,3
              if (btest(reverse, j - 1)) itemp(j) = iand(not(itemp(j)), 2_kind_key**(i) - 1)
            end do
          end do
        end function intcoord_to_key_hilbert3D


        !>
        !> calculates particle coordinate as vector from key
        !>
        subroutine key_to_coord(b, key, x)
          use treevars, only : idim, maxlevel
          use module_box
          use module_debug
          implicit none

          type(t_box), intent(in) :: b
          integer(kind_key), intent(in) :: key
          real(kind_physics), intent(inout) :: x(3)

          integer(kind_key) :: ic(idim)
          real(kind_physics) :: s(idim)

          ! construct particle coordiantes
          select case (idim)
            case (1) ! Z-curve
              call key_to_intcoord_morton(key, ic)
            case (2, 3) ! Hilbert curve (original pattern)
              call key_to_intcoord_hilbert(key, ic)
            case default
              DEBUG_ERROR(*, "Key generation implemented for 1D, 2D and 3D")
          end select

          s = b%boxsize(1:idim) / 2_kind_key**maxlevel       ! refinement length

          x(1:idim) = (real(ic, kind_physics) + 0.5_kind_physics) * s + b%boxmin(1:idim)
        end subroutine key_to_coord


        !>
        !> (Morton-)Z-curve
        !> keys were constructed by interleaving coord bits and adding placeholder bit
        !> input key must be right-adjusted and must contain a placeholder bit
        !>
        subroutine key_to_intcoord_morton(key, ic)
          use treevars
          use module_debug
          implicit none
          integer(kind_key), intent(in) :: key
          integer(kind_key), intent(out) :: ic(1)
          integer(kind_level) :: i, lev

          DEBUG_ASSERT(idim == 1)

          lev = level_from_key(key)

          ic = 0_kind_key

          do i = 0_kind_level, lev - 1_kind_level
            ic(1) = ior(ic(1), ishft(ibits(key, i, 1), maxlevel - lev + i))
          end do
        end subroutine key_to_intcoord_morton


        !>
        !> inverse Hilbert mapping
        !> input key as vector must be right-adjusted and must contain a placeholder bit
        !>
        subroutine key_to_intcoord_hilbert(key, ic)
          use treevars
          use module_debug
          implicit none
          integer(kind_key), intent(out) :: ic(idim)
          integer(kind_key), intent(in) :: key

          integer(kind_key) :: change, horder, cval
          integer :: exchange, reverse
          integer(kind_level) :: i, lev
          integer(kind_dim) :: j, k

          integer, parameter :: C2(0:3) = [0,2,3,1] ! 2D - hilbert cell
          integer, parameter :: G2(0:3,0:1) = reshape([3,0,0,3,0,0,0,3],shape(G2)) ! 2D - hilbert gene
          integer, parameter :: C3(0:7)    = [0,1,3,2,6,7,5,4] ! 3D - hilbert cell
          integer, parameter :: G3(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G3)) ! 3D - hilbert gene

          DEBUG_ASSERT(idim == 2 .or. idim == 3)

          lev = level_from_key(key)

          ic = 0_kind_key

          do i=0_kind_level,lev-1_kind_level
            horder = ibits(key,idim*i,idim)

            if (idim == 2) then
              cval = C2(horder)
              reverse = G2(horder,1)
              exchange = G2(horder,0)
            else
              cval = C3(horder)
              reverse = G3(horder,1)
              exchange = G3(horder,0)
            end if

            if (i>0) then
              ! reverse
              do j=1_kind_dim,idim
                if (btest(reverse, j - 1)) ic(j) = iand(not(ic(j)), 2_kind_key**(i) - 1)
              end do

              ! exchange
              do j=1_kind_dim,idim-1_kind_dim
                if (.not. btest(exchange, j - 1)) cycle

                do k=j+1_kind_dim,idim
                  if (btest(exchange, k - 1)) then
                    change = ic(j)
                    ic(j) = ic(k)
                    ic(k) = change
                  end if
                end do
              end do
            end if

            do j=1_kind_dim,idim
              ic(j) = ior(ishft(ibits(cval, j-1, 1), i), ic(j))
            end do
          end do

          do j=1_kind_dim,idim
            ic(j) = ishft(ic(j), maxlevel-lev)
          end do
        end subroutine key_to_intcoord_hilbert
end module module_spacefilling
