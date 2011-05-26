!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains mapper functions for space-filling curves
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_spacefilling
      use treevars
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: curve_type = 0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates keys form local particles (faster than per-particle call to coord_to_key())
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine compute_particle_keys(local_key)
          implicit none
          integer*8, intent(out) :: local_key(nppm)
          integer*8, dimension(nppm) :: ix, iy, iz
          real*8 :: s
          integer :: j

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

          ix(1:npp) = int(( x(1:npp) - xmin )/s)           ! partial keys
          iy(1:npp) = int(( y(1:npp) - ymin )/s)           !
          iz(1:npp) = int(( z(1:npp) - zmin )/s)

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              do j = 1,npp
                 local_key(j) = intcoord_to_key_morton(ix(j), iy(j), iz(j))
              end do

            case (1) ! Hilbert curve (original pattern)

             ! for all particles
             do j=1,npp
                 local_key(j) = intcoord_to_key_hilbert(ix(j), iy(j), iz(j))
             end do

          end select

          if (domain_debug) then
             write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8,2f12.4))') 'Particle list before key sort:', &
                  '  key,             label   coords     q ', &
                  (local_key(j),pelabel(j),x(j),y(j),z(j),ix(j),iy(j),iz(j),q(j),work(j),j=1,npp)
             !    write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8))') '(last 10):', &
             !         '  key,                  label        coords              q ', &
             !         (local_key(i),pelabel(i),x(i),y(i),z(i),ix(i),iy(i),iz(i),q(i),work(i),i=max(1,npp-10),npp)

             write(ipefile,'(/)')
          endif

        end subroutine compute_particle_keys


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates key from particle coordiante
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function coord_to_key(x, y, z)
          implicit none
          integer*8 :: coord_to_key
          real*8, intent(in) :: x, y, z
          integer*8 :: ix, iy, iz
          real*8 :: s

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          ix = int(( x - xmin )/s)           ! partial keys
          iy = int(( y - ymin )/s)           !
          iz = int(( z - zmin )/s)

          ! construct particle keys
          select case (curve_type)
            case (0) ! Z-curve
              coord_to_key = intcoord_to_key_morton(ix, iy, iz)
            case (1) ! Hilbert curve (original pattern)
              coord_to_key = intcoord_to_key_hilbert(ix, iy, iz)
          end select

        end function coord_to_key


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> calculates particle coordiante from key
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_coord(key, x, y, z)
          implicit none
          integer*8, intent(in) :: key
          real*8, intent(out) :: x, y, z
          integer*8 :: ix, iy, iz, tmp
          real*8 :: s

          tmp = key
          ! shift left until placeholder bit is in front
          do while ((iand(tmp, iplace) .eq. 0) .and. (tmp .ne. 0))
            tmp = ishft(tmp, 3)
          end do
          ! eliminate placeholder bit
          tmp = iand(tmp, not(iplace))

          ! construct particle coordiantes
          select case (curve_type)
            case (0) ! Z-curve
              call key_to_intcoord_morton(tmp, ix, iy, iz)
            case (1) ! Hilbert curve (original pattern)
              call key_to_intcoord_hilbert(tmp, ix, iy, iz)
          end select

          s=boxsize/2**nlev       ! refinement length

          ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)
          x = (real(ix,kind(1._8)) + 0.5_8) * s + xmin
          y = (real(iy,kind(1._8)) + 0.5_8) * s + ymin
          z = (real(iz,kind(1._8)) + 0.5_8) * s + zmin

          !write(*,'(2O30, 3G15.5)') key, coord_to_key(x, y, z), x, y, z

        end subroutine key_to_coord


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> construct keys by interleaving coord bits and add placeholder bit
        !> note use of 64-bit constants to ensure correct arithmetic
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function intcoord_to_key_morton(ix, iy, iz)
          implicit none
          integer*8, intent(in) :: ix, iy, iz
          integer*8 :: intcoord_to_key_morton
          integer :: i

          !     local_key(j) = iplace + &
          !          SUM( (/ (8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) ),i=0,nlev) /) )
          intcoord_to_key_morton = iplace

          do i=0,nlev
            intcoord_to_key_morton = intcoord_to_key_morton &
                 + 8_8**i*(4_8*ibits( iz,i,1) + 2_8*ibits( iy,i,1 ) + 1_8*ibits( ix,i,1) )
          end do
        end function intcoord_to_key_morton


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> (Morton-)Z-curve
        !> keys were constructed by interleaving coord bits and adding placeholder bit
        !> input key must be adjusted to lowest level and may not contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_intcoord_morton(key, ix, iy, iz)
          implicit none
          integer*8, intent(out) :: ix, iy, iz
          integer*8, intent(in) :: key
          integer :: i

          ix = 0
          iy = 0
          iz = 0

          do i=0,nlev
            ix = ior(ix, ishft(ibits(key, 3*i + 0, 1), i))
            iy = ior(iy, ishft(ibits(key, 3*i + 1, 1), i))
            iz = ior(iz, ishft(ibits(key, 3*i + 2, 1), i))
          end do

          !write(*,'(O24.24," ",O24.24,3(/,B24.24),/)') key, intcoord_to_key_morton(ix, iy, iz), ix, iy, iz

        end subroutine key_to_intcoord_morton


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function intcoord_to_key_hilbert(ix, iy, iz)
          implicit none
          integer*8, intent(in) :: ix, iy, iz
          integer*8 :: intcoord_to_key_hilbert
          integer :: i

		  integer*8, parameter :: CI(0:7)    = [0,1,3,2,7,6,4,5] ! 3D - inverse hilbert cell
          integer*8, parameter :: G(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G))     ! 3D - hilbert gene
		  integer*8 :: horder           ! order of the hilbert cell C
		  integer*8 :: xtemp,ytemp,ztemp,change
		  integer*8 :: cval

	        ! copy, because construction alters original values
	        xtemp=ix
	        ytemp=iy
	        ztemp=iz

	        ! set placeholder bit
	        intcoord_to_key_hilbert = 0!1

	        ! key generation
	        do i=nlev-1,0,-1

	          cval = 4_8*ibits( ztemp,i, 1_8 ) + 2_8*ibits( ytemp,i, 1_8 ) + 1_8*ibits( xtemp,i, 1_8 )
              xtemp = ibclr(xtemp, i)
              ytemp = ibclr(ytemp, i)
              ztemp = ibclr(ztemp, i)

	           ! get H-order (gathering upper bits as in z-mapping and combine to binary number)
	           horder = CI( cval )

	           ! appending H-order to hkey
	           intcoord_to_key_hilbert = ior(ishft(intcoord_to_key_hilbert, 3), horder)

	           ! transform partial curve with the gene rules for the next level step
	           ! exchange
	           select case (G(horder,0))
	             case (5) ! (= 101[zyx]) --> change z and x
	                change = ztemp
	                ztemp  = xtemp
	                xtemp  = change
                 case (6) ! (= 110[zyx]) --> change z and y
                    change = ztemp
                    ztemp  = ytemp
                    ytemp  = change
	           end select

	           ! reverse
	           select case (G(horder,1))
	             case (5) ! (= 101[zyx]) --> reverse z and x
	               ztemp = iand(not(ztemp),2**(i)-1)
	               xtemp = iand(not(xtemp),2**(i)-1)
                 case (6) ! (= 110[zyx]) --> reverse z and y
                   ztemp = iand(not(ztemp),2**(i)-1)
                   ytemp = iand(not(ytemp),2**(i)-1)
	           end select
	        end do

        end function intcoord_to_key_hilbert



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> inverse Hilbert mapping
        !> input key must be adjusted to lowest level and may not contain a placeholder bit
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine key_to_intcoord_hilbert(key, ix, iy, iz)
          implicit none
          integer*8, intent(out) :: ix, iy, iz
          integer*8, intent(in) :: key

		  integer*8 :: change, horder
		  integer*8 :: i

          integer*8, parameter :: C(0:7)    = [0,1,3,2,6,7,5,4] ! 3D - hilbert cell
          integer*8, parameter :: G(0:7,0:1) = reshape([5,6,0,5,5,0,6,5,0,0,0,5,0,0,6,5],shape(G)) ! 3D - hilbert gene
          !integer*8, parameter :: G(0:7,0:1) = reshape([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],shape(G)) ! 3D - hilbert gene

		  iz = 0
		  iy = 0
		  ix = 0

		  do i=0,nlev-1
		     horder = ibits(key,3*i,3)

             select case(G(horder,1))
             case(5)
                ix=iand(not(ix),2**(i)-1)
                iz=iand(not(iz),2**(i)-1)
             case(6)
                iy=iand(not(iy),2**(i)-1)
                iz=iand(not(iz),2**(i)-1)
             end select

             select case(G(horder,0))
             case(5)
                change=ix
                ix=iz
                iz=change
             case(6)
                change=iy
                iy=iz
                iz=change
             end select

             iz=ior(ishft(ibits(C(horder),2,1),i),iz)
             iy=ior(ishft(ibits(C(horder),1,1),i),iy)
             ix=ior(ishft(ibits(C(horder),0,1),i),ix)
		  end do

          !write(*,'(O24.24," ",O24.24,3(/,B24.24),/)') key, intcoord_to_key_hilbert(ix, iy, iz), ix, iy, iz

        end subroutine key_to_intcoord_hilbert



end module module_spacefilling
