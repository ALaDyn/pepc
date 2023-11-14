module module_math

   use encap
!  use module_pepc_kinds
!  use module_pepc_types
   implicit none

contains

   subroutine div(field_grid, divergence, flag)
      implicit none
      type(field_grid_t), intent(in)                   :: field_grid
      integer(kind_particle), intent(in)               :: flag
      integer(kind_particle)                          :: nx, ny, ip, jp, rc, it
      real(kind_particle)                             :: dx, dy, L2x, L2y
      real(kind_particle), allocatable                 :: Ax(:, :), Ay(:, :), derX(:, :), derY(:, :), vx(:), vy(:)!, derZ(:,:)
      real(kind_particle), allocatable, intent(out)     :: divergence(:, :)

      nx = field_grid%n(1)
      ny = field_grid%n(2)
      dx = field_grid%dx(1)
      dy = field_grid%dx(2)

      if (.not. allocated(divergence)) allocate (divergence(nx, ny), stat=rc)
      allocate (Ax(nx, ny), stat=rc)
      allocate (Ay(nx, ny), stat=rc)
!    allocate(Az(nx,ny), stat=rc)
      allocate (derX(nx, ny), stat=rc)
      allocate (derY(nx, ny), stat=rc)
!    allocate(derZ(nx,ny), stat=rc)
      allocate (vx(nx * ny), stat=rc)
      allocate (vy(nx * ny), stat=rc)

!    derX = 0
!    derY = 0
!    derZ = 0
      divergence = 0

      select case (flag)
      case (2)  !  divergence of A
!              vx(:) =   field_grid%p(:)%results%pot!E(1)
!              vy(:) =   field_grid%p(:)%results%pot!E(2)

         do ip = 1, nx
            do jp = 1, ny
               Ax(ip, jp) = cos(field_grid%p(ip + (jp - 1) * ny)%x(1)) * cos(field_grid%p(ip + (jp - 1) * ny)%x(2))!field_grid%p(ip +(jp-1)*ny)%results%pot !vx( ip +(jp-1)*ny )
               Ay(ip, jp) = cos(field_grid%p(ip + (jp - 1) * ny)%x(2)) * cos(field_grid%p(ip + (jp - 1) * ny)%x(1))!field_grid%p(ip +(jp-1)*ny)%results%pot!vy( ip +(jp-1)*ny )
            end do
         end do
      case (3)  !  divergence of B
         vx = field_grid%p(:)%results%B(1)
         vy = field_grid%p(:)%results%B(2)
      case (4)  !  divergence of J
         vx = field_grid%p(:)%results%J(1)
         vy = field_grid%p(:)%results%J(2)

      case default
         !  divergence of A
         vx = field_grid%p(:)%results%A(1)
         vy = field_grid%p(:)%results%A(2)

      end select

      do ip = 2, nx - 1
         do jp = 2, ny - 1
            derX(ip, jp) = .5 * (Ax(ip + 1, jp) - Ax(ip - 1, jp)) / dx
            derY(ip, jp) = .5 * (Ay(ip, jp + 1) - Ay(ip, jp - 1)) / dy
         end do
      end do

      derY(1, 2:ny - 1)  = .5 * (Ay(1, 3:ny)      -     Ay(1, 1:ny - 2)) / dy                          !&
      derY(nx, 2:ny - 1) = .5 * (Ay(nx, 3:ny)     -     Ay(nx, 1:ny - 2)) / dy                         !&
      derY(1:nx, 1)      = .5 * (4 * Ay(1:nx, 2)  - 3 * Ay(1:nx, 1)      - Ay(1:nx, 3)) / dy           !&
      derY(1:nx, ny)     = .5 * (3 * Ay(1:nx, ny) - 4 * Ay(1:nx, ny - 1) + Ay(1:nx, ny - 2)) / dy      !&

      derX(1, 1:ny)      = .5 * (4 * Ax(2, 1:ny)  - 3 * Ax(1, 1:ny)      - Ax(3, 1:ny)) / dx           !&
      derX(nx, 1:ny)     = .5 * (3 * Ax(nx, 1:ny) - 4 * Ax(nx - 1, 1:ny) + Ax(nx - 2, 1:ny)) / dx      !&
      derX(2:nx - 1, 1)  = .5 * (Ax(3:nx, 1)      -     Ax(1:nx - 2, 1)) / dx                          !&
      derX(2:nx - 1, ny) = .5 * (Ax(3:nx, ny)     -     Ax(1:nx - 2, ny)) / dx                         !&

      L2x = 0
      L2y = 0

      divergence = derX + derY

      it = 1

      open (unit=1, file="divergence.dat", action="write", status="replace")
      do ip = 1, nx
         do jp = 1, ny

!           L2x   = L2x + ( derX(ip,jp) + field_grid%p(it)%results%E(1) )**2
!           L2y   = L2y + ( derY(ip,jp) + field_grid%p(it)%results%E(2) )**2
            L2x = L2x + (derX(ip, jp) + sin(field_grid%p(ip + (jp - 1) * ny)%x(1)) * cos(field_grid%p(ip + (jp - 1) * ny)%x(2)))**2
            L2y = L2y + (derY(ip, jp) + sin(field_grid%p(ip + (jp - 1) * ny)%x(2)) * cos(field_grid%p(ip + (jp - 1) * ny)%x(1)))**2
!           L2x   = L2x + ( derX(ip,jp) - 2*vx( ip +(jp-1)*nx ) )**2
!           L2y   = L2y + ( derY(ip,jp) - 2*vy( ip +(jp-1)*nx ) )**2
!            Az(ip,jp)   = field_grid%p( jp +(ip-1)*nx )%results%A(3)
            write (1, *) field_grid%p(ip + (jp - 1) * ny)%results%E(1), field_grid%p(ip + (jp - 1) * ny)%results%E(2), -derX(ip, jp), -derY(ip, jp), divergence(ip, jp)
            it = it + 1

         end do
      end do

      close (unit=1)

      write (*, *) " ==== Gradient Norm L2 Error", L2x, L2y

      deallocate (Ax, Ay, derX, derY, vx, vy)

   end subroutine

   subroutine curl(field_grid, curlX, curlY, curlZ, flag)
      implicit none
      type(field_grid_t), intent(in)                  :: field_grid
      integer(kind_particle), intent(in)              :: flag
      integer(kind_particle)                          :: nx, ny, ip, jp, rc, it
      real(kind_particle)                             :: dx, dy, L2x, L2y
      real(kind_particle), allocatable                :: Ax(:, :), Ay(:, :), Az(:, :), derX(:, :), derY(:, :), vx(:), vy(:), vz(:), B(:, :)!, derZ(:,:)
      real(kind_particle), allocatable, intent(out)   :: curlX(:, :), curlY(:, :), curlZ(:, :)

      nx = field_grid%n(1)
      ny = field_grid%n(2)
      dx = field_grid%dx(1)
      dy = field_grid%dx(2)

      if (.not. allocated(curlX)) allocate (curlX(nx, ny), stat=rc)
      if (.not. allocated(curlY)) allocate (curlY(nx, ny), stat=rc)
      if (.not. allocated(curlZ)) allocate (curlZ(nx, ny), stat=rc)
      allocate (Ax(nx, ny), stat=rc)
      allocate (Ay(nx, ny), stat=rc)
      allocate (Az(nx, ny), stat=rc)
      allocate (derX(nx, ny), stat=rc)
      allocate (derY(nx, ny), stat=rc)
      allocate (B(nx, ny), stat=rc)
!    allocate(derZ(nx,ny), stat=rc)
      allocate (vx(nx), stat=rc)
      allocate (vy(nx), stat=rc)
      allocate (vz(nx), stat=rc)

      derX = 0
      derY = 0
!    derZ = 0
      curlX = 0
      curlY = 0
      curlZ = 0

      select case (flag)
      case (2)  !  divergence of A
         vx = field_grid%p(:)%results%A(1)
         vy = field_grid%p(:)%results%A(2)

         do ip = 1, nx
            do jp = 1, ny

               Ax(ip, jp) = field_grid%p(ip + (jp - 1) * ny)%results%A(1)!cos( field_grid%p( ip +(jp-1)*ny )%x(2) )
               Ay(ip, jp) = field_grid%p(ip + (jp - 1) * ny)%results%A(2)!cos( field_grid%p( ip +(jp-1)*ny )%x(1) )
               B(ip, jp) = field_grid%p(ip + (jp - 1) * ny)%results%B(3)! cos( field_grid%p( ip +(jp-1)*ny )%x(1) )

            end do
         end do

                              !!! x - derivative  2nd order
         do ip = 2, nx - 1
            do jp = 2, ny - 1
               derX(ip, jp) = .5 * (Ay(ip + 1, jp) - Ay(ip - 1, jp)) / dx
               derY(ip, jp) = .5 * (Ax(ip, jp + 1) - Ax(ip, jp - 1)) / dy
            end do
         end do

!                do jp = 2,ny-1
!                    derX(1,jp)           =  .5* ( 4*Ay( 2, jp ) - 3*Ay( 1, jp) - Ay( 3, jp ) )/dx
!                    derX(nx,jp)          =  .5* ( 3*Ay( nx, jp ) - 4*Ay( nx-1, jp ) + Ay( nx-2, jp ) )/dx
!                    derY(1,jp)           =  .5* ( Ax( 1, jp+1 ) - Ax( 1, jp-1 ) )/dy
!                    derY(nx,jp)          =  .5* ( Ax( nx, jp+1 ) - Ax( nx, jp-1 ) )/dy
!                enddo

!                do ip = 2,nx-1
!                    derY(ip,1)           =  .5* ( 4*Ax( ip, 2 ) - 3*Ax( ip, 1 ) - Ax( ip, 3 ) )/dy
!                    derY(ip,ny)          =  .5* ( 3*Ax( ip, ny ) - 4*Ax( ip, ny-1 ) + Ax( ip, ny-2 ) )/dy
!                    derX(ip,1)           =  .5* ( Ay( ip+1, 1 ) - Ay( ip-1, 1 ) )/dx
!                    derX(ip,ny)          =  .5* ( Ay( ip+1, ny ) - Ay( ip-1, ny ) )/dx
!                enddo

         derY(1, 2:ny - 1)  = .5 * (Ax(1, 3:ny)      -     Ax(1, 1:ny - 2)) / dy                       !&
         derY(nx, 2:ny - 1) = .5 * (Ax(nx, 3:ny)     -     Ax(nx, 1:ny - 2)) / dy                      !&
         derY(1:nx, 1)      = .5 * (4 * Ax(1:nx, 2)  - 3 * Ax(1:nx, 1)      - Ax(1:nx, 3)) / dy        !&
         derY(1:nx, ny)     = .5 * (3 * Ax(1:nx, ny) - 4 * Ax(1:nx, ny - 1) + Ax(1:nx, ny - 2)) / dy   !&

         derX(1, 1:ny)      = .5 * (4 * Ay(2, 1:ny)  - 3 * Ay(1, 1:ny)      - Ay(3, 1:ny)) / dx        !&
         derX(nx, 1:ny)     = .5 * (3 * Ay(nx, 1:ny) - 4 * Ay(nx - 1, 1:ny) + Ay(nx - 2, 1:ny)) / dx   !&
         derX(2:nx - 1, 1)  = .5 * (Ay(3:nx, 1)      -     Ay(1:nx - 2, 1)) / dx                       !&
         derX(2:nx - 1, ny) = .5 * (Ay(3:nx, ny)     -     Ay(1:nx - 2, ny)) / dx                      !&

!                derX(1,1)           =  .5* ( 4*Ay( 2, 1 ) - 3*Ay( 1, 1) - Ay( 3, 1 ) )/dx
!                derX(1,ny)          =  .5* ( 4*Ay( 2, ny ) - 3*Ay( 1, ny) - Ay( 3, ny ) )/dx
!                derX(nx,1)          =  .5* ( 3*Ay( nx, 1 ) - 4*Ay( nx-1, 1 ) + Ay( nx-2, 1 ) )/dx
!                derX(nx,ny)         =  .5* ( 3*Ay( nx, ny ) - 4*Ay( nx-1, ny ) + Ay( nx-2, ny ) )/dx

!                derY(1,1)           =  .5* ( 4*Ax( 1, 2 ) - 3*Ax( 1, 1 ) - Ax( 1, 3 ) )/dy
!                derY(1,ny)          =  .5* ( 3*Ax( 1, ny ) - 4*Ax( 1, ny-1 ) + Ax( 1, ny-2 ) )/dy
!                derY(nx,1)          =  .5* ( 4*Ax( nx, 2 ) - 3*Ax( nx, 1 ) - Ax( nx, 3 ) )/dy
!                derY(nx,ny)         =  .5* ( 3*Ax( nx, ny ) - 4*Ax( nx, ny-1 ) + Ax( nx, ny-2 ) )/dy

         curlX = 0
         curlY = 0
         curlZ(1:nx, 1:ny) = derX(1:nx, 1:ny) - derY(1:nx, 1:ny)

         it = 1
         L2x = 0
         open (unit=1, file="curl.dat", action="write", status="replace")
         do jp = 1, ny
            do ip = 1, nx

               L2x = L2x + (curlZ(ip, jp) - B(ip, jp))**2
!                       L2x   = L2x + ( curlZ(ip,jp)  + sin( field_grid%p( ip +(jp-1)*ny )%x(1) ) )**2
!                       L2y   = L2y + ( derY(ip,jp) + field_grid%p(it)%results%E(2) )**2
               !           L2x   = L2x + ( derX(ip,jp) - 2*vx( ip +(jp-1)*nx ) )**2
               !           L2y   = L2y + ( derY(ip,jp) - 2*vy( ip +(jp-1)*nx ) )**2
               !            Az(ip,jp)   = field_grid%p( jp +(ip-1)*nx )%results%A(3)
               write (1, *) B(ip, jp), curlZ(ip, jp)!field_grid%p(ip +(jp-1)*ny)%results%B(3)
               it = it + 1

            end do
         end do

         close (unit=1)

         write (*, *) "Absolute Error Curl Norm L2 ", L2x

      case (3)  !  divergence of B
         vx = field_grid%p(:)%results%B(1)
         vy = field_grid%p(:)%results%B(2)

         do ip = 1, nx
            do jp = 1, ny

!                        Ax(ip,jp)   = vx( jp +(ip-1)*nx )
!                        Ay(ip,jp)   = vy( jp +(ip-1)*nx )

            end do
         end do

                                                !!! x - derivative  2nd order

         do ip = 2, nx - 1
            do jp = 2, ny - 1
               derX(ip, jp) = .5 * (Az(ip + 1, jp) - Az(ip - 1, jp)) / dx
               derY(ip, jp) = .5 * (Az(ip, jp + 1) - Az(ip, jp - 1)) / dy
            end do
         end do

         do jp = 2, ny - 1
            derX(1, jp)  = .5 * (4 * Az(2, jp)  - 3 * Az(1, jp)      - Az(3, jp)) / dy         !&
            derX(nx, jp) = .5 * (3 * Az(nx, jp) - 4 * Az(nx - 1, jp) + Az(nx - 2, jp)) / dy    !&
            derY(1, jp)  = .5 * (Az(1, jp + 1)  -     Az(1, jp - 1)) / dx                      !&
            derY(nx, jp) = .5 * (Az(nx, jp + 1) -     Az(nx, jp - 1)) / dx                     !&
         end do

         do ip = 2, nx - 1
            derY(ip, 1)  = .5 * (4 * Az(ip, 2)  - 3 * Az(ip, 1)      - Az(ip, 3)) / dx         !&
            derY(ip, ny) = .5 * (3 * Az(ip, ny) - 4 * Az(ip, ny - 1) + Az(ip, ny - 2)) / dx    !&
            derX(ip, 1)  = .5 * (Az(ip + 1, 1)  -     Az(ip - 1, 1)) / dy                      !&
            derX(ip, ny) = .5 * (Az(ip + 1, ny) -     Az(ip - 1, ny)) / dy                     !&
         end do

         derX(1, 1)   = .5 * (4 * Az(2, 1)   - 3 * Az(1, 1)       - Az(3, 1)) / dy             !&
         derX(1, ny)  = .5 * (4 * Az(2, ny)  - 3 * Az(1, ny)      - Az(3, ny)) / dy            !&
         derX(nx, 1)  = .5 * (3 * Az(nx, 1)  - 4 * Az(nx - 1, 1)  + Az(nx - 2, 1)) / dy        !&
         derX(nx, ny) = .5 * (3 * Az(nx, ny) - 4 * Az(nx - 1, ny) + Az(nx - 2, ny)) / dy       !&

         derY(1, 1)   = .5 * (4 * Az(1, 2)   - 3 * Az(1, 1)       - Az(1, 3)) / dx             !&
         derY(1, ny)  = .5 * (3 * Az(1, ny)  - 4 * Az(1, ny - 1)  + Az(1, ny - 2)) / dx        !&
         derY(nx, 1)  = .5 * (4 * Az(nx, 2)  - 3 * Az(nx, 1)      - Az(nx, 3)) / dx            !&
         derY(nx, ny) = .5 * (3 * Az(nx, ny) - 4 * Az(nx, ny - 1) + Az(nx, ny - 2)) / dx       !&

         curlX = derY
         curlY = -derX
         curlZ = 0

      case (4)  !  divergence of J
         vx = field_grid%p(:)%results%J(1)
         vy = field_grid%p(:)%results%J(2)

         do ip = 1, nx
            do jp = 1, ny

               Ax(ip, jp) = vx(jp + (ip - 1) * nx)
               Ay(ip, jp) = vy(jp + (ip - 1) * nx)

            end do
         end do

                              !!! x - derivative  2nd order
         derX(2:nx - 1, :) = .5 * (Ay(3:, :)     - Ay(1:nx - 2, :)) / dx                       !&
         derX(1, :)        = .5 * (4 * Ay(2, :)  - 3 * Ay(1, :)      - Ay(3, :)) / dx          !&
         derX(nx, :)       = .5 * (3 * Ay(nx, :) - 4 * Ay(nx - 1, :) + Ay(nx - 2, :)) / dx     !&
                !!! y - derivative   2nd order
         derY(:, 2:ny - 1) = .5 * (Ax(:, 3:)     - Ax(:, 1:ny - 2)) / dy                       !&
         derY(:, 1)        = .5 * (4 * Ax(:, 2)  - 3 * Ax(:, 1)      - Ax(:, 3)) / dy          !&
         derY(:, ny)       = .5 * (3 * Ax(:, ny) - 4 * Ax(:, ny - 1) + Ax(:, ny - 2)) / dy     !&

         curlX = 0
         curlY = 0
         curlZ = derX - derY

      case default
         !  divergence of A
         vx = field_grid%p(:)%results%A(1)
         vy = field_grid%p(:)%results%A(2)
         do ip = 1, nx
            do jp = 1, ny

               Ax(ip, jp) = vx(jp + (ip - 1) * nx)
               Ay(ip, jp) = vy(jp + (ip - 1) * nx)

            end do
         end do

                              !!! x - derivative  2nd order
         derX(2:nx - 1, :) = .5 * (Ay(3:, :)     -     Ay(1:nx - 2, :)) / dx                   !&
         derX(1, :)        = .5 * (4 * Ay(2, :)  - 3 * Ay(1, :)      - Ay(3, :)) / dx          !&
         derX(nx, :)       = .5 * (3 * Ay(nx, :) - 4 * Ay(nx - 1, :) + Ay(nx - 2, :)) / dx     !&
                !!! y - derivative   2nd order
         derY(:, 2:ny - 1) = .5 * (Ax(:, 3:)     -     Ax(:, 1:ny - 2)) / dy                   !&
         derY(:, 1)        = .5 * (4 * Ax(:, 2)  - 3 * Ax(:, 1)      - Ax(:, 3)) / dy          !&
         derY(:, ny)       = .5 * (3 * Ax(:, ny) - 4 * Ax(:, ny - 1) + Ax(:, ny - 2)) / dy     !&

         curlX = 0
         curlY = 0
         curlZ = derX - derY

      end select

      deallocate (Ax, Ay, derX, derY, vx, vy, vz)

   end subroutine

   subroutine verifyRing(field_grid, flag)
      implicit none
      type(field_grid_t), intent(in)                  :: field_grid
      integer(kind_particle), intent(in)              :: flag
      integer(kind_particle)                          :: nx, ny, ip, jp, rc, it
      real(kind_particle)                             :: dx, dy, L2x, L2y

!    select case (flag)
!          case (2)  !  ring on lies on yz, so it is a check on Ax
!
!
!          case (3)  !  ring on lies on xz, so it is a check on Ay
!
!
!          case default
!
!    end select

   end subroutine

end module
