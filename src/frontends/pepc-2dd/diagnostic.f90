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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!


module module_diagnostic

  use module_pepc_kinds
  use module_pepc_types
  use module_interaction_Specific_types
  use module_globals!, only: root,np,nt,my_rank, n_ranks
  implicit none

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)

  contains


      subroutine test_particles()

        use module_pepc_types
        use module_directsum
        use helper, only: get_time
        implicit none
        include 'mpif.h'

        integer(kind_particle), allocatable   :: tindx(:)
        real(kind_particle), allocatable      :: trnd(:)
        type(t_particle_results), allocatable :: trslt(:)
        integer(kind_particle)                :: tn, tn_global, ti
        integer                               :: rc
        real(kind_particle)                   :: Ex,ExTilde,Ey,EyTilde,E_norm_loc,E_norm_global,E_global,E_local
        real(kind_particle)                   :: Bz,BzTilde,B_norm,B_norm_global,B_local,B_global,B_norm_loc
        real(kind_particle)                   :: Ax,AxTilde,Ay,AyTilde,A_norm,A_norm_global,A_local,A_global,A_norm_loc
        real(kind_particle)                   :: phi,phiTilde,phi_global,phi_norm_loc,phi_norm_global,ta,tb,phi_local
        real(kind_particle)                   :: rho,rhoTilde,rho_norm_global,rho_global,rho_local,rho_norm_loc, Jx,JxTilde,Jy,JyTilde
        real(kind_particle)                   :: J_norm_global,J_global,J_local,J_norm_loc,devE,devPhi
        real(kind_particle)                   :: F_dar_loc,F_dar_glo,F_dar_den_loc,F_dar_den_glo,F_el_loc,F_el_glo,F_el_den_loc,F_el_den_glo,vx,vy,x,y,z,m,q


        ta = get_time()

    !    if(allocated(direct_L2)) then
    !      deallocate(direct_L2)
    !    end if
    !    if(allocated(direct_EL2A)) then
    !      deallocate(direct_EL2A)
    !    end if
    !    if(allocated(direct_EL2B)) then
    !      deallocate(direct_EL2B)
    !    end if
    !    if(allocated(direct_EL2Pot)) then
    !      deallocate(direct_EL2Pot)
    !    end if
    !    if(allocated(direct_EL2rho)) then
    !      deallocate(direct_EL2rho)
    !    end if
    !
    !    allocate(direct_L2(np))
    !    allocate(direct_EL2rho(np))
    !    allocate(direct_EL2Pot(np))
    !    allocate(direct_EL2B(np))
    !    allocate(direct_EL2A(np))
    !
    !    direct_L2 = -1.0_8
    !    direct_EL2rho = -1.0_8
    !    direct_EL2Pot = -1.0_8
    !    direct_EL2A = -1.0_8
    !    direct_EL2A = -1.0_8

        E_local       = 0.0_8
        E_norm_loc    = 0.0_8
        E_global      = 0.0_8
        E_norm_global = 0.0_8

        phi_local     = 0.0_8
        phi_norm_loc  = 0.0_8
        phi_global    = 0.0_8

        B_global      = 0.0_8
        B_norm_global = 0.0_8
        B_local       = 0.0_8
        B_norm_loc    = 0.0_8

        A_global      = 0.0_8
        A_norm_global = 0.0_8
        A_local       = 0.0_8
        A_norm_loc    = 0.0_8

        J_global      = 0.0_8
        J_norm_global = 0.0_8
        J_local       = 0.0_8
        J_norm_loc    = 0.0_8

        rho_global      = 0.0_8
        rho_norm_global = 0.0_8
        rho_local       = 0.0_8
        rho_norm_loc    = 0.0_8

        F_dar_loc           = 0.0_8
        F_dar_glo           = 0.0_8
        F_dar_den_loc       = 0.0_8
        F_dar_den_glo       = 0.0_8

        F_el_loc           = 0.0_8
        F_el_glo           = 0.0_8
        F_el_den_loc       = 0.0_8
        F_el_den_glo       = 0.0_8

        tn = np!particle_direct / n_ranks
        if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)

        allocate(tindx(tn), trnd(tn), trslt(tn))

        call random(trnd)

        tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1

        call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

        do ti = 1, tn

          vx          = particles(tindx(ti))%data%v(1)
          vy          = particles(tindx(ti))%data%v(2)

          x           = particles(tindx(ti))%x(1)
          y           = particles(tindx(ti))%x(2)

          m           = particles(tindx(ti))%data%m
          q           = particles(tindx(ti))%data%q

          Ex          = trslt(ti)%e(1)
          Ey          = trslt(ti)%e(2)

          ExTilde     = particles(tindx(ti))%results%e(1)
          EyTilde     = particles(tindx(ti))%results%e(2)

          Ax          = trslt(ti)%A(1)
          Ay          = trslt(ti)%A(2)

          AxTilde     = particles(tindx(ti))%results%A(1)
          AyTilde     = particles(tindx(ti))%results%A(2)

          Jx          = trslt(ti)%J(1)
          Jy          = trslt(ti)%J(2)

          JxTilde     = particles(tindx(ti))%results%J(1)
          JyTilde     = particles(tindx(ti))%results%J(2)

          Bz          = trslt(ti)%B(3)

          BzTilde     = particles(tindx(ti))%results%B(3)

          rho         = trslt(ti)%rho

          rhoTilde    = particles(tindx(ti))%results%rho

          E_local     = E_local + ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2
          E_norm_loc  = E_norm_loc + Ex**2 + Ey**2

          phi         = trslt(ti)%pot
          phiTilde    = particles(tindx(ti))%results%pot

          phi_local        = phi_local +  ( phi - phiTilde )**2
          phi_norm_loc     = phi_norm_loc +  phi**2


          A_local     = A_local + ( AxTilde -  Ax )**2 + ( AyTilde - Ay )**2
          A_norm_loc  = A_norm_loc + Ax**2 + Ay**2

          B_local     = B_local + ( BzTilde -  Bz )**2
          B_norm_loc  = B_norm_loc + Bz**2

          J_local     = J_local + ( JxTilde -  Jx )**2 + ( JyTilde - Jy )**2
          J_norm_loc  = J_norm_loc + Jx**2 + Jy**2

          rho_local        = rho_local +  ( rho - rhoTilde )**2
          rho_norm_loc     = rho_norm_loc +  rho**2

          F_dar_loc       = F_dar_loc + (q)**2*( ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2 &
                        +   ( vx**2 + vy**2 )*( BzTilde -  Bz )**2 + 2*( BzTilde -  Bz )*( vx*( ExTilde -  Ex ) + vy* ( EyTilde -  Ey ) )  )

          F_dar_den_loc  = F_dar_den_loc + (q)**2*( Ex**2 +  Ey**2 &
                        +   ( vx**2 + vy**2 )*Bz**2 + 2*Bz*( vx*Ex+ vy*Ey  )  )

          F_el_loc       = F_el_loc + (q)**2*( ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2 )

          F_el_den_loc  = F_el_den_loc + (q)**2*( Ex**2 +  Ey**2 )



    !      direct_L2(tindx(ti)) = L2

    !      EL2A          = &
    !                    (particles(tindx(ti))%results%A(1) - trslt(ti)%A(1))**2+ &
    !                    (particles(tindx(ti))%results%A(2) - trslt(ti)%A(2))**2
    !      EL2Asum_local = EL2Asum_local + EL2A
    !      direct_EL2A(tindx(ti)) = EL2A

    !      EL2B          = (particles(tindx(ti))%results%B(3) - trslt(ti)%B(3))**2
    !      EL2Bsum_local = EL2Bsum_local + EL2B
    !      direct_EL2B(tindx(ti)) = EL2B

    !      direct_EL2Pot(tindx(ti)) = EL2Pot

    !      EL2rho      = &
    !                    (particles(tindx(ti))%results%rho - trslt(ti)%rho)**2
    !      EL2rhosum_local = EL2rhosum_local + EL2rho
    !      direct_EL2rho(tindx(ti)) = EL2rho
        end do

        call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_local, phi_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_norm_loc, phi_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_local, E_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_norm_loc, E_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_norm_loc, B_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_local, B_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_local, A_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_norm_loc, A_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_local, J_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_norm_loc, J_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(rho_local, rho_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(rho_norm_loc, rho_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_loc, F_dar_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_den_loc, F_dar_den_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_loc, F_el_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_den_loc, F_el_den_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)


        devE                 = sqrt(E_global)/(tn_global-1.0_8)
        devPhi               = sqrt(phi_global)/(tn_global-1.0_8)

        phi_global           = sqrt(phi_global) / sqrt(phi_norm_global)
        rho_global           = sqrt(rho_global) / sqrt(rho_norm_global)
        E_global             = sqrt(E_global) / sqrt(E_norm_global)
        A_global             = sqrt(A_global) / sqrt(A_norm_global)
        B_global             = sqrt(B_global) / sqrt(B_norm_global)
        J_global             = sqrt(J_global) / sqrt(J_norm_global)
        F_dar_glo            = sqrt(F_dar_glo  / F_dar_den_glo )
        F_el_glo             = sqrt(F_el_glo  / F_el_den_glo )


        tb = get_time()
        if(root) then
          write(*,'(a,i12)')    " == [direct test] number tested particles         : ", tn
    !      write(*,'(a,es12.4)') " == [direct test] l2 A                            : ", A_norm_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 A error                      : ", A_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 B                            : ", B_norm_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 B error                      : ", B_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in El Pot        : ", phi_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in E             : ", E_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in A             : ", A_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in B             : ", B_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in J             : ", J_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in rho           : ", rho_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in F dar         : ", F_dar_glo
          write(*,'(a,es12.4)') " == [direct test] Relative error in F el          : ", F_el_glo
          write(*,'(a,es12.4)') " == [direct test] L2 error in E                   : ", devE
          write(*,'(a,es12.4)') " == [direct test] L2 error in El Pot              : ", devPhi
          write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", tb - ta


        end if

        deallocate(tindx)
        deallocate(trnd)
        deallocate(trslt)

      end subroutine test_particles
!
!
!
      subroutine write_field(p)
        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle) :: indexvec(7), tmp
        character (len=12), dimension(7) :: stringvec

        integer(kind_particle) :: ifield, ip
        indexvec = (/ 1, 2, 3, 4, 5, 6, 7/)
        stringvec = (/ "data/xyz.dat", "data/phi.dat", "data/EEE.dat","data/BBB.dat", "data/rho.dat","data/JJJ.dat","data/AAA.dat"/)

        do ifield = 1,7

            open (unit=indexvec(ifield),file=stringvec(ifield),action="write",status="replace")

        enddo


        write(*,*) " ==  fields written "
        do tmp = 1, np

                do ip = 1, np

                    if ( p(ip)%label .eq. tmp ) then

                        write (1,*) p(ip)%x(1),p(ip)%x(2),0.0_8
                        write (2,*) p(ip)%results%pot
                        write (3,*) p(ip)%results%E(1),p(ip)%results%E(2),0.0_8
                        write (4,*) 0.0_8,0.0_8,p(ip)%results%B(3)
                        write (5,*) p(ip)%results%rho
                        write (6,*) p(ip)%results%J(1),p(ip)%results%J(2),0.0_8
                        write (7,*) p(ip)%results%A(1),p(ip)%results%A(2),0.0_8

                    endif

                enddo

        enddo


        do ifield = 1,7

            close ( unit=indexvec(ifield) )

        enddo

      end subroutine

      subroutine write_field_rectilinear_grid(p)
        use module_vtk
        use helper, only: get_time
        implicit none
        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle)  :: i
        integer, dimension(2,3) ::  nnp
        type(vtkfile_rectilinear_grid) :: vtk
        integer :: vtk_step
        real(kind_particle) :: time
        real(kind_particle) :: ta, tb,Ex(np),Ey(np),Ez(np),x(np),y(np),z(np)


        ta = get_time()
        time = dt * step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        x  = p(:)%x(1)
        y  = p(:)%x(2)
        z  = p(:)%x(3)

        Ex = p(:)%results%E(1)
        Ey = p(:)%results%E(2)
        Ez = p(:)%results%E(3)
    !    nnp(1) = np
    !    nnp(2) = np
    !    nnp(3) = 1


        !nnp = np

        call vtk%create_parallel("fields_on_grid", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(nnp, nnp)
        call vtk%startcoordinates()
        call vtk%write_data_array("x_coordinate", x )
        call vtk%write_data_array("y_coordinate", y )
        call vtk%write_data_array("z_coordinate", 0 )

        call vtk%finishcoordinates()
        call vtk%startpointdata()

        call vtk%write_data_array("E", Ex,  Ey,  Ez )
              ! no point data here
        call vtk%finishpointdata()
        call vtk%startcelldata()
              ! no cell data here
        call vtk%finishcelldata()
        call vtk%write_final()
        call vtk%close()



        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta



           !call vtk%write_data_array(vectorname, vectorvalues(:,:,:,1), vectorvalues(:,:,:,2), vectorvalues(:,:,:,3))

      end subroutine
!
      subroutine write_particles(p)
        use module_vtk
        use helper, only: get_time
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer(kind_particle) :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time,vect(np)
        real*8 :: ta, tb

        ta = get_time()
        time = dt * step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        vect = 0.0_8

        write(*,*) "particles"
        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0_kind_particle)
        call vtk%startpoints()
        call vtk%write_data_array("x", p(1:np)%x(1), p(1:np)%x(2), vect(1:np) )
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("v", p(1:np)%data%v(1), &
                                       p(1:np)%data%v(2), &
                                       vect(1:np) )

        call vtk%write_data_array("E", p(1:np)%results%e(1), &
                                       p(1:np)%results%e(2), &
                                       vect(1:np) )

        call vtk%write_data_array("A", p(1:np)%results%A(1), &
                                       p(1:np)%results%A(2), &
                                       vect )

        call vtk%write_data_array("rho", p(1:np)%results%rho )

        call vtk%write_data_array("B",   vect(1:np),&
                                         vect(1:np),&
                                         p(1:np)%results%B(3) )

        call vtk%write_data_array("J",   p(1:np)%results%J(1),&
                                         p(1:np)%results%J(2),&
                                         vect(1:np) )

        call vtk%write_data_array("phi", p(1:np)%results%pot)
        call vtk%write_data_array("q", p(1:np)%data%q)
        call vtk%write_data_array("m", p(1:np)%data%m)
        call vtk%write_data_array("work", p(1:np)%work)
        call vtk%write_data_array("pelabel", p(1:np)%label)
        call vtk%write_data_array("local index", [(i,i=1,np)])
        call vtk%write_data_array("processor", int(np, kind = 4), my_rank)
!        if(particle_test) call vtk%write_data_array("L2 error", direct_L2(1:np))
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

      end subroutine write_particles
!
      subroutine write_domain(p)

        use module_vtk
        use module_vtk_helpers
        use module_pepc, only: global_tree
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer :: vtk_step

        ! output of tree diagnostics
        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif
        call vtk_write_branches(step,  dt * step, vtk_step, global_tree)
        call vtk_write_spacecurve(step, dt * step, vtk_step, p)

      end subroutine write_domain



end module
