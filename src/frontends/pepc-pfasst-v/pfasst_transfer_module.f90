! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

module pfasst_transfer_module


contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restrict(yF, yG, NvarF, NvarG, levelF)
    implicit none

    integer,      intent(in   ) :: NvarF, NvarG, levelF
    real(kind=8), intent(in   ) :: yF(NvarF)
    real(kind=8), intent(inout) :: yG(NvarG)

    yG = yF

  end subroutine restrict

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restrict_time_space_fas(t0, Dt, &
       ySDC_F, ySDC_G, fSDC_F, fSDC_G, NvarF, NvarG, levelF, tau_G)
    !  Restrict the values of ySDC_F to fit G and return the FAS
    !  correction.
    !
    !  The function values are restricted too (by re-evaluation).
    !  This does the time and space direction.

    use pfasst_calc_module, only: eval_f1, eval_f2
    use pfasst_helper_module, only: NnodesG, NnodesF, SmatF, SmatG, qnodesG

    implicit none

    integer,      intent(in ) :: NvarF, NvarG, levelF
    real(kind=8), intent(in ) :: t0, Dt
    real(kind=8), intent(in ) :: ySDC_F(NvarF,NnodesF)
    real(kind=8), intent(out) :: ySDC_G(NvarG,NnodesG)
    real(kind=8), intent(in ) :: fSDC_F(NvarF,NnodesF,2)
    real(kind=8), intent(out) :: fSDC_G(NvarG,NnodesG,2)
    real(kind=8), intent(out) :: tau_G(NvarG,NnodesG-1)

    ! Local Stuff
    integer :: m
    real(kind=8) :: tm(NnodesG)
    real(kind=8) :: Ic_of_G(NvarG,1:NnodesG-1)   ! Coarse integral of f_G
    real(kind=8) :: If_of_F(NvarF,1:NnodesF-1)   ! Fine integral of f_F
    real(kind=8) :: If_of_Fx(NvarG,1:NnodesF-1)  ! Fine integral of f_F coarsened in space
    real(kind=8) :: Ic_of_F(NvarG,1:NnodesG-1)   ! Fine integral of f_F coarsened in both
    real(kind=8) :: f1(NvarG)
    real(kind=8) :: f2(NvarG)
    real(kind=8) :: fsumSDC_F(NvarF,NnodesF)
    real(kind=8) :: fsumSDC_G(NvarG,NnodesG)

    ! Note: even if the number of variables and nodes is the same, we
    ! should still compute the FAS correction since the function
    ! evaluations may be different (eg, lower order operator for G)

    !  Restrict solution by simple sampling
    call restrict_time_space(ySDC_F, ySDC_G, NvarF, NvarG, NnodesF, NnodesG, levelF)

    tm = t0+Dt*qnodesG
    do m = 1,NnodesG
       call eval_f1(ySDC_G(:,m), tm(m), NvarG, levelF+1, f1)
       fSDC_G(:,m,1)=f1
       call eval_f2(ySDC_G(:,m), tm(m), NvarG, levelF+1, f2)
       fSDC_G(:,m,2)=f2
    end do

    fsumSDC_G = fSDC_G(:,:,1) + fSDC_G(:,:,2)
    fsumSDC_F = fSDC_F(:,:,1) + fSDC_F(:,:,2)

    !  Coarse integral of coarse function
    Ic_of_G = dt*matmul(fsumSDC_G, transpose(SmatG))
    !  Fine integral of fine function
    If_of_F = dt*matmul(fsumSDC_F, transpose(SmatF))

    call restrict_time_space(If_of_F, If_of_Fx, NvarF, NvarG, NnodesF-1, NnodesF-1, levelF)
    if (NnodesF == NnodesG) then
       Ic_of_F = If_of_Fx
    else
       !  Coarse value of fine Integral
       Ic_of_F = If_of_Fx(:,1:NnodesF-2:2) + If_of_Fx(:,2:NnodesF-1:2)
    end if

    !  FAS correction
    tau_G = -Ic_of_G + Ic_of_F

  end subroutine restrict_time_space_fas

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine restrict_time_space(ySDC_F, ySDC_G, NvarF, NvarG, NnodesF, NnodesG, levelF)
    !  Restrict the values of ySDC_F to fit G

    implicit none

    integer,      intent(in   )  :: NvarF, NvarG, NnodesF, NnodesG, levelF
    real(kind=8), intent(in   )  :: ySDC_F(NvarF,NnodesF)
    real(kind=8), intent(inout)  :: ySDC_G(NvarG,NnodesG)

    integer :: m, trat
    real(kind=8) :: yF(NvarF), yG(NvarG)

    if (NnodesG > 1) then
       trat = (NnodesF-1)/(NnodesG-1)
    else
       trat = 1
    end if

    do m = 1, NnodesG
       yF = ySDC_F(:,trat*(m-1)+1)
       call restrict(yF, yG, NvarF, NvarG, levelF)
       ySDC_G(:,m) = yG
    end do

  end subroutine restrict_time_space

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Interpolate correction of G and add to F.
  subroutine interpolate(yendF, yendG, NvarF, NvarG, levelF)
    implicit none
    integer,      intent(in   ) :: NvarF, NvarG, levelF
    real(kind=8), intent(inout) :: yendF(NvarF)
    real(kind=8), intent(in   ) :: yendG(NvarG)

    yendF = yendG
  end subroutine interpolate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate_time_space(t0, Dt, ySDC_F, ySDC_G, fSDC_F, fSDC_G, levelF)
    !  Interp correction of ySDC_G and add to ySDC_F in time and
    !  space.  Hence ySDC_F should match ySDC_G at coarse nodes.  Also
    !  interpolate fSDC_F (by re-evaluation).

    use pfasst_calc_module, only: eval_f1, eval_f2
    use pfasst_helper_module, only: NnodesG, NnodesF, NvarG, NvarF, qnodesF, qnodesG, InterpMat

    implicit none

    integer,      intent(in   ) :: levelF
    real(kind=8), intent(in   ) :: t0, Dt
    real(kind=8), intent(inout) :: ySDC_F(NvarF,NnodesF)
    real(kind=8), intent(in   ) :: ySDC_G(NvarG,NnodesG)
    real(kind=8), intent(inout) :: fSDC_F(NvarF,NnodesF,2)
    real(kind=8), intent(in   ) :: fSDC_G(NvarG,NnodesG,2)

    integer :: m, trat
    real(kind=8) :: tm(NnodesF)
    real(kind=8) :: ySDC_Fx(NvarG,NnodesG)  !  Restriction of ySDC_F
    real(kind=8) :: delG(NvarG,1:NnodesG)   !  Coarse increment of ySDC_F
    real(kind=8) :: delGF(NvarF,1:NnodesG)  !  Fine space Coarse time increment
    real(kind=8) :: f1(NvarF)               !  F space C time increment
    real(kind=8) :: f2(NvarF)               !  F space C time increment

    !  Make sure there is something to do
    if (NvarF == NvarG .and. NnodesF == NnodesG) then
       ySDC_F = ySDC_G
       fSDC_F = fSDC_G
       return
    endif

    trat = (NnodesF-1) / (NnodesG-1)

    !  We are going to interpolate the increments
    if (NvarF == NvarG) then
       delGF = ySDC_G - ySDC_F(:,::trat)
    else
       call restrict_time_space(ySDC_F, ySDC_Fx, NvarF, NvarG, NnodesF, NnodesG, levelF)
       delG = ySDC_G - ySDC_Fx

       delGF = 0.0
       call interpolate_space(delGF, delG, NvarF, NvarG, NnodesG, levelF)
    end if

    ySDC_F(:,::trat) = ySDC_F(:,::trat) + delGF

    if (NnodesF > NnodesG) then  !  Need to do time interp
       ySDC_F(:,2:NnodesF:trat) = ySDC_F(:,2:NnodesF:trat) &
            + matmul(delGF, transpose(InterpMat))
    end if

    !  Recompute fine function values
    tm = t0 + Dt*qnodesF
    do m = 1,NnodesF
       call eval_f1(ySDC_F(:,m), tm(m), NvarF, levelF, f1)
       call eval_f2(ySDC_F(:,m), tm(m), NvarF, levelF, f2)
       fSDC_F(:,m,1) = f1
       fSDC_F(:,m,2) = f2
    end do

  end subroutine interpolate_time_space

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate_space(ySDC_F, ySDC_G, NvarF, NvarG, Nnodes, levelF)
    !  Interp correction of ySDC_G and add to ySDC_F in space only.

    implicit none

    integer,      intent(in   ) :: NvarF, NvarG, Nnodes, levelF
    real(kind=8), intent(inout) :: ySDC_F(NvarF,Nnodes)
    real(kind=8), intent(in   ) :: ySDC_G(NvarG,Nnodes)

    real(kind=8) :: yF(NvarF)
    real(kind=8) :: yG(NvarG)

    integer :: i

    do i = 1, Nnodes
       yG = ySDC_G(:,i)
       yF = ySDC_F(:,i)
       call interpolate(yF, yG, NvarF, NvarG, levelF)
       ySDC_F(:,i) = yF
    end do

  end subroutine interpolate_space


end module pfasst_transfer_module
