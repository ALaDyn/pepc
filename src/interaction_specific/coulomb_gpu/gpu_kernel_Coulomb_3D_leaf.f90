!!   case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
subroutine kernel1(particle, eps2, WORKLOAD_PENALTY_INTERACTION)
   use module_pepc_types
   use module_interaction_specific_types
   use module_interaction_specific, only: MAX_IACT_PARTNERS
   implicit none

   type(t_particle_thread), intent(inout) :: particle
   real*8, intent(in) :: eps2, WORKLOAD_PENALTY_INTERACTION

   integer :: idx
   real*8 :: exyz(3), phic
   real*8 :: delta(3), dist2, charge
   real*8 :: rd,rd3charge
   real*8 :: e_1, e_2, e_3, pot

   real*8 :: delta_(1:3,1:MAX_IACT_PARTNERS)
   real*8 :: charge_(1:MAX_IACT_PARTNERS)

   !!write(*,*) 'kernel will start with ',particle%queued_l,' iterations'

   do idx = 1, particle%queued_l
      charge_(idx)    = particle%partner_l(idx)%charge
      delta_(1:3,idx) = particle%partner_l(idx)%delta(1:3)
   enddo

   !!!!!$acc kernels  !! with kernels, PGI will try and make loops GPU-able. could be at start of code...
   !!! code...
   !!!!!$acc end kernels

   !$acc parallel loop reduction(+: e_1, e_2, e_3, pot) copyin(delta_(1:3,1:particle%queued_l), charge_(1:particle%queued_l)) private(rd, rd3charge, phic, exyz)
   do idx = 1, particle%queued_l

      dist2     =         delta_(1,idx) * delta_(1,idx)
      dist2     = dist2 + delta_(2,idx) * delta_(2,idx)
      dist2     = dist2 + delta_(3,idx) * delta_(3,idx)
      !dist2     = dot_product(delta_(1:3,idx),delta_(1:3,idx))
      !dist2     = sum(delta_(:,idx) * delta_(:,idx))

      rd        = 1/sqrt(dist2+eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd3charge = charge_(idx)*rd*rd*rd

      phic = charge_(idx)*rd
      exyz = rd3charge*delta_(1:3,idx)

      e_1  = e_1  + exyz(1)
      e_2  = e_2  + exyz(2)
      e_3  = e_3  + exyz(3)
      pot  = pot  + phic

   end do
   !$acc end parallel loop

   particle%results%e   = particle%results%e + [e_1, e_2, e_3]
   particle%results%pot = particle%results%pot + pot
   particle%work        = particle%work + particle%queued_l * WORKLOAD_PENALTY_INTERACTION

   return

end subroutine kernel1
