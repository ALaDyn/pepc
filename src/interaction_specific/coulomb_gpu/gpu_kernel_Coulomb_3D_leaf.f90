!!   case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
subroutine kernel1(particle, eps2, WORKLOAD_PENALTY_INTERACTION)
   use module_pepc_types
   use module_interaction_specific_types
   use module_interaction_specific, only: MAX_IACT_PARTNERS
#ifdef __OPENACC
   use openacc
#endif
   implicit none

   type(t_particle_thread), intent(inout) :: particle
   real*8, intent(in) :: eps2, WORKLOAD_PENALTY_INTERACTION

   integer :: idx
   real*8 :: exyz(3), phic, exyz1, exyz2, exyz3
   real*8 :: delta(3), dist2, charge
   real*8 :: rd,rd3charge
   real*8 :: e_1, e_2, e_3, pot

   !real*8 :: delta_(1:3,1:MAX_IACT_PARTNERS)
   !real*8 :: charge_(1:MAX_IACT_PARTNERS)
   type chargedelta
      real*8 :: delta1(1:MAX_IACT_PARTNERS)
      real*8 :: delta2(1:MAX_IACT_PARTNERS)
      real*8 :: delta3(1:MAX_IACT_PARTNERS)
      real*8 :: charge(1:MAX_IACT_PARTNERS)
   end type chargedelta
   type(chargedelta) :: gpu

   !!write(*,*) 'kernel will start with ',particle%queued_l,' iterations'

   !!!call acc_init(acc_device_nvidia)
   !!write(*,*) 'num devices ', acc_get_num_devices(acc_device_nvidia)

   do idx = 1, particle%queued_l
!      charge_(idx)    = particle%partner_l(idx)%charge
!      delta_(1:3,idx) = particle%partner_l(idx)%delta(1:3)
       gpu%delta1(idx) = particle%partner_l(idx)%delta(1)
       gpu%delta2(idx) = particle%partner_l(idx)%delta(2)
       gpu%delta3(idx) = particle%partner_l(idx)%delta(3)
       gpu%charge(idx) = particle%partner_l(idx)%charge
   enddo

   !!!!!$acc kernels  !! with kernels, PGI will try and make loops GPU-able. could be at start of code...
   !!! code...
   !!!!!$acc end kernels

   !!call acc_set_device_num(0, acc_device_nvidia)

#ifdef __OPENACC
   !!!!$acc parallel loop reduction(+: e_1, e_2, e_3, pot) copyin(delta_(1:3,1:particle%queued_l), charge_(1:particle%queued_l)) private(dist2, rd, rd3charge, phic, exyz)
   !$acc parallel loop reduction(+: e_1, e_2, e_3, pot) copyin(gpu) private(dist2, rd, rd3charge, phic, exyz1, exyz2, exyz3)
   !! exyz(3) vs. exyz[123] == 109secs vs. 78secs!
#endif
   do idx = 1, particle%queued_l

      dist2     =         gpu%delta1(idx) * gpu%delta1(idx)
      dist2     = dist2 + gpu%delta2(idx) * gpu%delta2(idx)
      dist2     = dist2 + gpu%delta3(idx) * gpu%delta3(idx)
      !dist2     =         delta_(1,idx) * delta_(1,idx)
      !dist2     = dist2 + delta_(2,idx) * delta_(2,idx)
      !dist2     = dist2 + delta_(3,idx) * delta_(3,idx)
      !dist2     = dot_product(delta_(1:3,idx),delta_(1:3,idx))
      !dist2     = sum(delta_(:,idx) * delta_(:,idx))

      rd        = 1/sqrt(dist2+eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd3charge = gpu%charge(idx)*rd*rd*rd
      !rd3charge = charge_(idx)*rd*rd*rd

      phic  = gpu%charge(idx)*rd
      exyz1 = rd3charge*gpu%delta1(idx)
      exyz2 = rd3charge*gpu%delta2(idx)
      exyz3 = rd3charge*gpu%delta3(idx)
      !exyz = rd3charge*[gpu%delta1(idx), gpu%delta2(idx), gpu%delta3(idx)]
      !phic = charge_(idx)*rd
      !exyz = rd3charge*delta_(1:3,idx)

      e_1  = e_1  + exyz1
      e_2  = e_2  + exyz2
      e_3  = e_3  + exyz3
      !e_1  = e_1  + exyz(1)
      !e_2  = e_2  + exyz(2)
      !e_3  = e_3  + exyz(3)
      pot  = pot  + phic

   end do
#ifdef __OPENACC
   !$acc end parallel loop
#endif

   particle%results%e   = particle%results%e + [e_1, e_2, e_3]
   particle%results%pot = particle%results%pot + pot
   particle%work        = particle%work + particle%queued_l * WORKLOAD_PENALTY_INTERACTION

   !!!call acc_shutdown(acc_device_nvidia)

   return

end subroutine kernel1
