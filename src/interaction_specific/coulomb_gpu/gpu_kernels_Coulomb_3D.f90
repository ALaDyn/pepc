!!   case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
subroutine kernel_leaf(particle, eps2, WORKLOAD_PENALTY_INTERACTION)
   use module_pepc_types
   use module_interaction_specific_types
#ifdef __OPENACC
   use openacc
#endif
   implicit none

   type(t_particle_thread), intent(inout) :: particle
   real*8, intent(in) :: eps2, WORKLOAD_PENALTY_INTERACTION

   integer :: idx
   real*8 :: dist2
   real*8 :: rd,rd3charge
   real*8 :: e_1, e_2, e_3, pot

   !!write(*,*) 'kernel will start with ',particle%queued_l,' iterations'

   do idx = 1, particle%queued_l
       gpu_l%delta1(idx) = particle%partner_l(idx)%delta(1)
       gpu_l%delta2(idx) = particle%partner_l(idx)%delta(2)
       gpu_l%delta3(idx) = particle%partner_l(idx)%delta(3)
       gpu_l%charge(idx) = particle%partner_l(idx)%charge
   enddo

   e_1 = 0.d0
   e_2 = 0.d0
   e_3 = 0.d0
   pot = 0.d0

#ifdef __OPENACC
   !$acc parallel loop reduction(+: e_1, e_2, e_3, pot) present_or_copyin(gpu_l) private(dist2, rd, rd3charge)
#endif
   do idx = 1, particle%queued_l

      dist2     =         gpu_l%delta1(idx) * gpu_l%delta1(idx)
      dist2     = dist2 + gpu_l%delta2(idx) * gpu_l%delta2(idx)
      dist2     = dist2 + gpu_l%delta3(idx) * gpu_l%delta3(idx)
      !dist2     =         delta_(1,idx) * delta_(1,idx)
      !dist2     = dist2 + delta_(2,idx) * delta_(2,idx)
      !dist2     = dist2 + delta_(3,idx) * delta_(3,idx)
      !dist2     = dot_product(delta_(1:3,idx),delta_(1:3,idx))
      !dist2     = sum(delta_(:,idx) * delta_(:,idx))

      rd        = 1/sqrt(dist2+eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd3charge = gpu_l%charge(idx)*rd*rd*rd

      pot = pot + gpu_l%charge(idx)*rd
      e_1 = e_1 + rd3charge*gpu_l%delta1(idx)
      e_2 = e_2 + rd3charge*gpu_l%delta2(idx)
      e_3 = e_3 + rd3charge*gpu_l%delta3(idx)

   end do
#ifdef __OPENACC
   !$acc end parallel loop
#endif

   particle%results%e   = particle%results%e + [e_1, e_2, e_3]
   particle%results%pot = particle%results%pot + pot
   particle%work        = particle%work + particle%queued_l * WORKLOAD_PENALTY_INTERACTION

   return

end subroutine kernel_leaf

subroutine kernel_node(particle, eps2, WORKLOAD_PENALTY_INTERACTION)
   use module_pepc_types
   use module_interaction_specific_types
#ifdef __OPENACC
   use openacc
#endif
   implicit none

   type(t_particle_thread), intent(inout) :: particle
   real*8, intent(in) :: eps2, WORKLOAD_PENALTY_INTERACTION

   integer :: idx
   real*8 :: dist2
   real*8 :: e_1, e_2, e_3, pot

   real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6

   do idx = 1, particle%queued
       gpu%delta1(idx) = particle%partner(idx)%delta(1)
       gpu%delta2(idx) = particle%partner(idx)%delta(2)
       gpu%delta3(idx) = particle%partner(idx)%delta(3)
       gpu%charge(idx) = particle%partner(idx)%node%charge
       gpu%dip1(idx)   = particle%partner(idx)%node%dip(1)
       gpu%dip2(idx)   = particle%partner(idx)%node%dip(2)
       gpu%dip3(idx)   = particle%partner(idx)%node%dip(3)
       gpu%quad1(idx)  = particle%partner(idx)%node%quad(1)
       gpu%quad2(idx)  = particle%partner(idx)%node%quad(2)
       gpu%quad3(idx)  = particle%partner(idx)%node%quad(3)
       gpu%xyquad(idx) = particle%partner(idx)%node%xyquad
       gpu%yzquad(idx) = particle%partner(idx)%node%yzquad
       gpu%zxquad(idx) = particle%partner(idx)%node%zxquad
   enddo

   e_1 = 0.d0
   e_2 = 0.d0
   e_3 = 0.d0
   pot = 0.d0

#ifdef __OPENACC
   !$acc update device(gpu)
   !$acc parallel loop reduction(+: e_1, e_2, e_3, pot) present(gpu) private(dist2,rd,dx,dy,dz,r,dx2,dy2,dz2,dx3,dy3,dz3,rd2,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6)
#endif
   do idx = 1, particle%queued

      dist2     =         gpu%delta1(idx) * gpu%delta1(idx)
      dist2     = dist2 + gpu%delta2(idx) * gpu%delta2(idx)
      dist2     = dist2 + gpu%delta3(idx) * gpu%delta3(idx)

      dx = gpu%delta1(idx)
      dy = gpu%delta2(idx)
      dz = gpu%delta3(idx)

      r  = sqrt(dist2+eps2) ! eps2 is added in calling routine to have plummer instead of coulomb here
      rd = 1.d0/r
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

      fd1 = 3.d0*dx2*rd5 - rd3
      fd2 = 3.d0*dy2*rd5 - rd3
      fd3 = 3.d0*dz2*rd5 - rd3
      fd4 = 3.d0*dx*dy*rd5
      fd5 = 3.d0*dy*dz*rd5
      fd6 = 3.d0*dx*dz*rd5

      pot = pot + gpu%charge(idx)*rd                                                 &  !  monopole term
            + (dx*gpu%dip1(idx) + dy*gpu%dip2(idx) + dz*gpu%dip3(idx))*rd3           &  !  dipole
            + 0.5d0*(fd1*gpu%quad1(idx)  + fd2*gpu%quad2(idx)  + fd3*gpu%quad3(idx)) &  !  quadrupole
            +        fd4*gpu%xyquad(idx) + fd5*gpu%yzquad(idx) + fd6*gpu%zxquad(idx)

      e_1 = e_1 + gpu%charge(idx)*dx*rd3                                           &  ! monopole term
                + fd1*gpu%dip1(idx) + fd4*gpu%dip2(idx) + fd6*gpu%dip3(idx)        &  ! dipole term
                + 3.d0   * (                                                       &  ! quadrupole term
                   0.5d0 * (                                                       &
                       ( 5.d0*dx3   *rd7 - 3.d0*dx*rd5 )*gpu%quad1(idx)            &
                     + ( 5.d0*dx*dy2*rd7 -      dx*rd5 )*gpu%quad2(idx)            &
                     + ( 5.d0*dx*dz2*rd7 -      dx*rd5 )*gpu%quad3(idx)            &
                   )                                                               &
                   + ( 5.d0*dy*dx2  *rd7 - dy*rd5 )*gpu%xyquad(idx)                &
                   + ( 5.d0*dz*dx2  *rd7 - dz*rd5 )*gpu%zxquad(idx)                &
                   + ( 5.d0*dx*dy*dz*rd7          )*gpu%yzquad(idx)                &
                  )

      e_2 = e_2 + gpu%charge(idx)*dy*rd3                                           &
                + fd2*gpu%dip2(idx) + fd4*gpu%dip1(idx) + fd5*gpu%dip3(idx)        &
                + 3 * (                                                            &
                   0.5d0 * (                                                       &
                       ( 5*dy3*rd7    - 3*dy*rd5 )*gpu%quad2(idx)                  &
                     + ( 5*dy*dx2*rd7 -   dy*rd5 )*gpu%quad1(idx)                  &
                     + ( 5*dy*dz2*rd7 -   dy*rd5 )*gpu%quad3(idx)                  &
                   )                                                               &
                   + ( 5*dx*dy2  *rd7 - dx*rd5 )*gpu%xyquad(idx)                   &
                   + ( 5*dz*dy2  *rd7 - dz*rd5 )*gpu%yzquad(idx)                   &
                   + ( 5*dx*dy*dz*rd7          )*gpu%zxquad(idx)                   &
                  )

      e_3 = e_3 + gpu%charge(idx)*dz*rd3                                           &
                + fd3*gpu%dip3(idx) + fd5*gpu%dip2(idx) + fd6*gpu%dip1(idx)        &
                + 3 * (                                                            &
                   0.5d0 * (                                                       &
                     + ( 5*dz3   *rd7 - 3*dz*rd5 )*gpu%quad3(idx)                  &
                     + ( 5*dz*dy2*rd7 -   dz*rd5 )*gpu%quad2(idx)                  &
                     + ( 5*dz*dx2*rd7 -   dz*rd5 )*gpu%quad1(idx)                  &
                                  )                                                &
                   + ( 5*dx*dz2  *rd7 - dx*rd5 )*gpu%zxquad(idx)                   &
                   + ( 5*dy*dz2  *rd7 - dy*rd5 )*gpu%yzquad(idx)                   &
                   + ( 5*dx*dy*dz*rd7          )*gpu%xyquad(idx)                   &
                  )
   end do
#ifdef __OPENACC
   !$acc end parallel loop
#endif

   particle%results%e   = particle%results%e + [e_1, e_2, e_3]
   particle%results%pot = particle%results%pot + pot
   particle%work        = particle%work + particle%queued * WORKLOAD_PENALTY_INTERACTION

   return

end subroutine kernel_node
