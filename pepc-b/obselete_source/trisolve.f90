! -------------------------------------------------
!
!   Triadiagonal matrix solver.
!
!   Solves equation system:
!
!        alpha_i x_(i-1) + beta_i x_i + gamma_i x_(i+1) = y_i
!
! ------------------------------------------------

subroutine trisolve(alpha,beta,gamma,y,x,n,nmax)

 implicit none
 integer, intent(in) :: n,nmax
 complex, dimension(nmax) :: alpha, beta, gamma,y,x, q
 integer :: i

  q(1) = beta(1)
  x(1) = y(1)/q(1)

  !  forward elimination

  do i = 2,n
     q(i) = beta(i) - alpha(i)*gamma(i-1)/q(i-1)
     x(i) = (y(i) - alpha(i)*x(i-1))/q(i)
  end do

  !  back substitution

  do i=n-1,1,-1
     x(i) = x(i) - gamma(i)/q(i)*x(i+1)
  end do
end subroutine trisolve

