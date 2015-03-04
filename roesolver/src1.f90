subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
    ! Need to import gamma1, area
    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    !  
 
    implicit double precision (a-h,o-z)
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)
    ! added pressure storage
    real(kind=8)   press(1-mbc:mx+mbc)
    real(kind=8)   pressstar(1-mbc:mx+mbc)
    real(kind=8)   Qstarstar(1-mbc:mx+mbc)
    real(kind=8)   darea
    common /cparam/  gamma

!   Define area
    do 150 i=1,mx
       xcell = xlower + (i-0.5d0)*dx
       darea = 0.2776d0*(1.0d0 - tanh(0.8d0*xcell - 4.0d0)**2.0d0)
    150    continue

!   Calculate P from q. P=gamma1*(e + 0.5 rho*u^2)
    do 11 i = 1,mx+mbc
       press(i) =(gamma-1.d0)*(q(3,i) + 0.5d0* (q(2,i)**2.d0)/(q(1,i)))
    11 end do

!   Solve the ODE: (rho*u*A)_t=P dA/dx using Runge-Kutta
    do 20 i=2-mbc,mx+mbc
!      Qstarstar is the 2nd conserved variable at the mid point, rho u A
       Qstarstar(i)=q(2,i)+ (dt/2) * press(i) * darea
       pressstar=(gamma-1.d0)*(q(3,:) + 0.5d0*(Qstarstar(:)**2.d0)/(q(1,:)))
       q(2,i)=q(2,i)+ dt * pressstar(i) * darea
    20 end do

end subroutine src1
