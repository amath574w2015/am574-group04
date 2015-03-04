subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
    ! Need to import gamma1, area
    ! Called to update q by solving source term equation 
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)
    ! added pressure storage
    dimension   press(1,1-mbc:maxmx+mbc)
    dimension   pressstar(1,1-mbc:maxmx+mbc)
    dimension   Qstarstar(1,1-mbc:maxmx+mbc)
    dimension   area(1,1-mbc:maxmx+mbc)
    common /cparam/  gamma

! Define area
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
	 area = 1.398 + 0.347 * tanh(0.8*xcell - 0.4)
  150    continue

! Calculate P from q. P=gamma1*(e + 0.5 rho*u^2)
press=(gamma-1.0)*(q(3,:) + 0.5* (q(2,:)**2)/(q(1,:)))

! Solve the ODE: (rho*u*A)_t=P dA/dx using Runge-Kutta
do 20 i=2-mbc,mx+mbc
	Qstarstar(i)=q(2,i)+ (dt/2) * press(i) * (area(i+1)-area(i))/dx
	pressstar=(gamma-1.0)*(Qstarstar(3,:) + 0.5* (Qstarstar(2,:)**2)/(Qstarstar(1,:)))
	q(2,i)=q(2,i)+ dt * pressstar(i) * (area(i+1)-area(i))/dx
20 end do

end subroutine src1
