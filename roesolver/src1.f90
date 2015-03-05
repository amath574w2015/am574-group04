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
    real(kind=8)   psistar(1-mbc:mx+mbc)
    real(kind=8)   psistarstar(1-mbc:mx+mbc)
    real(kind=8)   Qstarstar(1-mbc:mx+mbc)
    real(kind=8)   darea(1-mbc:mx+mbc)
    real(kind=8)   area(1-mbc:mx+mbc)
    common /cparam/  gambn ma

!   Define area
    do 150 i=1,mx
       xcell = xlower + (i-0.5d0)*dx
       area(i) = 1.398 + 0.347 * tanh(0.8*xcell - 0.4)
       darea(i) = 0.2776d0*(1.0d0 - tanh(0.8d0*xcell - 4.0d0)**2.0d0)
    150    continue

!   Calculate P from q. P=gamma1*(e + 0.5 rho*u^2)
    do 11 i = 1,mx+mbc
       press(i) =(gamma-1.d0)*(q(3,i) + 0.5d0* (q(2,i)**2.d0)/(q(1,i)))
    11 end do

!   Solve the coupled ODEs, Eq 17.40 Levque
    do 20 i=2-mbc,mx+mbc
!      Calculate psi* based on the Q*
       psistar(1,i)= -1.d0/area(i) * darea(i) * (q(2,i))
       psistar(2,i)= -1.d0/area(i) * darea(i) * (q(2,i)**2/q(1,i))
       psistar(3,i)= -1.d0/area(i) * darea(i) * (q(2,i)/q(1,i))*(q(3,i)+press(i))

!      Calculate Q** based on the psi*
       Qstarstar(1,i)=q(1,i)+ (dt/2) * psistar(1,i)
       Qstarstar(2,i)=q(2,i)+ (dt/2) * psistar(2,i)
       Qstarstar(3,i)=q(2,i)+ (dt/2) * psistar(3,i)

!      Calculate psi** based on Q**
       psistarstar(1,i)=-1.d0/area(i) * darea(i) * (Qstarstar(2,i))
       psistarstar(2,i)=-1.d0/area(i) * darea(i) * (Qstarstar(2,i)**2/Qstarstar(1,i))
       psistarstar(3,i)=-1.d0/area(i) * darea(i) &
        * (Qstarstar(2,i)/Qstarstar(1,i))*(Qstarstar(3,i)+press(i))

!      Calculate the updated Q based on psi**
       q(1,i)=q(1,i)+ dt * psistarstar(1,i) * darea(i)
       q(2,i)=q(2,i)+ dt * psistarstar(2,i) * darea(i)
       q(3,i)=q(3,i)+ dt * psistarstar(3,i) * darea(i)
    20 end do

end subroutine src1
