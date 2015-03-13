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
    real(kind=8)   press,pressstarstar
    real(kind=8)   psistar(3)
    real(kind=8)   psistarstar(3)
    real(kind=8)   Qstarstar(3)
    real(kind=8)   darea(1-mbc:mx+mbc)
    real(kind=8)   area(1-mbc:mx+mbc)
    real(kind=8)   radius(1-mbc:mx+mbc)
    real(kind=8)   pi,gamma,xcell
    common /cparam/  gamma

    pi=4d0*atan(1d0)

!   Define area
    area(:) = 0.d0
    darea(:) = 0.d0
    call def_area(dx,mx,xlower,mbc,area)
    call def_darea(dx,mx,xlower,mbc,darea)

!   Calculate P from q. P=gamma1*(e + 0.5 rho*u^2)

!   Solve the coupled ODEs, Eq 17.40 Levque
    do 20 i=2-mbc,mx+mbc
       press = (gamma-1.d0)*(q(3,i) - 0.5d0* q(2,i)**2.d0/(q(1,i)))

!      Calculate psi* based on the Q*
       psistar(1)= 0.0d0 !-1.d0/area(i) * darea(i) * (q(2,i))
       psistar(2)= press*darea(i) !  -1.d0/area(i) * darea(i) * (q(2,i)**2/q(1,i))
       psistar(3)= 0.0d0 !-1.d0/area(i) * darea(i) * (q(2,i)/q(1,i))*(q(3,i)+press)

!      Calculate Q** based on the psi*
       Qstarstar(1)=q(1,i)+ (dt/2) * psistar(1)
       Qstarstar(2)=q(2,i)+ (dt/2) * psistar(2)
       Qstarstar(3)=q(2,i)+ (dt/2) * psistar(3)

       pressstarstar = (gamma-1.d0)*(Qstarstar(3)-0.5d0*Qstarstar(2)**2.d0 &
              /Qstarstar(1))

!      Calculate psi** based on Q**
       psistarstar(1)= 0.d0 !-1.d0/area(i) * darea(i) * (Qstarstar(2))
       psistarstar(2)= pressstarstar*darea(i) !-1.d0/area(i) * darea(i) * ( &
!           Qstarstar(2)**2/Qstarstar(1))
       psistarstar(3)= 0.0d0 ! -1.d0/area(i) * darea(i) &
!        * (Qstarstar(2)/Qstarstar(1))*(Qstarstar(3)+pressstarstar)

!      Calculate the updated Q based on psi**
       q(1,i)=q(1,i)+ dt * psistarstar(1) 
       q(2,i)=q(2,i)+ dt * psistarstar(2) 
       q(3,i)=q(3,i)+ dt * psistarstar(3) 
    20 end do

end subroutine src1