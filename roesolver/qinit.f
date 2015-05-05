c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c

      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      real(kind=8) :: at,gamma,gamma_m,gamma_a
      real(kind=8) :: area(1-mbc:mx+mbc)
      real(kind=8) :: c(1-mbc:mx+mbc)
      real(kind=8) :: radius,rho,u
      real(kind=8) :: pi,patm,xcell,p
      real(kind=8) :: Mold,M,Ar,f,df,tol
      integer :: iter,i,j
      common /cparam/ gamma, po,patm,rho_o
c
      pi=4d0*atan(1d0)
      gamma_m = gamma - 1.d0
      gamma_a = gamma + 1.d0
      tol=1.d-7

      do i = 1-mbc,mx+mbc
        xcell = xlower + (i-0.5d0)*dx
        if (abs(xcell-15.d0).le.10.d0) then
           radius=0.75d0-0.25d0*sin(pi*(xcell-10)/10)
        else
           radius=1.d0
        end if
        area(i)=pi*radius**2.d0
      enddo
     
      at = minval(area)

      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx   
         Ar=area(i)/at
!        Calculate Mach
         if((area(i+1)-area(i)).le.0.d0) then
            M = 0.2d0
         else
            M = 1.2d0
         endif

         do j =1,150
            Mold=M
            f=(2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0)
     &         )**(gamma_a/(2.0d0*gamma_m))/Mold - Ar
            df=(2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0))**(
     &         gamma_a/(2.0d0*gamma_m))*(Mold**2.0d0*
     &         (2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0))**
     &         (-1.0d0) -1.0d0)/Mold**2.0d0
            M=Mold-f/df
            if (dabs(M-Mold)/M.le.tol) goto 923
         enddo
923      continue        
         

!         M = 1.4d0

!         if((area(i+1)-area(i)).eq.0.d0.and.i.le.mx/3) then
!            M = 0.d0
         if((area(i+1)-area(i)).eq.0.d0.and.i.ge.mx/3) then
            M = (q(2,i-1)/q(1,i-1))/c(i-1)
         endif        

         rho=rho_o*(1.d0+gamma_m/2.d0*M**2.d0)**(-1.d0/gamma_m)
         p=po*(1.d0+gamma_m/2.d0*M**2.d0)**(-gamma/gamma_m)
         c(i)=sqrt(p/rho*gamma)
         u=M*c(i)
   
!         rho = rho_o
!         u = 0.1d0

         q(1,i) = rho   ! density
         q(2,i) = rho*u ! momentum
!         p = po - (po - patm)/((mx-0.5d0)*dx)*xcell
         q(3,i) = (p/(gamma-1.d0) + 0.5d0*u**2.d0*rho) !energy
  150    continue
c

      return
      end

