c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c     =====================================================
c
c

      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      real(kind=8) :: at,gamma,gamma_m,gamma_a
      real(kind=8) :: area(1:mx)
      real(kind=8) :: c(1-mbc:mx+mbc)
      real(kind=8) :: radius,rho,u
      real(kind=8) :: pi,patm,xc,yc,p
      real(kind=8) :: Mold,M,Ar,f,df,tol
      integer :: iter,i,j,k

      common /cparam/ gamma, po,patm,rho_o     

      pi=4d0*atan(1d0)
      gamma_m = gamma - 1.d0
      gamma_a = gamma + 1.d0
      tol=1.d-7

      do i=1,mx
        xc = xlower + (i-0.5d0)*dx
        if (abs(xc-15.d0).le.10.d0) then
           radius=0.75d0-0.25d0*sin(pi*(xc-10)/10)
        else
           radius=1.d0
        end if
        area(i)=pi*radius**2.d0
      enddo

      at = minval(area)

      do i=1,mx
        xc = xlower + (i-0.5d0)*dx
        do j=1,my
            yc = ylower + (j-0.5d0)*dy
            call mapc2p(xc,yc,xp,yp)
            Ar=area(i)/at
!           Calculate Mach
            if((area(i+1)-area(i)).le.0.d0) then
               M = 0.2d0
            else
               M = 1.2d0
            endif

            do k =1,150
               Mold=M
               f=(2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0)
     &            )**(gamma_a/(2.0d0*gamma_m))/Mold - Ar
               df=(2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0))**(
     &            gamma_a/(2.0d0*gamma_m))*(Mold**2.0d0*
     &            (2.0d0/gamma_a*(1.0d0+gamma_m/2.0d0*Mold**2.0d0))**
     &            (-1.0d0) -1.0d0)/Mold**2.0d0
               M=Mold-f/df
               if (dabs(M-Mold)/M.le.tol) goto 923
            enddo

923         continue

!         if((area(i+1)-area(i)).eq.0.d0.and.i.le.mx/3) then
!            M = 0.d0
            if((area(i+1)-area(i)).eq.0.d0.and.i.ge.mx/3) then
               M = (q(2,i-1,j)/q(1,i-1,j))/c(i-1)
            endif

            rho=rho_o*(1.d0+gamma_m/2.d0*M**2.d0)**(-1.d0/gamma_m)
            p=po*(1.d0+gamma_m/2.d0*M**2.d0)**(-gamma/gamma_m)
            c(i)=sqrt(p/rho*gamma)
            u=M*c(i)

!            rho = rho_o
!            u = 0.1d0
!            p = po - (po - patm)/((mx-0.5d0)*dx)*xcell


!           Now set the conserved variables
            q(1,i,j) = rho
            q(2,i,j) = rho*u
            q(3,i,j) = 0.0
            q(4,i,j) = (p/(gamma-1.d0) + 0.5d0*u**2.d0*rho)

        enddo
      enddo

      return
      end
