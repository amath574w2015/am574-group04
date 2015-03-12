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
      real(kind=8) :: gamma, patm,xcell
      common /cparam/ gamma, po,patm

      real(kind=8) area(1-mbc:mx+mbc)

c
      area(:) = 0.d0
      call def_area(dx,mx,xlower,mbc,area)

      
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         q(1,i) = 0.8d0*area(i) !# density
         q(2,i) = q(1,i)*0.5d0 !# momentum
         q(3,i) = (patm/(gamma-1.d0) + 0.5d0*q(2,i)
     &       **2.d0/q(1,i))*area(i)
  150    continue
c

      return
      end

