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
c
      
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx
         q(1,i) = 0.5d0 !# density
         q(2,i) = 0.0d0 !# momentum
         q(3,i) = patm/(gamma-1.d0) + 0.5d0*q(2,i)**2.d0/q(1,i) !# energy
  150    continue
c

      return
      end

