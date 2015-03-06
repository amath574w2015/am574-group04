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
c
      do 150 i=1,mx
             q(1,i) = 1.d0 !# density
             q(2,i) = 0.5d0 !# momentum
             q(3,i) = 1.d0 !# energy
  150    continue
c
      return
      end

