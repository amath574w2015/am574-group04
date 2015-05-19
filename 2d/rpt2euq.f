c
c
c     =====================================================
      subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,
     &           ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision (a-h,o-z)
c
c     # Riemann solver in the transverse direction for the Euler 
c     # equations  on a curvilinear grid.
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c    
c     # Use the same idea as in rpn2 but now rotate into the direction 
c     # normal to the cell edge above or below this cell.
c
c     # Uses Roe averages 
c
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      dimension   aux1(meqn,1-mbc:maxm+mbc)
      dimension   aux2(meqn,1-mbc:maxm+mbc)
      dimension   aux3(meqn,1-mbc:maxm+mbc)
c
      parameter (maxm2 = 2002)  !# assumes at most 1000x1000 grid with mbc=2
      dimension delta(4)
      dimension u2v2(-1:maxm2),
     &       u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2),
     &       g1a2(-1:maxm2),euv(-1:maxm2) 
      dimension alf(-1:maxm2)
      dimension beta(-1:maxm2)
      dimension wave(meqn, mwaves,-1:maxm2)
      dimension    s(3,-1:maxm2)
      common /param/ gamma, gamma1
c
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpt2'
         stop
      endif
c
c
      if (ixy.eq.1) then
          inx = 4
          iny = 5
          ilenrat = 6
        else
          inx = 1
          iny = 2
          ilenrat = 3
        endif
c
c        # imp is used to flag whether wave is going to left or right,
c        # since states and grid orientation are different on each side.
c
         if (imp.eq.1) then
c            # asdq = amdq, moving to left
             ix1 = 2-mbc
	     ixm1 = mx+mbc
           else
c            # asdq = apdq, moving to right
             ix1 = 1-mbc
	     ixm1 = mx+mbc
           endif
c
c        --------------
c        # up-going:
c        --------------
c

c       # determine rotation matrix for interface above cell, using aux3
c               [ alf  beta ]
c               [-beta  alf ]
c

        do i=ix1,ixm1
c
         if (imp.eq.1) then
             i1 = i-1
           else
             i1 = i
           endif
c
           alf(i) = aux3(inx,i1)
           beta(i) = aux3(iny,i1)

!	   pres = gamma1*(ql(4,i1)  - 0.5d0*(ql(2,i1)**2 +
!     &            ql(3,i1)**2)/ql(1,i1))
           u(i) = (alf(i)*ql(2,i1) + beta(i)*ql(3,i1)) / ql(1,i1)
           v(i) = (-beta(i)*ql(2,i1) + alf(i)*ql(3,i1)) / ql(1,i1)
           pres = gamma1*(ql(4,i1)  - 0.5d0*(u(i1)**2 +
     &            v(i1)**2)*ql(1,i1))

	   enth(i) = (ql(4,i1)+pres) / ql(1,i1)
	   u2v2(i) = u(i)**2 + v(i)**2
           a2 = gamma1*(enth(i) - .5d0*u2v2(i))
           a(i) = dsqrt(a2)
	   g1a2(i) = gamma1 / a2
	   euv(i) = enth(i) - u2v2(i) 
           enddo
c
c
c     # now split asdq into waves:
c
      do 20 i = ix1,ixm1
         delta(1) = asdq(1,i) 
         delta(2) = alf(i)*asdq(2,i) + beta(i)*asdq(3,i)
         delta(3) = -beta(i)*asdq(2,i) + alf(i)*asdq(3,i)
         delta(4) = asdq(4,i) 

         a3 = g1a2(i) * (euv(i)*delta(1)
     &      + u(i)*delta(2) + v(i)*delta(3) - delta(4))
         a2 = delta(3) - v(i)*delta(1)
         a4 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a3) / (2.d0*a(i))
         a1 = delta(1) - a3 - a4

c
c        # Compute the waves.
c
         wave(1,1,i) = a1
         wave(2,1,i) = a1*(u(i)-a(i)) 
         wave(3,1,i) = a1*v(i)
         wave(4,1,i) = a1*(enth(i) - u(i)*a(i))
         s(1,i) = (u(i)-a(i))
c
         wave(1,2,i) = a3
         wave(2,2,i) = a3*u(i)
         wave(3,2,i) = a3*v(i) + a2
         wave(4,2,i) = a3*0.5d0*u2v2(i)  + a2*v(i)
         s(2,i) = u(i)
c
         wave(1,3,i) = a4
         wave(2,3,i) = a4*(u(i)+a(i))
         wave(3,3,i) = a4*v(i)
         wave(4,3,i) = a4*(enth(i)+u(i)*a(i))
         s(3,i) = (u(i)+a(i)) 

   20    continue
c
c
c    # compute flux difference bpasdq
c    --------------------------------
c
      do 40 m=1,meqn
         do 40 i=ix1,ixm1
	    bpasdq(i,m) = 0.d0
	    do 30 mw=1,mwaves
	       bpasdq(m,i) = bpasdq(m,i) + dmax1(s(mw,i),0.d0)
     &                        *wave(m,mw,i)*aux3(ilenrat,i)
   30          continue
   40       continue
c
c     # rotate momentum components:
      do 50 i=ix1,ixm1
	 bpasdq2 = alf(i)*bpasdq(2,i) - beta(i)*bpasdq(3,i)
	 bpasdq3 = beta(i)*bpasdq(2,i) + alf(i)*bpasdq(3,i)
	 bpasdq(2,i) = bpasdq2
	 bpasdq(3,i) = bpasdq3
   50    continue
c
c        --------------
c        # down-going:
c        --------------
c

c       # determine rotation matrix for interface below cell, using aux2
c               [ alf  beta ]
c               [-beta  alf ]
c
        do i=ix1,ixm1
c
         if (imp.eq.1) then
             i1 = i-1
           else
             i1 = i
           endif
c
           alf(i) = aux2(inx,i1)
           beta(i) = aux2(iny,i1)
	   pres = gamma1*(ql(4,i1)  - 0.5d0*(ql(2,i1)**2 +
     &            ql(3,i1)**2)/ql(1,i1))
           u(i) = (alf(i)*ql(2,i1) + beta(i)*ql(3,i1)) / ql(1,i1)
           v(i) = (-beta(i)*ql(2,i1) + alf(i)*ql(3,i1)) / ql(1,i1)
	   enth(i) = (ql(4,i1)+pres) / ql(1,i1)
	   u2v2(i) = u(i)**2 + v(i)**2
           a2 = gamma1*(enth(i) - .5d0*u2v2(i))
           a(i) = dsqrt(a2)
	   g1a2(i) = gamma1 / a2
	   euv(i) = enth(i) - u2v2(i) 
           enddo
c
c
c
c     # now split asdq into waves:
c
      do 80 i = ix1,ixm1
         delta(1) = asdq(1,i) 
         delta(2) = alf(i)*asdq(2,i) + beta(i)*asdq(3,i)
         delta(3) = -beta(i)*asdq(2,i) + alf(i)*asdq(3,i)
         delta(4) = asdq(4,i) 

         a3 = g1a2(i) * (euv(i)*delta(1)
     &      + u(i)*delta(2) + v(i)*delta(3) - delta(4))
         a2 = delta(3) - v(i)*delta(1)
         a4 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a3) / (2.d0*a(i))
         a1 = delta(1) - a3 - a4

c
c        # Compute the waves.

         wave(i,1,1) = a1
         wave(i,2,1) = a1*(u(i)-a(i)) 
         wave(i,3,1) = a1*v(i)
         wave(i,4,1) = a1*(enth(i) - u(i)*a(i))
         s(i,1) = (u(i)-a(i))
c
         wave(i,1,2) = a3
         wave(i,2,2) = a3*u(i)
         wave(i,3,2) = a3*v(i) + a2
         wave(i,4,2) = a3*0.5d0*u2v2(i)  + a2*v(i)
         s(i,2) = u(i)
c
         wave(i,1,3) = a4
         wave(i,2,3) = a4*(u(i)+a(i))
         wave(i,3,3) = a4*v(i)
         wave(i,4,3) = a4*(enth(i)+u(i)*a(i))
         s(i,3) = (u(i)+a(i)) 
c
   80    continue
c
c
c    # compute flux difference bmasdq
c    --------------------------------
c
      do 100 m=1,meqn
         do 100 i=ix1,ixm1
	    bmasdq(m,i) = 0.d0
	    do 90 mw=1,mwaves
	       bmasdq(m,i) = bmasdq(m,i) + dmin1(s(mw,i), 0.d0)
     &                        *wave(m,mw,i)*aux2(ilenrat,i1)
   90          continue
  100       continue
c
c     # rotate momentum components:
      do 150 i=ix1,ixm1
	 bmasdq2 = alf(i)*bmasdq(2,i) - beta(i)*bmasdq(3,i)
	 bmasdq3 = beta(i)*bmasdq(2,i) + alf(i)*bmasdq(3,i)
	 bmasdq(i,2) = bmasdq2
	 bmasdq(i,3) = bmasdq3
  150    continue
c
c
      return
      end
