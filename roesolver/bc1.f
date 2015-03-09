
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)
c     Local Variables
      real(kind=8) rho1,u1,p1,c1,rho0,u0,p0
      real(kind=8) patm,uN,pN,cN,rhoN,eN1,e0 
      real(kind=8) uN1,pN1,rhoN1
      real(kind=8) gamma_m,gamma_a,beta
      

      common /cparam/ gamma, p0, patm 

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     Subsonic inflow. Define p and u for the right going waves assuming
c     u is positive. po, reservoir pressure, is allocated in each ghost
c     cell
      gamma_m = gamma - 1.d0
      gamma_a = gamma + 1.d0
      beta = gamma_a/gamma_m

!     Define variables in 1st interior cell
      rho1 = q(1,1)
      u1 = q(2,1)/q(1,1)
      p1 = gamma_m*(q(3,1)-0.5d0*q(2,1)**2.d0/q(1,1))     
      c1 = sqrt(gamma*p1/rho1)

      do ibc=1,mbc
!         if u is positive, two right going waves
          if(u1.gt.0.0d0) then
             u0=u1-2.d0/sqrt(2.d0*gamma*gamma_m)   
     &          *c1*(1.d0-p0/p1)/sqrt(1.d0+beta*p0/p1)
             rho0 = (1.d0+beta*p0/p1)/(p0/p1+beta)*rho1 
             e0 = p0/gamma_m+0.5d0*u0**2.d0*rho0

             q(1,1-ibc) = rho0
             q(2,1-ibc) = rho0*u0
             q(3,1-ibc) = e0
!         if u is negative, one right going wave
          else 
             rho0 = rho1
             u0 = u1
             e0 = p0/gamma_m+0.5d0*u0**2.d0*rho0

             q(1,1-ibc) = rho0
             q(2,1-ibc) = rho0*u0
             q(3,1-ibc) = e0
          endif
      enddo
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
         q(2,1-ibc) = -q(2,ibc)
      end do
c     # negate the normal velocity:
      go to 199

  199 continue

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     #Outflow boundary with two cases: Subsonic or Supersonic
c      patm=0.5d0
      gamma_m = gamma-1.d0
      gamma_a = gamma + 1.d0
      beta = gamma_a/gamma_m
      

      do ibc=1,mbc
        uN = q(2,mx)/q(1,mx)
        rhoN = q(1,mx)
        cN = sqrt(gamma*pN/rhoN)
        pN = gamma_m*(q(3,mx)-0.5d0*uN
     &        **2.d0*rhoN)       

c       subsonic outflow
        if(abs(uN).lt.cN) then
!          velocity positive, one left going wave
           if(uN.gt.0.d0) then
             rhoN1 = (1.d0+beta*patm/pN)/
     &          (patm/pN+beta)*rhoN
             uN1 = uN
             eN1 = pN/gamma_m+0.5d0*uN1**2.d0*rhoN1

             q(1,mx+ibc)=rhoN1
             q(2,mx+ibc)=rhoN1*uN1
             q(3,mx+ibc)=eN1

!          velocity negative, two left going waves
           else
             rhoN1 = (1.d0+beta*patm/pN)/(patm/pN+beta)*
     &          rhoN
             uN1 = uN+2.d0/sqrt(2.d0*gamma*
     &          gamma_m)*cN*(1.d0-patm/pN)
     &          /sqrt(1.d0+beta*patm/pN)
             eN1 = patm/gamma_m+0.5d0*uN1**2.d0*rhoN1

             q(1,mx+ibc) = rhoN1
             q(2,mx+ibc) = rhoN1*uN1
             q(3,mx+ibc) = eN1
           endif
c       supersonic outflow, just extrapolate
        else
           do m=1,meqn
              q(m,mx+ibc)=q(m,mx)
           enddo
        endif
      enddo
      go to 299

  210 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end

