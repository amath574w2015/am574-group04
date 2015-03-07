
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
      real(kind=8) p,ubc,po,patm
      real(kind=8) gamma_m,gamma_a,beta
      

      common /cparam/ gamma, po, patm 

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
      do ibc=1,mbc
          p=gamma_m*(q(3,2-ibc)-0.5d0*q(2,2-ibc)**2.d0/q(1,2-ibc))
!         if u is positive, two right going waves
          if(q(2,2-ibc)/q(1,2-ibc).lt.0.0d0) then
             ubc=q(2,2-ibc)/q(1,2-ibc)-2.d0/sqrt(2.d0*gamma*gamma_m)   
     &          *sqrt(gamma*p/q(1,2-ibc))*(1.d0-po/p) 
     &          /sqrt(1.d0+beta*po/p)
             q(1,1-ibc) = (1.d0+beta*po/p)/(po/p+beta)*q(1,2-ibc)
             q(2,1-ibc) = q(1,2-ibc)*ubc
             q(3,1-ibc) = po/gamma_m+0.5d0*q(2,1-ibc)**2.d0/q(1,1-ibc)
!         if u is negative, one right going wave
          else 
             q(1,1-ibc) = q(1,2-ibc)
             q(2,1-ibc) = q(2,2-ibc)
             q(3,1-ibc) = po/gamma_m+0.5d0*q(2,1-ibc)**2.d0/q(1,1-ibc)
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
      do ibc=1,mbc
        p=gamma_m*(q(3,mx+ibc-1)-0.5d0*q(2,mx+ibc-1)
     &        **2.d0/q(1,mx+ibc-1))
c       subsonic outflow
        if(abs(q(2,mx+ibc-1)/q(1,mx+ibc-1)).lt.sqrt(gamma*p
     &        /q(1,mx+ibc-1))) then
!          velocity positive, one left going wave
           if(q(2,mx+ibc-1)/q(1,mx+ibc-1).lt.0.d0) then
             q(1,mx+ibc)=q(1,mx+ibc-1)
             q(2,mx+ibc)=q(2,mx+ibc-1)
             q(3,mx+ibc)=patm/gamma_m+0.5d0*q(2,mx+ibc-1)**
     &         2.d0/q(1,mx+ibc-1)
!          velocity negative, two left going waves
           else
             ubc=q(2,mx+ibc-1)/q(1,mx+ibc-1)+2.d0/sqrt(2.d0*gamma*
     &          gamma_m)*sqrt(gamma*p/q(1,mx+ibc-1))*(1.d0-patm/p)
     &          /sqrt(1.d0+beta*patm/p)
             q(1,mx+ibc) = (1.d0+beta*patm/p)/(patm/p+beta)*
     &          q(1,mx+ibc-1)
             q(2,mx+ibc) = q(1,mx+ibc-1)*ubc
             q(3,mx+ibc) = patm/gamma_m+0.5d0*q(2,mx+ibc-1)**
     &          2.d0/q(1,mx+ibc-1)
           endif
c       supersonic outflow
        else
           do m=1,meqn
              q(m,mx+ibc)=q(m,mx+ibc-1)
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

