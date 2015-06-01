c
c
c     =====================================================
      subroutine bc2(meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c
c     #    aux(i,j,1)  = ax
c     #    aux(i,j,2)  = ay   where (ax,ay) is unit normal to left face
c     #    aux(i,j,4)  = bx
c     #    aux(i,j,5)  = by   where (bx,by) is unit normal to bottom face
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension mthbc(4)

c     Local Variables
      real(kind=8) rho1,u1,p1,c1,rho0,u0,p0
      real(kind=8) v0,v1,vn0,vn1,umag0,umag1
      real(kind=8) patm,uN,pN,cN,rhoN,eN1,e0
      real(kind=8) uN1,pN1,rhoN1,umagN0,umagN1
      real(kind=8) gamma_m,gamma_a,beta,theta
      integer :: i,j,ibc
      common /cparam/ gamma,gamma1, pto, patm,rho_o
      common /comic/ rhol,rhor,rhoul,rhour,el,er
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     Subsonic inflow. Define p, u, and v for the right going waves 
C     assuming u is positive. po, reservoir pressure, is allocated
c     in each ghost cell
      gamma_m = gamma - 1.d0
      gamma_a = gamma + 1.d0
      beta = gamma_a/gamma_m
      do 105 ibc=1,mbc
          do 105 j = 1-mbc,my+mbc
!            Define variables in 1st interior cells
             rho1 = q(1,1,j)
             u1 = q(2,1,j)/q(1,1,j)
             v1 = q(3,1,j)/q(1,1,j)
             theta = tan(v1/u1)
             umag1 = sqrt(u1**2.0d0+v1**2.0d0)
             p1 = gamma_m*(q(4,1,j)-0.5d0*rho1*(u1**2.0d0
     &             +v1**2.0d0))
             c1 = sqrt(gamma*p1/rho1)
           
!            U positive, three right going waves
             if(u1.gt.0.0d0) then
                p0 = pto
                umag0 = umag1-2.d0/sqrt(2.d0*gamma*gamma_m)
     &          *c1*(1.d0-p0/p1)/sqrt(1.d0+beta*p0/p1)
                u0 = umag0 !*cos(theta)
                v0 = 0.0d0 !umag0*sin(theta)
                rho0  = (1.d0+beta*p0/p1)/(p0/p1+beta)*rho1
                e0 = p0/gamma_m+0.5d0*rho0*(u0**2.d0 + v0**2.d0)

!                write(*,*) 'p1 = ',p1,'p0 = ',p0,'po = ',po
                q(1,1-ibc,j) = rho0
                q(2,1-ibc,j) = rho0*u0
                q(3,1-ibc,j) = rho0*v0
                q(4,1-ibc,j) = e0
!            U is negative, one right going wave
             else
                p0 = po
                rho0 = rho1
                u0 = u1
                v0 = v1
                e0 = p0/gamma_m+0.5d0*rho0*(u0**2.d0 + v0**2.d0)

                q(1,1-ibc,j) = rho0
                q(2,1-ibc,j) = rho0*u0
                q(3,1-ibc,j) = rho0*v0
                q(3,1-ibc,j) = e0  
              endif
  105     continue
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 j = 1-mbc, my+mbc
         do 115 ibc=1,mbc
            do 115 m=1,meqn
               q(m,1-ibc,j) = q(m,1,j)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 j = 1-mbc, my+mbc
         do 125 ibc=1,mbc
            do 125 m=1,meqn
               q(m,1-ibc,j) = q(m,mx+1-ibc,j)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 j = 1-mbc, my+mbc
         do 135 ibc=1,mbc
            do 135 m=1,meqn
               q(m,1-ibc,j) = q(m,ibc,j)
  135       continue
c     # negate the normal velocity:
      do 136 j = 1-mbc, my+mbc
         do 136 ibc=1,mbc
            q(2,1-ibc,j) = -q(2,ibc,j)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output

      gamma_m = gamma - 1.d0
      gamma_a = gamma + 1.d0
      beta = gamma_a/gamma_m

      do 205 ibc=1,mbc
          do 205 j = 1-mbc,my+mbc
!            Define variables in 1st interior cells
             rhoN = q(1,mx,j)
             uN = q(2,mx,j)/q(1,mx,j)
             vN = q(3,mx,j)/q(1,mx,j)
             theta = tan(vN/uN)
             umagN = sqrt(uN**2.0d0+vN**2.0d0)
             pN = gamma_m*(q(4,mx,j)-0.5d0*rhoN*(uN**2.0d0
     &             +vN**2.0d0))
             cN = sqrt(gamma*pN/rhoN)

!            U positive, three right going waves
             if(u1.gt.0.0d0) then
                rhoN1 = (1.d0 + beta*patm/pN)/(patm/pN+beta)*rhoN
                umagN1 = umagN+2.d0/sqrt(2.d0*gamma*gamma_m)
     &          *cN*(1.d0-patm/pN)/sqrt(1.d0+beta*patm/pN)
                uN1 = umagN1
                vN1 = 0.0d0 
                eN1 = patm/gamma_m+0.5d0*rhoN1*(uN1**2.d0 + vN1**2.d0)

                q(1,mx+ibc,j) = rhoN1
                q(2,mx+ibc,j) = rhoN1*uN1
                q(3,mx+ibc,j) = rhoN1*vN1
                q(4,mx+ibc,j) = eN1
!            U is negative, one right going wave
             else
                rhoN1 = (1.d0+beta*patm/pN)/(patm/pN+beta)*
     &              rhoN
                umagN1 = umagN+2.d0/sqrt(2.d0*gamma*
     &              gamma_m)*cN*(1.d0-patm/pN)
     &              /sqrt(1.d0+beta*patm/pN)
                uN1 = umagN1
                vN1 = 0.d0
                eN1 = patm/gamma_m+0.5d0*uN1**2.d0*rhoN1

                q(1,mx+ibc,j) = rhoN1
                q(2,mx+ibc,j) = rhoN1*uN1
                q(3,mx+ibc,j) = rhoN1*vN1
                q(3,mx+ibc,j) = eN1
              endif
  205    continue
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 j = 1-mbc, my+mbc
         do 215 ibc=1,mbc
            do 215 m=1,meqn
               q(m,mx+ibc,j) = q(m,mx,j)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 j = 1-mbc, my+mbc
         do 225 ibc=1,mbc
            do 225 m=1,meqn
               q(m,mx+ibc,j) = q(m,ibc,j)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 j = 1-mbc, my+mbc
         do 235 ibc=1,mbc
            do 235 m=1,meqn
               q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
  235       continue
c     # negate the normal velocity:
      do 236 j = 1-mbc, my+mbc
         do 236 ibc=1,mbc
            q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 jbc=1,mbc
         do 315 i = 1-mbc, mx+mbc
            do 315 m=1,meqn
               q(m,i,1-jbc) = q(m,i,1)
  315       continue
      go to 399

  320 continue
c     # periodic:  
      do 325 jbc=1,mbc
         do 325 i = 1-mbc, mx+mbc
            do 325 m=1,meqn
               q(m,i,1-jbc) = q(m,i,my+1-jbc)
  325       continue
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 jbc=1,mbc
             do 335 i = 1-mbc, mx+mbc            
               q(m,i,1-jbc) = q(m,i,jbc)
  335       continue

c     # negate the normal velocity:
c     # (for a general quadrilateral grid)
c
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
             alf = aux(4,i,1)
             beta = aux(5,i,1)
             unorm = alf*q(2,i,jbc) + beta*q(3,i,jbc)
             utang = -beta*q(2,i,jbc) + alf*q(3,i,jbc)
             unorm = -unorm
             q(2,i,1-jbc) = alf*unorm - beta*utang
             q(3,i,1-jbc) = beta*unorm + alf*utang 
  336       continue

      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 jbc=1,mbc
         do 415 i = 1-mbc, mx+mbc
            do 415 m=1,meqn
               q(m,i,my+jbc) = q(m,i,my)
  415       continue
      go to 499

  420 continue
c     # periodic:  
      do 425 jbc=1,mbc
         do 425 i = 1-mbc, mx+mbc
            do 425 m=1,meqn
               q(m,i,my+jbc) = q(m,i,jbc)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               q(m,i,my+jbc) = q(m,i,my+1-jbc)
  435       continue
c     # negate the normal velocity:
c     # (for a general quadrilateral grid)
c
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            alf = aux(4,i,my+1)
            beta = aux(5,i,my+1)
            unorm = alf*q(2,i,my+1-jbc) + beta*q(3,i,my+1-jbc)
            utang = -beta*q(2,i,my+1-jbc) + alf*q(3,i,my+1-jbc)
            unorm = -unorm
            q(2,i,my+jbc) = alf*unorm - beta*utang
            q(3,i,my+jbc) = beta*unorm + alf*utang
  436    continue
      go to 499

  499 continue

      return
      end

