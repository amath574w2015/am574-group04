!
!
!     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr, &
                       wave,s,amdq,apdq)
!     =====================================================
!
!     # Roe-solver for the Euler equations on a curvilinear grid
!     # mwaves = 3
!
!     # solve Riemann problems along one slice of data.
!
!     # On input, ql contains the state vector at the left edge of each cell
!     #           qr contains the state vector at the right edge of each cell
!
!     # This data is along a slice in the x-direction if ixy=1 
!     #                            or the y-direction if ixy=2.
!     # On output, wave contains the waves, s the speeds, 
!     # and amdq, apdq the decomposition of the flux difference
!     #   f(qr(i-1)) - f(ql(i))  
!     # into leftgoing and rightgoing parts respectively.
!     # With the Roe solver we have   
!     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
!     # where A is the Roe matrix.  An entropy fix can also be incorporated
!     # into the flux differences.
!
!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     # From the basic clawpack routines, this routine is called with ql = qr
!
!
      implicit double precision (a-h,o-z)
!
    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension s(mwaves, 1-mbc:maxm+mbc)
    dimension ql(meqn, 1-mbc:maxm+mbc)
    dimension qr(meqn, 1-mbc:maxm+mbc)
    dimension apdq(meqn, 1-mbc:maxm+mbc)
    dimension amdq(meqn, 1-mbc:maxm+mbc)
    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)
!
!     local arrays -- common block comroe is passed to rpt2eu
!     ------------
      parameter (maxm2 = 2002)  !# assumes at most 1000x1000 grid with mbc=2
      dimension delta(4)
      logical efix
      dimension u2v2(-1:maxm2), &
            u(-1:maxm2),v(-1:maxm2),enth(-1:maxm2),a(-1:maxm2), &
            g1a2(-1:maxm2),euv(-1:maxm2) 
      dimension q2l(-1:maxm2), q2r(-1:maxm2)
      dimension q3l(-1:maxm2), q3r(-1:maxm2)
      dimension alf(-1:maxm2)
      dimension beta(-1:maxm2)
      

      common /cparam/ gamma, gamma1
!
      data efix /.false./    !# use entropy fix for transonic rarefactions
!     #### not yet working on curvilinear grid!
!
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpn2'
         stop
         endif
!
!     # rotate the velocities q(2) and q(3) so that it is aligned with grid
!     # normal.  The normal vector for the face at the i'th Riemann problem
!     # is stored in the aux array
!     # in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
!     # length of the cell side to the length of the computational cell
!     # is stored in aux(3) or aux(6) respectively.


      if (ixy.eq.1) then
          inx = 1
          iny = 2
          ilenrat = 3
        else
          inx = 4
          iny = 5
          ilenrat = 6
        endif

!       # determine rotation matrix
!               [ alf  beta ]
!               [-beta  alf ]

!       # note that this reduces to identity on standard cartesian grid

        do i=2-mbc,mx+mbc
           alf(i) = auxl(inx,i)
           beta(i) = auxl(iny, i)
!           write(*,*) 'alf = ',alf(i),'beta = ',beta(i)

           if(alf(i)/=alf(i)) write(*,*) 'alf(',i,') is NaN'
           if(beta(i)/=beta(i)) write(*,*) 'beta(',i,') is NaN'

           if(ql(2,i)/=ql(2,i)) write(*,*) 'ql(2,',i,') is NaN= ',ql(2,i)
           if(qr(2,i)/=qr(2,i)) write(*,*) 'qr(2,',i,') is NaN= ',qr(2,i)
           if(ql(3,i)/=ql(3,i)) write(*,*) 'ql(3,',i,') is NaN= ',ql(3,i)
           if(qr(3,i)/=qr(3,i)) write(*,*) 'qr(3,',i,') is NaN= ',qr(3,i)

           q2l(i) = alf(i)*ql(2, i) + beta(i)*ql(3, i)
           q2r(i-1) = alf(i)*qr(2,i-1) + beta(i)*qr(3,i-1)
           q3l(i) = -beta(i)*ql(2,i) + alf(i)*ql(3,i)
           q3r(i-1) = -beta(i)*qr(2,i-1) + alf(i)*qr(3,i-1)

           if(q2l(i)/=q2l(i)) write(*,*) 'q2l(',i,') is Nan'
           if(q2r(i-1)/=q2r(i-1)) write(*,*) 'q2r(',i-1,') is Nan'
           if(q3l(i)/=q3l(i)) write(*,*) 'q3l(',i,') is Nan'
           if(q3r(i-1)/=q3r(i-1)) write(*,*) 'q3r(',i-1,') is Nan'    
        enddo


      do 10 i = 2-mbc, mx+mbc
         rhsqrtl = dsqrt(qr(1,i-1))
         rhsqrtr = dsqrt(ql(1,i-1))
         pl = gamma1*(qr(4,i-1) - 0.5d0*(q3r(i-1)**2.d0 + &
             q3r(i-1)**2.d0)/qr(1,i-1))
         pr = gamma1*(ql(4,i) - 0.5d0*(q3l(i)**2.d0 +  &
             q3l(i)**2.d0)/ql(1,i))
         if(pr/=pr) write(*,*) 'pr NaN'
         if(pl/=pl) write(*,*) 'pl NaN'
         rhsq2 = rhsqrtl + rhsqrtr
         if(rhsq2/=rhsq2) write(*,*) 'rho neg'
         u(i) = (q2l(i)/rhsqrtr + q2r(i-1)/rhsqrtl) / rhsq2
         v(i) = (q3l(i)/rhsqrtr + q3r(i-1)/rhsqrtl) / rhsq2
         enth(i) = (((qr(4,i-1)+pl)/rhsqrtl       &
                  + (ql(4,i)+pr)/rhsqrtr)) / rhsq2
         u2v2(i) = u(i)**2.d0 + v(i)**2.d0
         a2 = gamma1*(enth(i) - .5d0*u2v2(i))
         if(a2 <= 0) then
             write(*,*) 'a2 = ',a2,'gamma1 = ',gamma1
             write(*,*) 'enth(i) = ',enth(i), 'u2v2 = ',u2v2(i)
             write(*,*) 'u = ',u(i),'v = ',v(i)
         endif
         a(i) = dsqrt(a2)
         g1a2(i) = gamma1 / a2
         euv(i) = enth(i) - u2v2(i) 

         if(g1a2(i)/=g1a2(i)) write(*,*) 'g1a2 nan'
         if(euv(i)/=euv(i)) write(*,*) 'euv NaN'
   10    continue
!
!
!     # now split the jump in q at each interface into waves
!
!     # find a1 thru a4, the coefficients of the 4 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(1,i) - qr(1,i-1)
         delta(2) = q2l(i) - q2r(i-1)
         delta(3) = q3l(i) - q3r(i-1)
         delta(4) = ql(4,i) - qr(4,i-1)
         a3 = g1a2(i) * (euv(i)*delta(1)   &
           + u(i)*delta(2) + v(i)*delta(3) - delta(4))
         a2 = delta(3) - v(i)*delta(1)
         a4 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a3) / (2.d0*a(i))
         a1 = delta(1) - a3 - a4
!
!        # Compute the waves.
!        # Note that the 2-wave and 3-wave travel at the same speed and 
!        # are lumped together in wave(2,.,.).  The 4-wave is then stored in
!        # wave(.,3,.).
!

         wave(1,1,i) = a1
         wave(2,1,i) = a1*(u(i)-a(i)) 
         wave(3,1,i) = a1*v(i)
         wave(4,1,i) = a1*(enth(i) - u(i)*a(i))
         s(1,i) = (u(i)-a(i))
!
         wave(1,2,i) = a3
         wave(2,2,i) = a3*u(i)
         wave(3,2,i) = a3*v(i) + a2
         wave(4,2,i) = a3*0.5d0*u2v2(i)  + a2*v(i)
         s(2,i) = u(i)
!
         wave(1,3,i) = a4
         wave(2,3,i) = a4*(u(i)+a(i))
         wave(3,3,i) = a4*v(i)
         wave(4,3,i) = a4*(enth(i)+u(i)*a(i))
         s(3,i) = (u(i)+a(i)) 

!	 Check for NaNs in waves
         do k = 1,4
            do j = 1,3
               if (wave(k,j,i)/=wave(k,j,i)) then
                  write(*,*) 'wave(',k,',',j,',',i,') is NaN'
               endif
            enddo
         enddo

   20    continue


!
!
!    # compute flux differences amdq and apdq.
!    ---------------------------------------
!
      if (efix) go to 110
!
!     # no entropy fix
!     ----------------
!
      do 80 i=2-mbc, mx+mbc
         do 80 mw=1,mwaves
!
!           # scale wave speeds by ratio of cell side length to dxc:
            s(mw,i) = s(mw,i) * auxl(ilenrat,i)
!
!           # rotate momentum components of waves back to x-y:
            wave2 = alf(i)*wave(2,mw,i) - beta(i)*wave(3,mw,i)
            wave3 = beta(i)*wave(2,mw,i) + alf(i)*wave(3,mw,i)
            wave(2,mw,i) = wave2
            wave(3,mw,i) = wave3
   80       continue
! 
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves
!
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
                 if (s(mw,i) .lt. 0.d0) then
                   amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
   90          continue
  100       continue
      go to 900   
!
!-----------------------------------------------------
!
!  Need to fix for curvilinear!
!
  110 continue
!
!     # With entropy fix
!     ------------------
!
!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.
!
      do 200 i = 2-mbc, mx+mbc
!
!        # check 1-wave:
!        ---------------
!
         rhoim1 = qr(1,i-1)
         pim1 = gamma1*(qr(4,i-1) - 0.5d0*u2v2(i)*rhoim1)
         cim1 = dsqrt(gamma*pim1/rhoim1)
!      # speed of left-most signal: s0 = u-c in left state (cell i-1)
         s0 = q2r(i-1)/rhoim1 - cim1

!        # check for fully supersonic case:
         if (s0.ge.0.d0 .and. s(1,i).gt.0.d0)  then
!            # everything is right-going
             do 60 m=1,meqn
                amdq(m,i) = 0.d0
   60           continue
             go to 200 
             endif
!
         rho1 = qr(1,i-1) + wave(1,1,i)
         rhou1 = q2r(i-1) + wave(2,1,i)
         rhov1 = q3r(i-1) + wave(3,1,i)
         en1 = qr(4,i-1) + wave(4,1,i)
         p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2)/rho1)
         c1 = dsqrt(gamma*p1/rho1)
         s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
!            # transonic rarefaction in the 1-wave
             sfract = s0 * (s1-s(1,i)) / (s1-s0)
         else if (s(1,i) .lt. 0.d0) then
!	     # 1-wave is leftgoing
             sfract = s(1,i)
         else
!	     # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
         endif
         do 120 m=1,4
            amdq(i,m) = sfract*wave(m,1,i)
  120       continue
!
!        # check 2-wave:
!        ---------------
!
         if (s(2,i) .ge. 0.d0) go to 200  !# 2- and 3- waves are rightgoing
         do 140 m=1,4
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
  140       continue
!
!        # check 3-wave:
!        ---------------
!
         rhoi = ql(1,i)
         pi = gamma1*(ql(4,i) - 0.5d0*u2v2(i)*rhoi)
         ci = dsqrt(gamma*pi/rhoi)
         s3 = q2l(i)/rhoi + ci     !# u+c in right state  (cell i)
!
         rho2 = ql(1,i) - wave(1,3,i)
         rhou2 = q2l(i) - wave(2,3,i)
         rhov2 = q3l(i) - wave(3,3,i)
         en2 = ql(4,i) - wave(4,3,i)
         p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2)/rho2)
         c2 = dsqrt(gamma*p2/rho2)
         s2 = rhou2/rho2 + c2   !# u+c to left of 3-wave
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0) then
!            # transonic rarefaction in the 3-wave
             sfract = s2 * (s3-s(3,i)) / (s3-s2)
         else if (s(3,i) .lt. 0.d0) then
!            # 3-wave is leftgoing
             sfract = s(3,i)
         else 
!            # 3-wave is rightgoing
             go to 200
         endif
!
         do 160 m=1,4
            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
  160       continue
  200    continue
!
      do 190 i=2-mbc, mx+mbc
         do 180 mw=1,mwaves
!
!           # scale wave speeds by ratio of cell side length to dxc:
            s(mw,i) = s(mw,i) * auxl(ilenrat,i)
!
!           # rotate momentum components of waves back to x-y:
            wave2 = alf(i)*wave(2,mw,i) - beta(i)*wave(3,mw,i)
            wave3 = beta(i)*wave(2,mw,i) + alf(i)*wave(3,mw,i)
            wave(2,mw,i) = wave2
            wave(3,mw,i) = wave3
  180       continue
!
!        # flux difference must also be rotated and scaled:
         amdq2 = (alf(i)*amdq(2,i) - beta(i)*amdq(3,i)) 
         amdq3 = (beta(i)*amdq(2,i) + alf(i)*amdq(3,i))
!
         amdq(1,i) = amdq(1,i) * auxl(ilenrat,i)
         amdq(2,i) = amdq2     * auxl(ilenrat,i)
         amdq(3,i) = amdq3     * auxl(ilenrat,i)
         amdq(4,i) = amdq(4,i) * auxl(ilenrat,i)
  190    continue
! 
!
!     # compute the rightgoing flux differences:
!     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
!
      do 220 m=1,4
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
               df = df + s(mw,i)*wave(m,mw,i)
  210          continue
            apdq(m,i) = df - amdq(m,i)
  220       continue
!
  900 continue
      return
      end
