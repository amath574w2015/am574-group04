c
c     =====================================================
      subroutine mapc2p(xc,yc,xp,yp)
c     =====================================================
c
c     # on input,  (xc,yc) is a computational grid point
c     # on output, (xp,yp) is corresponding point in physical space
c
      implicit double precision (a-h,o-z)
c
      pi = acos(-1.d0)
      rad1 = 0.75d0-0.25d0*sin(pi*(xc-10)/10)
C      rad1 = 1.d0 - 0.3d0*(1.d0 + cos(pi*(xc-16.d0)/8.d0))
      if (abs(xc-15.d0) <= 10.d0) then
           radius = rad1
        else
           radius = 1.00d0
        endif
      xp = xc
      yp = yc * radius
      return
      end
