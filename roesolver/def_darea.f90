!     Derivative of area function
      subroutine def_darea(dx,mx,xlower,mbc,darea)
      implicit none

      integer:: i,mx,mbc
      real(kind=8):: pi,xlower,xcell,dx
      real(kind=8):: darea(1-mbc:mx+mbc)
      
      pi = 4.d0*atan(1.d0)
!     Define the area
      do 151 i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        if (abs(xcell-16d0) .lt. 8d0) then
           darea(i)=0.74022d0 * (1d0-0.3d0 * (1d0+cos(1d0/8d0 &
               *pi*(xcell-16d0))))*sin(1d0/8d0*pi*(xcell-16d0))
        else
           darea(i)=0.d0
        end if

151   continue

      end subroutine
