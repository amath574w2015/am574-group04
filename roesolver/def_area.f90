!Area function of the nozzle

      subroutine def_area(dx,mx,xlower,mbc,area)
!      implicit none
 
      integer:: i,mx,mbc     
      real(kind=8):: pi,xlower,xcell,dx
      real(kind=8):: area(1-mbc:mx+mbc)
      real(kind=8):: radius(1-mbc:mx+mbc)
      
      pi = 4.d0*atan(1.d0)
!     Define the area
      do 150 i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        if (abs(xcell-16d0) .lt. 8d0) then
           radius(i)=1d0-0.3d0 * (1d0+cos(pi*(xcell-16d0)/8d0))
        else
           radius(i)=1.d0
        end if
        area(i)=pi*radius(i)**2.d0

150   continue
 
      end subroutine
