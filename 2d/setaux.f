c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dxc,dyc,maux,aux)
c     ============================================
c
c
c     #    aux(1,i,j)  = ax
c     #    aux(2,i,j)  = ay   where (ax,ay) is unit normal to left face
c     #    aux(3,i,j)  = ratio of length of left face to dyc
c
c     #    aux(4,i,j)  = bx
c     #    aux(5,i,j)  = by   where (bx,by) is unit normal to bottom face
c     #    aux(6,i,j)  = ratio of length of bottom face to dxc
c
c     #    aux(7,i,j)  = ratio of cell area to dxc*dyc
c     #                  (approximately Jacobian of mapping function)
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(7,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension xccorn(5),yccorn(5),xpcorn(5),ypcorn(5)
      integer i,j,k
c  

      dx2 = dxc/2.d0
      dy2 = dyc/2.d0

c
      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
c           # computational points (xc,yc) are mapped to physical
c           # coordinates (xp,yp) by mapc2p:
c
c           # lower left corner:
	    xccorn(1) = xlower + (i-1)*dxc
	    yccorn(1) = ylower + (j-1)*dyc
	    call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))

c           # upper left corner:
	    xccorn(2) = xccorn(1)
	    yccorn(2) = yccorn(1) + dyc
	    call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))
c
c           # upper right corner:
	    xccorn(3) = xccorn(1) + dxc
	    yccorn(3) = yccorn(1) + dyc
	    call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))
c
c           # lower right corner:
	    xccorn(4) = xccorn(1) + dxc
	    yccorn(4) = yccorn(1)
	    call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))
c
c           # compute normals to left and bottom side:
c
	    ax =  (ypcorn(2) - ypcorn(1))
	    ay = -(xpcorn(2) - xpcorn(1))
            anorm = dsqrt(ax*ax + ay*ay)
	    aux(1,i,j) = ax/anorm
	    aux(2,i,j) = ay/anorm
	    aux(3,i,j) = anorm/dyc
c
	    bx = -(ypcorn(4) - ypcorn(1))
	    by =  (xpcorn(4) - xpcorn(1))
            bnorm = dsqrt(bx*bx + by*by)
	    aux(4,i,j) = bx/bnorm
	    aux(5,i,j) = by/bnorm
	    aux(6,i,j) = bnorm/dxc
c
c           # compute area of physical cell from four corners:
            

	    xpcorn(5) = xpcorn(1)
	    ypcorn(5) = ypcorn(1)
	    area = 0.d0
	    do ic=1,4
	      area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) *
     &               (xpcorn(ic+1)-xpcorn(ic))
	    enddo
	    aux(7,i,j) = area / (dxc*dyc)
c
            do k = 1,7
               if(aux(k,i,j)/=aux(k,i,j)) then 
                   write(*,*) 'aux(',k,') NaN'
               endif
            enddo 
   20       continue
c
       return

       end
