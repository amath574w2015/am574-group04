      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      real(kind=8)::gamma,po,patm,rho_o

      common /cparam/ gamma, po, patm, rho_o

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      po=4.0d0
      rho_o=1.0d0
      patm=po*0.31d0
      
                



      return
      end

