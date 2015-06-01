      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      real(kind=8)::gamma,gamma1,pto,patm,rho_o

      common /cparam/ gamma, gamma1, pto, patm, rho_o

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      gamma1 = gamma - 1.d0

      pto=4.0d0
      rho_o=1.0d0
      patm=pto*0.4d0
               



      return
      end

