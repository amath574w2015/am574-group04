      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      real(kind=8)::gamma,po,patm

      common /cparam/ gamma, p0, patm

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      p0=5.0d0
      patm=4.5d0
      
                



      return
      end

