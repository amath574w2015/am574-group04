      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cparam/ gamma, p0, patm

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      p0=6.d0
      patm=4.d0
      
                



      return
      end

