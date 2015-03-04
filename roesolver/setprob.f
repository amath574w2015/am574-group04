      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cparam/ gamma

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      
                



      return
      end

