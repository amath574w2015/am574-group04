      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      real(kind=8)::gamma,po,patm

      common /cparam/ gamma, po, patm

c     # Set the common parameters used for the Euler Eqns
c     # Calorically perfect gas is assumed
c
c
      gamma = 1.4d0
      po=6.0d0
      patm=2.0d0
      
                



      return
      end

