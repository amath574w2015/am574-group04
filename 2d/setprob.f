      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /param/ gamma,gamma1
c

      gamma = 1.4d0
      gamma1 = gamma - 1.d0

c     # read data values for this problem
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

c     # These parameters are used in qinit.f
      read(7,*) A1

      close(unit=7)

      return
      end
