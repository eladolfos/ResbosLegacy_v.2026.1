      program add_2_hbooks
      implicit none
      integer nwpawc
      real*4 h
      character*80 infilename(2)
      character*80 outfilename
      parameter (nwpawc=50000000)
      common/pawc/h(nwpawc)

      character*80 string
      integer right


c      print *, 'hbook+: a small program to add two HBOOK files'
c      print *, 'Author: P. Nadolsky, March 2001'
c      print *, 'Usage: hbook+ in_hbook_1 in_hbook_2 out_hbook'

      if (iargc().ne.3) then
        print *,'Exiting hbook+: Wrong number of parameters =',iargc(),
     >       ' (should be three)'
        stop
      endif 

      call getarg (1,  infilename(1))
      call getarg (2,  infilename(2))
      call getarg (3,  outfilename)

      call hlimit(nwpawc)      
      call hmerge(2,infilename,outfilename)

      End

        
