      subroutine readSvnTemplateInfoLine(line,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c################################################################
c
c     This reads in a solvent template line
c
c     (C) Christopher Woods
c
c################################################################

      character*(*) line
      integer itmpl
      double precision adel,rdel
      
c     format is 'info (translate adel) (rotate rdel)'
      call readInfoLine(adel,rdel)
      
      SvnTemplateADel(itmpl) = adel
      SvnTemplateRDel(itmpl) = rdel
      
      return
      end
