      subroutine readGCSTemplateInfoLine(line,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c################################################################
c
c     This reads in a solvent template line
c
c     MSB 31/03/09
c
c################################################################

      character*(*) line
      integer itmpl
      double precision adel,rdel
      
c     format is 'info (translate adel) (rotate rdel)'
      call readInfoLine(adel,rdel)
      
      GCSTemplateADel(itmpl) = adel
      GCSTemplateRDel(itmpl) = rdel
      
      return
      end
