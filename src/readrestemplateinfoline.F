      subroutine readResTemplateInfoLine(line,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c######################################################
c
c     This routine reads in additional info about
c     the residue (amount to translate and rotate by)
c
c     (C) Christopher Woods
c
c######################################################

      character*(*) line
      integer itmpl
      double precision adel,rdel

c     format is "info (translate 0.1) (rotate 3.0)  - these are optional      
      call readInfoLine(adel,rdel)
      
      TemplateResADel(itmpl) = adel
      TemplateResRDel(itmpl) = rdel
      
      return
      end
