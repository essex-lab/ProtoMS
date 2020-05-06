      subroutine readSolTemplateInfoLine(line,itmpl)
      include 'dimensions.inc'
      include 'templates.inc'
c#############################################################
c
c     Read in the solute info (translate and rotation values)
c
c     (C) Christopher Woods
c 
c#############################################################

      character*(*) line
      integer itmpl
      double precision adel,rdel
      
c     format is 'info (translate adel) (rotate rdel)'      
      call readInfoLine(adel,rdel)
      
      MolTemplateAdel(itmpl) = adel
      MolTemplateRdel(itmpl) = rdel
      
      return 
      end
