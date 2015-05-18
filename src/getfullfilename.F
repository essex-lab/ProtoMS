      subroutine getFullFilename(filename,fullfilename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'parallel.inc'
      include 'sep.inc'
c################################################
c
c     This subroutine appends an output foldername
c     to a filename
c
c     Samuel Genheden, February 2015
c
c#################################################
      character*(*) filename
      character*(*) fullfilename
      character*255 partfilename
      integer lineLength


      if (doMultiJAWS2) then
        write(fullfilename,"(3A)")jaws2folder(1:lineLength(jaws2folder)),PATHSEP,
     .           filename(1:lineLength(filename))
      elseif (doMultiGCMC) then
        write(fullfilename,"(3A)")gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP,
     .           filename(1:lineLength(filename))
      elseif (doLambdaRE) then
c        write(6,*)"Task = ",myid,lrefolder,filename
        write(fullfilename,"(3A)")lrefolder(1:lineLength(lrefolder)),PATHSEP,
     .           filename(1:lineLength(filename))
      elseif (doTemperatureRE) then
c        write(6,*)"Task = ",myid,lrefolder,filename
        write(fullfilename,"(3A)")trefolder(1:lineLength(trefolder)),PATHSEP,
     .           filename(1:lineLength(filename))
      else      
        if (hasoutfolder) then
          write(fullfilename,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,
     .             filename(1:lineLength(filename))
        else
          fullfilename = filename  
        endif
      endif

      end subroutine
      
