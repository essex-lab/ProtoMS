      integer function openFile(filename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'parallel.inc'
c################################################
c
c     This function opens the file 'filename'
c     and returns the unit number.
c
c     If filename is already open then the 
c     unit number of the open file is returned
c
c     If the file cannot be opened then the
c     filehandle for STDOUT is returned and
c     an error is printed to the WARNING stream.
c
c     (C) Christopher Woods, October 2004
c
c#################################################

      character*(*) filename
      character*6 word
      integer lineLength,lgth,i,FILE,ierr
      character*255 fullfilename
      
      openFILE = STDOUT

      lgth = lineLength(filename)
      if (lgth.le.0) return
            
c     see if we have asked for off or err      
      if (lgth.eq.3) then
        word = filename(1:3)
        call lowercase(word)
        if (word.eq.'off') then
          openFILE = -1
          return
        else if (word.eq.'err') then
          openFILE = STDERR
          return
        endif
      endif
      
c     see if we have asked for on
      if (lgth.eq.2) then
        word = filename(1:2)
        call lowercase(word)
        if (word.eq.'on') then
          openFILE = STDOUT
          return
        endif
      endif      
      
c     see if we are asking for 'stdout' or 'stderr'
      if (lgth.eq.6) then
        word = filename(1:6)
        call lowercase(word)
        if (word.eq.'stdout') then
          openFILE = STDOUT
          return
        else if (word.eq.'stderr') then
          openFILE = STDERR
          return
        endif
      endif
      
      if (doLambdaRE) then
c        write(6,*)"Task = ",myid,lrefolder,filename
        write(fullfilename,"(3A)")lrefolder(1:lineLength(lrefolder)),"/",
     .           filename(1:lineLength(filename))
      elseif (doMultiGCMC) then
        write(fullfilename,"(3A)")gcmcfolder(1:lineLength(gcmcfolder)),"/",
     .           filename(1:lineLength(filename))
      else
        fullfilename = filename    
      endif

      lgth = lineLength(fullfilename)
c     see if this file has already been opened
      do i=1,nOpenFiles
        if (fullfilename(1:lgth).eq.openfiles(i)(1:lgth)) then
c         the file has already been opened - return the file number
          openFILE = i+40
          return
        endif
      enddo

c     this file has not been opened, so try to open it!
      FILE = nOpenFiles + 1
      if (FILE.gt.MAXOPENFILES) then
        call printLine(WARNING,"Exceeded maximum number of open files. Redirecting to STDOUT")
        call printLine(WARNING,"Increase MAXOPENFILES and recompile!")
        return
      endif
      
      nOpenFiles = FILE
      fileopened(FILE) = .false.     
      openfiles(FILE) = fullfilename
      openFILE = FILE+40
      return
      
      end
