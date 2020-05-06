c     These are all of the variables that are used by the print system     

c     the unit numbers of STDOUT and STDERR - these are probably not portable!
      integer STDOUT,STDERR
      parameter(STDOUT=6)
      parameter(STDERR=0)
     
c     now the different types of output...
      integer HEADER,PDB,DEBUG,INFO,WARNING,ENERGY,FATAL,MOVE
      integer SPENERGY,DETAIL,RESULTS,RESTART,RETI
      integer ACCPT,SPLITPOT,PPARAMS,MCPRO

c     FATAL and WARNING are special, so are 1 and 2
      parameter(FATAL=1)
      parameter(WARNING=2)
      parameter(HEADER = 3)
      parameter(INFO = 4)
      parameter(PDB = 5)
      parameter(DEBUG = 6)
      parameter(ENERGY = 7)
      parameter(MOVE = 8)
      parameter(SPENERGY = 9)      
      parameter(DETAIL = 10)                      
      parameter(RESULTS = 11) 
      parameter(RESTART = 12)
      parameter(ACCPT = 13) 
      parameter(RETI = 14)
      parameter(SPLITPOT = 15)
      parameter(PPARAMS = 16)
      parameter(MCPRO = 17)

c     the number of different streams - must match up with what is above!
      integer NSTREAMS
      parameter(NSTREAMS=17)

c     the tag for each stream
      character*8 streamtag(NSTREAMS)

c     the maximum number of opened files
      integer MAXOPENFILES
      parameter(MAXOPENFILES=100)
      
c     the names of all of the files that are open
      character*300 openfiles(MAXOPENFILES)
c     the number of open files
      integer nopenfiles            
c     whether or not the file has been opened yet
      logical fileopened(MAXOPENFILES)
      logical wasopened(MAXOPENFILES)

c     whether or not to append to the file
      logical fileappend(MAXOPENFILES),appendto
                  
c     the file handles open for each stream - 0 means no output!
      integer streamhandle(NSTREAMS)
      
c     a string that may be used to format lines to print
      character*300 printstring

c     whether or not to do pretty printing
      logical prettyprint

      character*300 outfolder
      logical hasoutfolder

c     If the streams should be appended to
      logical appendstreams
                  
      common /PROTOPRINT/ nopenfiles,streamhandle,printstring,outfolder,streamtag,openfiles,prettyprint,
     .                    fileopened,wasopened,fileappend,appendto,hasoutfolder,appendstreams

      save /PROTOPRINT/
      
