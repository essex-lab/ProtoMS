      subroutine initialisePrintSystem
      include 'dimensions.inc'
      include "printer.inc"
c##########################################
c
c     This routine needs to be called at 
c     the start of the program to
c     initialise the print system
c
c     (C) Christopher Woods, October 2004
c
c##########################################

      integer i

c     set the tag for each stream
      streamtag(WARNING) = 'WARNING '
      streamtag(FATAL)   = 'FATAL   '
      streamtag(HEADER)  = 'HEADER  '
      streamtag(PDB)     = 'PDB     '
      streamtag(DEBUG)   = 'DEBUG   '
      streamtag(INFO)    = 'INFO    '
      streamtag(ENERGY)  = 'ENERGY  '
      streamtag(MOVE)    = 'MOVE    '
      streamtag(SPENERGY)= 'SPENERGY'      
      streamtag(DETAIL)  = 'DETAIL  '                      
      streamtag(RESULTS) = 'RESULTS '
      streamtag(RESTART) = 'RESTART '
      streamtag(ACCPT)   = 'ACCEPT  '
      streamtag(RETI)    = 'RETI    '
      streamtag(SPLITPOT)= 'SPLITPOT'
      streamtag(PPARAMS) = 'PARAMS  '

c     now turn all of the streams on, printing to STDOUT
      do i=1,NSTREAMS     
        streamhandle(i) = -1
      enddo
      
c     set the warning and fatal streams pointing to STDERR
      streamhandle(WARNING) = STDERR
      streamhandle(FATAL) = STDERR

c     set the header, info, move, results streams pointing to STDOUT
      streamhandle(HEADER) = STDOUT
      streamhandle(INFO) = STDOUT
      streamhandle(MOVE) = STDOUT
      streamhandle(RESULTS) = STDOUT
      streamhandle(PPARAMS) = STDOUT                        

c     there are no open files
      nopenfiles = 0
      return
      end
