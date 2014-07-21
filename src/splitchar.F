      subroutine splitchar(line,char)
      implicit none
      include 'splitter.inc'
c####################################################################
c
c     This routine wraps up fullsplit, just as split does, but
c     this routine also allows you to specify the splitting 
c     character.
c
c     (C) Christopher Woods,   September 2004
c
c####################################################################

      character*(*) line
      character*1 char
       
      call fullsplit(line,char,1,SPLITTER_MAXSPLIT,SPLITTER_MAXLEN,
     .                            words,nwords)

      return
      end
