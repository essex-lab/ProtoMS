      logical function lambdaMove(updatenrg,movnum,movtyp)
      include 'dimensions.inc'
      include 'printer.inc'
c###############################################################
c
c     This function performs everything necessary for a 
c     successful lambda move.
c
c     (C) Christopher Woods   14 February 2002
c
c###############################################################

      integer updatenrg,movnum,movtyp

      call printLine(WARNING,"Lambda move not implemented yet!")
      call printLine(WARNING,"Maybe you would like to implement it?")

      write(printstring,*) "The user wanted to run movnum = ",movnum,
     .                     "updatenrg = ",updatenrg," and movtyp = ",movtyp
      call printLine(WARNING,printstring)

      lambdaMove = .false.
      
      return
      end
