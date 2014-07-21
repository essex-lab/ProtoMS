      subroutine parseRetiArgs
      include 'dimensions.inc'
      include 'reti.inc'
      include 'splitter.inc'
      include 'printer.inc'
c#######################################################
c
c     Parse commande line arguments
c
c     Samuel Genheden, November 2012
c     From parseargs.F
c
c#######################################################
      
      character*300 arg
      
      call getarg(1,arg)
      call lowercase(arg)
      call split(arg)

      if (nwords.gt.0) then
        if (words(1).eq.'-h' .or. words(1).eq.'--help') then
          print *,"reti 1.0"
          print *," "
          print *,"USAGE:"
          print *,"   reti -h / --help : View this help"
          print *,"   reti filename : Run using a command file"
          call exit(0)
        else
c         Parse the command file          
          call parseFile(arg)
        endif
      endif
      
      return
      end
