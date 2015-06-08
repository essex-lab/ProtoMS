      subroutine parseArgs
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'printer.inc'
      include 'version.inc'
      include 'parallel.inc'
#ifdef WITH_MPI_
      include 'mpif.h'
#endif

      integer ierr      
      character*300 arg

      call getarg(1,arg)
      call lowercase(arg)
      call split(arg)

      if (nwords.gt.0) then
        if (words(1).eq.'-h' .or. words(1).eq.'--help') then
          if (myid.eq.0) then
            print *,progname," ",progversion
            print *," "
            print *,"USAGE:"
            print *,"  ",exename," -h / --help : View this help"
            print *,"  ",exename," -v / --version : View version information"
            print *,"  ",exename," -l / --license : View license information"
            print *,"  ",exename," -c / --coders : List everyone who has helped program ",exename
            print *,"  ",exename," filename : Run using a command file"
          endif
#ifdef WITH_MPI_
          call MPI_FINALIZE(ierr)
#endif
          call exit(0)
        else if (words(1).eq.'-v' .or. words(1).eq.'--version') then
          if (myid.eq.0) print *,progname," ",progversion
#ifdef WITH_MPI_
          call MPI_FINALIZE(ierr)
#endif
          call exit(0)
        else if (words(1).eq.'-l' .or. words(1).eq.'--license') then
          if (myid.eq.0) print *,"GPL license..."
#ifdef WITH_MPI_
          call MPI_FINALIZE(ierr)
#endif
          call exit(0)
        else if (words(1).eq.'-c' .or. words(1).eq.'--coders') then
          if (myid.eq.0) then
            print *,"Huge thanks to everyone who has helped produce this code!"
            print *,"Thanks especially to the coding superstars;"
            print *,"Julien Michel, Michael Bodnarchuck, Samuel Genheden, Ana Isabel Cabedo Martinez, Gregory Ross, 
     . Chris Cave-Ayland, Richard Bradshaw, James Graham"
          endif
#ifdef WITH_MPI_
          call MPI_FINALIZE(ierr)
#endif
          call exit(0)
        else
c         Parse the command file          
          call parseFile(arg)
        endif
      endif
      
      return
      end
