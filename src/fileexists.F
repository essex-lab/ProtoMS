      logical function fileExists(filename)


      character*255 filename
      integer FILE
      parameter (FILE=88)

      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)
      fileExists = .true.
      return

 9600 fileExists = .false.
      return

      end function
