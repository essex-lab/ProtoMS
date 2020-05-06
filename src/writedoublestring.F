      subroutine writeDoubleString(str,val,reslt)
      implicit none

      character*(*) str
      character*(*) reslt
      double precision val,aval
      integer lineLength

      if (val.ge.0) then
        if (val.lt.10) then
          write(reslt,"(2A,F5.3)")str(1:lineLength(str)),"+",val
        elseif (val.lt.100) then
          write(reslt,"(2A,F6.3)")str(1:lineLength(str)),"+",val
        elseif (val.lt.1000) then
          write(reslt,"(2A,F7.3)")str(1:lineLength(str)),"+",val
        endif
      else
        aval = abs(val)
        if (aval.lt.10) then
          write(reslt,"(A,F6.3)")str(1:lineLength(str)),val
        elseif (aval.lt.100) then
          write(reslt,"(A,F7.3)")str(1:lineLength(str)),val
        elseif (aval.lt.1000) then
          write(reslt,"(A,F8.3)")str(1:lineLength(str)),val
        endif
      endif

      end subroutine
