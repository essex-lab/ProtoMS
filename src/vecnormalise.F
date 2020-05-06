      subroutine vecnormalise(vec)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      double precision vec(3)
      double precision length
      double precision veclength
      integer i
c##########################################

      length = veclength(vec)
      
      if (length.le.ZERO) then
        write(printstring,*) "Normalising vector of length 0.0!"
        call printLine(WARNING,printstring)
        write(printstring,*) (vec(i)," ",i=1,3)
        call printLine(WARNING,printstring)
        return
      endif        
      
      do i=1,3
        vec(i) = vec(i)/length
      enddo

      end
