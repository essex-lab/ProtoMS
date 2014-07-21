      integer function findDihedralIndex(at1,at2,at3,at4)
      include 'dimensions.inc'
      include 'parameters.inc'
      include 'printer.inc'
c#########################################################
c
c     This function returns the array index of the loaded
c     dihedral parameter for at1-at2-at3-at4, or 0 if it
c     has not been loaded
c
c     (C) Christopher Woods 
c
c#########################################################
      
      character*2 at1,at2,at3,at4
      integer i

      findDihedralIndex = 0

      if (NDihParam.gt.0) then
        call sortDihedral(at1,at2,at3,at4)
        do i=1,NDihParam
          if (at1.eq.DihName(i,1)) then
            if (at2.eq.DihName(i,2)) then
              if (at3.eq.DihName(i,3)) then
                if (at4.eq.DihName(i,4)) then
c                  write(printstring,*) at1,at2,at3,at4,i
c                  call printLine(DEBUG,printstring)
                  findDihedralIndex = i
                  return
                endif
              endif
            endif
          endif
        enddo
      endif

      return
      end
