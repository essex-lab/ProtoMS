      subroutine findDoubleInList(element,list,listlength,ind)
c      implicit none

      include 'dimensions.inc'
      include 'parallel.inc'

c#################################################
c
c     This function stores 0 in 'ind' if 'element'
c     is not in 'list'.
c     If 'element' is in 'list', it stores
c     the index of the element in the list.
c
c     Ana I. Cabedo Martinez, 2015
c
c#################################################

      integer i,listlength,ind
      double precision element, list(listlength)

      ind = 0
      do i=1,listlength
        if (list(i).eq.element) ind = i
      enddo
      
      return

      end
