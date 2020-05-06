      logical function getOption(option)
      include 'dimensions.inc'
      include 'keyset.inc'
      include 'splitter.inc'
c##########################################
c
c     This routine gets the value of
c     'option', splits it into words
c     and leaves it in the 'words' array.
c
c     It returns true if the option was
c     found or false if it was not found
c
c     (C) Christopher Woods, October 2004
c
c##########################################
      
      character*(*) option
      character*300 value
      
      integer lgth,lineLength,i
      
      getOption = .false.
      
c     get the length of the key string
      lgth = lineLength(option)

c     set the number of words to 0
      nwords = 0      

c     first see if the option has been set as an environmental variable
      call getenv(option(1:lgth),value)
      call split(value)
      getOption = nwords.gt.0

c     return if an environmental variable exists
      if (getOption) return      
            
c     search through all of the keys until we find it!
      do i=1,nkeys
         if (keys(i)(1:lgth) .eq. option(1:lgth)) then
c          we have found the key - split the value into words
           call split(values(i))
           getOption = (nwords.gt.0)
           return
         endif
      enddo
      
      return      
      
      end
