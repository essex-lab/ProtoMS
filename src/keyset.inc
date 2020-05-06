
C     Variables used to hold a dictionary of input options

c     maximum number of keys in the dictionary
      integer MAXKEYS
      parameter(MAXKEYS=10000)

c     how many keys are in the dictionary
      integer nkeys
      
c     array of all keys
      character*64 keys(MAXKEYS)
c     array of all values
      character*256 values(MAXKEYS)

      common /KEYVARS/ keys,values,nkeys 
      save /KEYVARS/
            
