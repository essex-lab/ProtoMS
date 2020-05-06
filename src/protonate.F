      logical function protonate(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c#########################################################
c
c     This function protonates iprot,ires, and returns 
c     whether it was successful...
c
c     (C) Christopher Woods      1 February 2002
c
c#########################################################

      integer iprot,ires
      character*4 inam
      logical convertRes
      
      protonate = .false.
      inam = ResNam(iprot,ires)
      
      if (inam.eq."ASP") then
c       this acid can be protonated
        protonate = convertRes(iprot,ires,"ASP","ASH")
      else if (inam.eq."GLU") then
        protonate = convertRes(iprot,ires,"GLU","GLH")
      else if (inam.eq."HIS") then
        protonate = convertRes(iprot,ires,"HIS","HIP")
      else if (inam.eq."HID") then
        protonate = convertRes(iprot,ires,"HID","HIP")
      endif

c     now that it has been protonated, move it from the base array
c     to the acid array
      call addAcid(iprot,ires)
      call takeBase(iprot,ires)
      
      return
      
      end
