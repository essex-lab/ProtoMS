      logical function deprotonate(iprot,ires)
      include 'dimensions.inc'
      include 'proteins.inc'
c############################################################
c
c     This function will attempt to deprotonate ires of iprot,
c     and will return true if successful!
c
c     (C) Christopher Woods
c
c############################################################

      integer iprot,ires
      integer iran,ranInt
      character*4 inam
      logical convertRes
      
      deprotonate = .false.
      
      inam = ResNam(iprot,ires)
      
c     you can add other transformations to this list below!
      
      if (inam.eq."ASH") then
c       can deprotonate this protonated acid
        deprotonate = convertRes(iprot,ires,"ASH","ASP")
      else if (inam.eq."GLH") then
c       can deprotonate this protonated acid
        deprotonate = convertRes(iprot,ires,"GLH","GLU")
      else if (inam.eq."HIP") then
c       can deprotonate this base
c       random choice - use an inlined function to generate a random
c       number between 1 and 2
        iran = ranInt(1,2)
        
        if (iran.eq.1) then
          deprotonate = convertRes(iprot,ires,"HIP","HIS")
        else
          deprotonate = convertRes(iprot,ires,"HIP","HID")
        endif
      endif      
      
c     now its been deprotonated, move the residue from the acid array
c     to the base array...
      call takeAcid(iprot,ires)
      call addBase(iprot,ires)
      
      return
      
      end
