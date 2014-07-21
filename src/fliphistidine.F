      logical function flipHistidine(iprot,ires)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
c###########################################################
c
c     This function flips a histidine from one
c     tautomer to another (HIS->HID, HID->HIS)
c
c     (C) Christopher Woods 12 July 2002
c
c###########################################################

      integer iprot,ires
      logical convertRes
      character*4 inam
      
      flipHistidine = .false.
      inam = ResNam(iprot,ires)

c     change a HIS into a HID, or a HID into a HIS      
      if (inam.eq."HIS") then
        flipHistidine = convertRes(iprot,ires,"HIS","HID")
      else if (inam.eq."HID") then
        flipHistidine = convertRes(iprot,ires,"HID","HIS")
      else
        write(printstring,*) "Can only flip HIS or HID, not ",inam
        call closeProgram(-1,printstring)
      endif
      
      return 
      end
