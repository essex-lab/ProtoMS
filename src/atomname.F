      character*4 function atomName(iprot,izid)
      include 'dimensions.inc'
      include 'proteins.inc'
c#########################################################
c
c     This returns the name of atom izid in protein
c     iprot
c
c     (C) Christopher Woods, November 2004
c
c#########################################################

      integer iprot,izid
      integer iat,ires
      
      iat = Zid2ResID(iprot,izid,2)
      ires = Zid2ResID(iprot,izid,1)
      
      if (iat.gt.0) then
        atomName = SCNam(iprot,ires,iat)
      else if (iat.lt.0) then
        atomName = BBNam(iprot,ires,-iat)
      else
        atomName = "????"
      endif
      
      return
      end
