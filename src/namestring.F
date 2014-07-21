      character*13 function nameString(iprot,izat)
      include 'dimensions.inc'
      include 'proteins.inc'
c########################################################
c
c     Return a pretty formatted name string identifying
c     the atom izat in protein iprot
c     (in the form 'ATM-GLY(  12)')
c
c     (C) Christopher Woods    12 March 2002
c
c#######################################################

      integer iprot,izat
      integer ires,iatm,iresnum
      character*4 anam,rnam
      
      ires = Zid2ResID(iprot,izat,1)
      iatm = Zid2ResID(iprot,izat,2)
      
      if (iatm.lt.0) then
        anam = BBNam(iprot,ires,-iatm)
      else if (iatm.eq.0) then
        anam = "???"
      else
        anam = SCNam(iprot,ires,iatm)
      endif
      
      iresnum = ResNum(iprot,ires)
      rnam = ResNam(iprot,ires)
      
      write(nameString,10) anam,rnam,iresnum
10    format(a3,"-",a3,"(",i4,")")
      return
      end
