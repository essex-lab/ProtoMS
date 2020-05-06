      character*16 function getAtomName(iprot,izid)
      include 'dimensions.inc'
      include 'proteins.inc'
c##################################################
c
c     Returns the atom and residue name of atom with
c     zid 'izid' in protein iprot
c
c     (C) Christopher Woods 25 November 2002
c
c##################################################

      integer iprot,izid
      integer ires,iat

      ires = Zid2ResID(iprot,izid,1)
      iat = Zid2ResID(iprot,izid,2)
      
      if (iat.gt.0) then
        write(getAtomName,10) SCNam(iprot,ires,iat),ResNam(iprot,ires),
     .                      ResNum(iprot,ires)
      else if (iat.lt.0) then
        write(getAtomName,10) BBNam(iprot,ires,-iat),ResNam(iprot,ires),
     .                      ResNum(iprot,ires)
      else
        write(getAtomName,10) "????",0,"???",0
      endif
      
10    format(a4," ",a4,"(",i5,")")

      return
      end
