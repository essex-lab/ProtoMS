      character*10 function getSolAtomName(isol,izid)
      include 'dimensions.inc'
      include 'solutes.inc'
c####################################################
c
c     Returns the atom and residue name of atom with
c     zid 'izid' in solute isol
c
c     (C) Christopher Woods, October 2004
c
c####################################################

      integer isol,izid
      
      if (izid.le.0) then
        if (izid.eq.-1) then
          write(getSolAtomName,10) " dm1"," dum"
        else if (izid.eq.-2) then
          write(getSolAtomName,10) " dm2"," dum"
        else if (izid.eq.-3) then
          write(getSolAtomName,10) " dm3"," dum"
        else
          write(getSolAtomName,10) "????","????"
        endif
      else
        write(getSolAtomName,10) SolAtomNam(isol,izid),SolResNam(isol,SolResidue(isol,izid))
      endif
      
10    format(a4,"(",a4,")")
      
      return
      end
