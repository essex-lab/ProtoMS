      
      subroutine readChainTemplateBBAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
c##################################################################
c
c     This routine reads in a chain template 'BBAtom' line
c
c     (C) Christopher Woods
c
c##################################################################

      character*(*) line
      integer itmpl
      character*4 at1
      integer iatm,ipar0,ipar1
      integer readInt,findChainAtom,itmp
      
c     format it 'bbatom atom-number atom-name ipar0 ipar1'
      if (nwords.ge.5) then
        ipar0 = readInt(words(4))
        ipar1 = readInt(words(5))
      else if (nwords.eq.4) then
c       format is 'bbatom atom-number atom-name ipar' ipar0=ipar1=ipar
        ipar0 = readInt(words(4))
        ipar1 = ipar0
      else if (nwords.eq.3) then
c       format is 'bbatom atom-number atom-name' ipar0=ipar1=0
        ipar0 = 0
        ipar1 = 0
      else
        call printLine(WARNING,"Incorrectly formatted chain bbatom template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     get the parameters
      iatm = readInt(words(2))
      at1 = words(3)(1:4)

c     make sure clj parameters less than 0 are interpreted as dummies
      if (ipar0.lt.0 .or. ipar0.gt.MAXCLJ) ipar0 = 0
      if (ipar1.lt.0 .or. ipar1.gt.MAXCLJ) ipar1 = 0
            
c     see if this atom is valid, or already exists
      itmp = findChainAtom(itmpl,at1)
      
      if (iatm.le.0 .or. iatm.gt.4) then
        call printLine(WARNING,"Trying to specify an invalid chain bbatom?")
        call printLine(WARNING,line)
        return
      else if (ipar0.lt.0 .or. ipar1.lt.0) then
        call printLine(WARNING,"Cannot have negative clj parameters!")
        call printLine(WARNING,line)
        return
      else if (ChainBBAtom(itmpl,iatm).ne."   ") then
        call printLine(WARNING,"Overwriting chain bbatom?")
        call printLine(WARNING,line)
      else if (itmp.ne.0) then
        call printLine(WARNING,"Cannot have a bbatom with the same name as another atom!")
        call printLine(WARNING,line)
      endif
      
      ChainBBAtom(itmpl,iatm) = at1
      ChainBBPar(itmpl,iatm,1) = ipar0
      ChainBBPar(itmpl,iatm,2) = ipar1
      
      return
      end
