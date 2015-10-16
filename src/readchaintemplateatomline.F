      
      subroutine readChainTemplateAtomLine(line,itmpl)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'templates.inc'
      include 'splitter.inc'
c#####################################################################
c
c     This routine reads in a chain template atom line
c
c     (C) Christopher Woods
c
c#####################################################################

      character*(*) line
      integer itmpl
      integer ipar0,ipar1,readInt,istrt
      character*4 iatm,ibnd,iang,idih
      integer nat,findChainAtom
      
c     format is "atom name ipar0 ipar1 bond angle dihedral"
      if (nwords.ge.7) then
        ipar0 = readInt(words(3))
        ipar1 = readInt(words(4))
        istrt = 5
      else if (nwords.eq.6) then
c       format is "atom name ipar bond angle dihedral" ipar0=ipar1=ipar
        ipar0 = readInt(words(3))
        ipar1 = ipar0
        istrt = 4
      else if (nwords.eq.5) then
c       format is "atom name bond angle dihedral" ipar0=ipar1=0
        ipar0 = 0
        ipar1 = 0
        istrt = 3
      else
        call printLine(WARNING,"Incorrectly formatted chain atom template line?")
        call printLine(WARNING,line)
        return
      endif
      
c     read in the data 
      iatm = words(2)(1:4)
      ibnd = words(istrt)(1:4)
      iang = words(istrt+1)(1:4)
      idih = words(istrt+2)(1:4)

c     make sure clj parameters less than 0 are interpreted as dummies
      if (ipar0.lt.0 .or. ipar0.gt.MAXCLJ) ipar0 = 0
      if (ipar1.lt.0 .or. ipar1.gt.MAXCLJ) ipar1 = 0      
            
c     now see if this atom has been used before
      nat = findChainAtom(itmpl,iatm)
      if (nat.lt.0) then
        call printLine(WARNING,"Cannot have a chain atom with the same name as a chain bbatom")
        call printLine(WARNING,line)
        return
      else if (nat.ne.0) then
        write(printstring,*) "Overwriting chain template atom ",iatm
        call printLine(WARNING,printstring)
        call printLine(WARNING,line)
      else
        nat = NChainAtoms(itmpl) + 1
        if (nat.gt.MAXCHAINATOMS) then
          write(printstring,*) WARNING,"Exceeded maximum number of chain template atoms (",
     .                           MAXCHAINATOMS,")"
          call printLine(FATAL,printstring)
          write(printstring,*)"Increase MAXCHAINATOMS and recompile"
          call closeProgram(-1,printstring)
          return
        endif
        NChainAtoms(itmpl) = nat
      endif
      
      ChainAtm(itmpl,nat,1) = iatm
      ChainAtm(itmpl,nat,2) = ibnd
      ChainAtm(itmpl,nat,3) = iang
      ChainAtm(itmpl,nat,4) = idih
      ChainAtmPar(itmpl,nat,1) = ipar0
      ChainAtmPar(itmpl,nat,2) = ipar1

      return
      end
