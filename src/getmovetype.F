      integer function getMoveType()
      include 'dimensions.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'printer.inc'
c######################################################
c
c     This will randomly select which type of move
c     should be made based on the probabilities loaded
c     up by the user
c
c     (C) Christopher Woods   12 March 2002
c
c######################################################

      integer ranInt
      integer ranval

      if (probTotal.le.0) then
        write(printstring,*) "Cannot run a simulation without specifying ",
     .              "the move probability ratios."
        call printLine(WARNING,printstring)
        getMoveType = 0
        return
      endif

c     generate a random number from 1 to probTotal 
      ranval = ranInt(1,probTotal)
      if (ranval.le.probProtein) then
        getMoveType = TRYPROTEIN
      else if (ranval.le.probSolute) then
        getMoveType = TRYSOLUTE
      else if (ranval.le.probSolvent) then
        getMoveType = TRYSOLVENT
      else if (ranval.le.probTitrate) then
        getMoveType = TRYTITRATE
      else if (ranval.le.probVolume) then
        getMoveType = TRYVOLUME
      else if (ranval.le.probTheta) then
        getMoveType = TRYTHETA
      else if (ranval.le.probInsertion) then
        getMoveType = TRYINSERTION
      else if (ranval.le.probDeletion) then
        getMoveType = TRYDELETION
      else if (ranval.le.probGCSolute) then
        getMoveType = TRYGCSOLUTE
      else if (ranval.le.probSample) then
        getMoveType = TRYSAMPLE
      else if (ranval.le.probSolTheta) then
        getMoveType = TRYSOLTHETA
      else
        write(printstring,*) "Could not work out the move type from ",getMoveType
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug")
        getMoveType = 0
      endif
      return
      end
