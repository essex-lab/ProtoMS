
      subroutine setEquilibriumGeometry(ityp,imol,ires,igeom)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'printer.inc'
      include 'parameters.inc'
c##################################################################
c
c     This routine sets the internal degree of freedom
c     (ityp, BONDS or ANGLES) of molecule (ityp, PROTEIN,SOLUTE)
c     number imol (residue ires if a protein) to the 
c     equilibrium values taken from the harmonic energy terms
c     (r0 or theta0)
c
c     (C) Christopher Woods, January 2005
c
c##################################################################

      integer ityp,imol,ires,igeom
      integer i,nats,ibnd,iang,iparam
      integer specifiedBond,specifiedAngle
      
      if (ityp.eq.PROTEIN) then
        call printLine(WARNING,"Not implemented yet!")
        return
      
      else if (ityp.eq.SOLUTE) then
        if (imol.le.0 .or. imol.gt.NSolutes) then
          call printLine(WARNING,"Requesting an invalid solute!")
          return
        endif
        
        if (igeom.eq.BONDS) then
c         run through each atom in the solute and reset its zmatrix bond to
c         its equilibrium value
          nats = NSolAtoms(imol)
          
          do i=1,nats
c           see if there is a bond specified between this atom and its z-matrix bond atom
            ibnd = specifiedBond(SOLUTE,imol,i,ZSolMat(imol,i,1))
            if (ibnd.gt.0) then
c             there is! Set the bond to the equilibrium value
              iparam = SolBond(imol,ibnd,3)
              if (iparam.gt.0) then
                ZSolVal(imol,i,1) = BndParam(iparam,2)
              else if (iparam.lt.0) then
                call prettyLine(WARNING,"Variable BONDS")
              endif
            endif
          enddo
        else if (igeom.eq.ANGLES) then
c         run through each atom in the solutes and reset its zmatrix angle to its
c         equilibrium value
          nats = NSolAtoms(imol)
          
          do i=1,nats
c           see if there is an angle specified between this atom and its z-matrix atoms
            iang = specifiedAngle(SOLUTE,imol,i,ZSolMat(imol,i,1),ZSolMat(imol,i,2))
            if (iang.gt.0) then
              iparam = SolAngle(imol,iang,4)
              if (iparam.gt.0) then
                ZSolVal(imol,i,2) = AngParam(iparam,2)
              else if (iparam.lt.0) then
                call prettyLine(WARNING,"Variable ANGLES")
              endif
            endif
          enddo
        endif

c       rebuild the solute coordinates                
        call buildSolute(imol)
        
c       sync the FEP solutes...
c       call syncFEPSolutes()        
      endif
      
      return
      end
