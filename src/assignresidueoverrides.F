
      subroutine assignResidueOverrides(iprot,ires,iovr)
      include 'dimensions.inc'
      include 'templates.inc'
      include 'proteins.inc'
      include 'printer.inc'
c##############################################################
c
c     This function assigns override parameters to residue
c     ires in protein iprot from residue override template iovr      
c
c     (C) Christopher Woods, November 2004
c
c##############################################################

      integer iprot,ires,iovr
      integer i,iatm,findAtom
      integer ilgth,lineLength
      integer ipar0,ipar1,addCLJ,getCLJ
      
      if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0) return
      if (ires.gt.NResidues(iprot)) return
      
      if (NOverrideParams(iovr).le.0) return
      
c     loop over each parameter in this override and see if an atom of this
c     name is present in this residue
      do i=1,NOverrideParams(iovr)
c       find the atom with this name
        iatm = findAtom(iprot,ires,OverrideAtm(iovr,i))
        if (iatm.le.0) then
          ilgth = lineLength(ProteinName(iprot))
          write(printstring,10) OverrideAtm(iovr,i),ProteinName(iprot)(1:ilgth),ResNam(iprot,ires),
     .                         ResNum(iprot,ires)
10        format("Cannot find override atom ",a," in protein ",a," residue ",a," (",i5,")")
          call printLine(WARNING,printstring)
        else
c         overwrite the atom's parameter
          ipar0 = addCLJ(OverrideParam(iovr,i,1))
          ipar1 = addCLJ(OverrideParam(iovr,i,2))
          
          ZPar(iprot,iatm) = getCLJ(ipar0,ipar1)
        endif
      enddo 
      
      return
      end
