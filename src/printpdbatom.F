      subroutine printPDBAtom(anum,anam,rnam,rnum,ix,iy,iz,ipar,showdummies,printstd)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'parameters.inc'
c####################################################
c
c     This subroutine prints the atom given in the arguments 
c     (atomnum,atomname,resname,resnum,x,y,z,cljparameter)
c     as a pdb atom, to iunit
c
c     (C) Christopher Woods  21 January 2002
c
c####################################################
      
      integer anum,rnum,ipar,ipar1,ipar2
      integer chkres,chkatm
      character*4 anam,rnam
      double precision ix,iy,iz
      double precision charge
      logical showdummies,printstd
      integer cljID

      if (ipar.eq.0) then
c       skip the atom if we are not printing dummies        
        if (.not.showdummies) return
        ipar1 = 0
        ipar2 = 0
      else if (ipar.lt.0) then
        ipar1 = cljID(VariableCLJ(-ipar,1))
        ipar2 = cljID(VariableCLJ(-ipar,2))
      else
        ipar1 = cljID(ipar)
        ipar2 = ipar1
      endif

c     ensure that the residue and atom numbers are bounded between
c     1 and 9999 (and loop between these values)      
      chkres = rnum
100   continue
      if (chkres.gt.9999) then
        chkres = chkres - 10000
        goto 100
      else if (chkres.le.0) then
        chkres = chkres + 9999
        goto 100
      endif
      
      chkatm = anum
101   continue
      if (chkatm.gt.9999) then
        chkatm = chkatm - 10000
        goto 101
      else if (chkatm.le.0) then
        chkatm = chkatm + 9999
        goto 101
      endif
      

      if (printstd) then
        call uppercase(anam)
        call uppercase(rnam)
        write(printstring,11) chkatm,anam,rnam,chkres,ix,iy,iz
            
        call printLine(PDB,printstring)
11      format("ATOM",3x,i4,1x,a4,1x,a4,1x,i4,4x,3f8.3)
      else
        write(printstring,10) chkatm,anam,rnam,chkres,ix,iy,iz,ipar1,ipar2,charge(ipar)
            
        call printLine(PDB,printstring)
10      format("ATOM",3x,i4,1x,a4,1x,a4,1x,i4,4x,3f8.3,2x,i5,1x,i5,1x,f8.3)
      endif

      return 
      end
