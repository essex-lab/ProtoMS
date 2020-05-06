      subroutine printPDB(ityp,ival,svndist,showdummies,showhidden,printstd,numb)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
c#####################################################################
c
c     This routine prints out the current system in PDB format.
c     
c     'ityp' controls what is printed;
c 
c     ityp = 0 : all proteins and solutes
c     ityp = PROTEIN : Protein ival is printed (all if ival <= 0)
c     ityp = SOLUTE : Solute ival is printed (all if ival <= 0)
c
c     svndist controls how much of the solvent is printed.
c     If svndist is < 0 then all of the solvent is printed
c     If svndist == ZERO then none of the solvent is printed
c     otherwise only solvents with svndist A of any printed
c     protein or solute will be printed.
c
c     Dummy atoms are printed if printdummies is true
c     Hidden solutes (ignored solutes) are printed if
c     showhidden is true
c
c     (C) Christopher Woods, October 2004
c
c####################################################################

      integer ityp,ival
      double precision svndist
      logical showdummies,showhidden,printstd
      double precision minSolventProteinDist2,minSolventSoluteDist2
      double precision minSolventGCSoluteDist2
      integer i,j,k,numb
      
      if (NSolvents.gt.0) then
        if (svndist.eq.ZERO) then
          call printLine(PDB,"REMARK No solvent molecules are shown")
        else if (svndist.gt.ZERO) then
          write(printstring,*) "REMARK Only solvent molecules within ",svndist,
     .                         " A of the system are shown"
          call leftJustify(printstring)
          call printLine(PDB,printstring)
        endif
        
        if (svndist.ne.ZERO) then            
c         write information about the boundary conditions for the solvent
          if (SolventBox) then
c           write the box dimensions to the file            
            write(printstring,"('HEADER ',a3,'  ',6f12.4)") 
     &            "box",(BoxOrigin(i),i=1,3),(BoxTop(i),i=1,3)
            call printLine(PDB,printstring)
          else if (SolventCap) then
c           write out the solvent cap information to the file
            write(printstring,"('HEADER ',a3,'  ',5f12.4)") 
     &            "cap",(cap(i),i=1,3),capradius,capconstant
            call printLine(PDB,printstring)
          endif
        endif
      endif
      
      if (ityp.eq.PROTEIN .or. ityp.le.0) then
        if (ival.le.0 .and. NProteins.gt.0) then
          do i=1,NProteins
            call printPDBProtein(i,showdummies,printstd)
          enddo
        else if (ival.gt.0 .and. ival.le.NProteins) then
          call printPDBProtein(ival,showdummies,printstd)
        endif
      endif
      
      if (ityp.eq.SOLUTE .or. ityp.le.0) then
        if (ival.le.0 .and. NSolutes.gt.0) then
          do i=1,NSolutes
            if (showhidden) then
              call printPDBSolute(i,showdummies,printstd)
            else if ((.not.ignoreSolute(i)).and.(SoluteTheta(i).ge.thres))  then
              call printPDBSolute(i,showdummies,printstd)
            endif
          enddo
        else if (ival.gt.0 .and. ival.le.NSolutes) then
          call printPDBSolute(i,showdummies,printstd)
        endif
      endif
          
       if (ityp.eq.GCSOLUTE .or. ityp.le.0) then
        if (ival.le.0 .and. NGCSolutes.gt.0) then
          do i=1,NGCSolutes
             call printPDBGCSolute(i,showdummies,printstd)
          enddo
        endif
       endif

c     now do the solvent...
      if (NSolvents.gt.0) then
        if (svndist.lt.ZERO) then
c         print all of the solvent      
          do k=1,NSolvents
            call printPDBSolvent(k,showdummies,printstd)
          enddo
        else if (svndist.gt.ZERO) then
c         print all solvents within svndist of whatever has been printed

c         square the distance to remove sqrts
          svndist = svndist**2

          do i=1,NSolvents
            if (ityp.eq.PROTEIN .or. ityp.le.0) then
              if (ival.le.0 .and. NProteins.gt.0) then
                do j=1,NProteins
                  if (minSolventProteinDist2(j,i).le.svndist) then
                    call printPDBSolvent(i,showdummies,printstd)
                    goto 100
                  endif
                enddo
              else if (ival.gt.0 .and. ival.le.NProteins) then
                if (minSolventProteinDist2(ival,i).le.svndist) then
                  call printPDBSolvent(i,showdummies,printstd)
                  goto 100
                endif
              endif
            endif
            
            if (ityp.eq.SOLUTE .or. ityp.le.0) then
              if (ival.le.0 .and. NSolutes.gt.0) then
                do j=1,NSolutes
                  if (minSolventSoluteDist2(i,j).le.svndist) then
                    call printPDBSolvent(i,showdummies,printstd)
                    goto 100
                  endif
                enddo
              else if (ival.gt.0 .and. ival.le.NSolutes) then
                if (minSolventSoluteDist2(i,j).le.svndist) then
                  call printPDBSolvent(i,showdummies,printstd)
                  goto 100
                endif
              endif
            endif

           if (ityp.eq.GCSOLUTE .or. ityp.le.0) then
              if (ival.le.0 .and. NGCSolutes.gt.0) then
                do j=1,NGCSolutes
                  if (minSolventGCSoluteDist2(i,j).le.svndist) then
                    call printPDBSolvent(i,showdummies,printstd)
                    goto 100
                  endif
                enddo
              else if (ival.gt.0 .and. ival.le.NGCSolutes) then
                if (minSolventGCSoluteDist2(i,j).le.svndist) then
                  call printPDBSolvent(i,showdummies,printstd)
                  goto 100
                endif
              endif
            endif
 
100         continue
          enddo

        endif
      endif
        
c     finish off the pdb file with an 'END' line
      if (numb.eq.0) then
        call printLine(PDB,"END")
      endif

      return
      end
