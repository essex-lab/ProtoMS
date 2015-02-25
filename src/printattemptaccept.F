      subroutine printAttemptAccept
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'enums.inc'
      include 'moves.inc'
      include 'pressure.inc'
c#############################################################
c
c     This prints out the attempted and accepted moves for
c     each part of the system!
c
c     (C) Christopher Woods     19 August 2002
c
c#############################################################

      integer i,j,itry,iacp,igcs,isol
      integer lgth,lineLength
      double precision iprcnt

      if (NEwaldTry.gt.0) then
        write(printstring,*) "Attempted ",NEwaldTry," Ewald moves. Accepted ",NEwaldAcp," (",
     .                      100.0*dble(NEwaldAcp)/dble(NEwaldTry),"%)"
        call prettyLine(ACCPT,printstring)
      endif

      if (NVolTry.gt.0) then
        write(printstring,*) "Attempted ",NVolTry," volume moves. Accepted ",NVolAcp," (",
     .                      100.0*dble(NVolAcp)/dble(NVolTry),"%)"
        call printLine(ACCPT,printstring)
      endif
      
      if (NProteins.gt.0) then
        do i=1,NProteins
          lgth = lineLength(ProteinName(i))
          write(printstring,*) "Attempt/accept/% moves for protein ",i," ",ProteinName(i)(1:lgth)
          call prettyLine(ACCPT,printstring)
          itry = 0
          iacp = 0
          do j=1,NResidues(i)
            if (NResTry(i,j).gt.0) then
              iprcnt = 100.0 * dble(NResAcp(i,j))/dble(NResTry(i,j))
              itry = itry + NResTry(i,j)
              iacp = iacp + NResAcp(i,j)
              write(printstring,10) j,ResNam(i,j),NResAcp(i,j),NResTry(i,j),iprcnt
10            format('Residue ',i4,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')            
              call printLine(ACCPT,printstring)
            endif
          enddo

          if (itry.gt.0) then
            iprcnt = 100.0 * dble(iacp)/dble(itry)
            write(printstring,11) iacp,itry,iprcnt     
            call printLine(ACCPT,printstring)
11          format('Total ',i8,'/',i8,'(',f8.3,'%)')
          endif
        enddo
      endif
      
      if (NSolutes.gt.0) then
        do i=1,NSolutes
          if (.not.IgnoreSolute(i)) then
            lgth = lineLength(SoluteName(i))
            write(printstring,*) "Attempt/accept/% moves for solute ",i," ",SoluteName(i)(1:lgth)
            call prettyLine(ACCPT,printstring)     
            itry = 0
            iacp = 0
            do j=1,NSolResidues(i)
              if (NSolTry(i,j).gt.0) then
                iprcnt = 100.0 * dble(NSolAcp(i,j))/dble(NSolTry(i,j))
                itry = itry + NSolTry(i,j)
                iacp = iacp + NSolAcp(i,j)
                write(printstring,20) j,SolResNam(i,j),NSolAcp(i,j),NSolTry(i,j),iprcnt
20              format('Residue ',i4,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')
                call printLine(ACCPT,printstring)
              endif
            enddo
            if (itry.gt.0) then
              iprcnt = 100.0 * dble(iacp)/dble(itry)
              write(printstring,21) iacp,itry,iprcnt
              call printLine(ACCPT,printstring)
21            format('Total ',i8,'/',i8,'(',f8.3,'%)')
            endif
          endif
        enddo
      endif
      
      if (NSolvents.gt.0) then
        write(printstring,*) "Attempt/accept/% moves for the solvent"
        call prettyLine(ACCPT,printstring)
        itry = 0
        iacp = 0
        do i=1,NSolvents   
          if (NSvnTry(i).gt.0) then
            iprcnt = 100.0 * dble(NSvnAcp(i))/dble(NSvnTry(i))
            itry = itry + NSvnTry(i)
            iacp = iacp + NSvnAcp(i)
            write(printstring,30) i,SvnName(i),NSvnAcp(i),NSvnTry(i),iprcnt
30          format('Solvent ',i5,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')        
            call printLine(ACCPT,printstring)
          endif
        enddo
        if (itry.gt.0) then
          iprcnt = 100.0 * dble(iacp)/dble(itry)
          write(printstring,31) iacp,itry,iprcnt
          call printLine(ACCPT,printstring)
31        format('Total ',i8,'/',i8,'(',f8.3,'%)')
        endif
      endif
      
      if (NGCSolutes.gt.0) then
        write(printstring,*) "Attempt/accept/% moves for insertions"
       call prettyLine(ACCPT,printstring)
        itry = 0
        iacp = 0
      do igcs=1,NGCSolutes
        if (NInsTry(igcs).gt.0) then
        iprcnt = 100.0 * dble(NInsAcp(igcs))/dble(NInsTry(igcs))
        itry = itry + NInsTry(igcs)
        iacp = iacp + NInsAcp(igcs)
       write(printstring,40) igcs,GCSName(igcs),NInsAcp(igcs),NInsTry(igcs),iprcnt
c        PRINT*, "iprcnt" ,iprcnt
40          format('Insertn ',i5,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')
           
            call printLine(ACCPT,printstring)
        else
        endif
      enddo
        if (itry.gt.0) then
          iprcnt = 100.0 * dble(iacp)/dble(itry)
           write(printstring,41) iacp,itry,iprcnt
          call printLine(ACCPT,printstring)
c           PRINT*, "%" ,iprcnt
41        format('Total ',i8,'/',i8,'(',f8.3,'%)')
        endif
      endif 

      if (NGCSolutes.gt.0) then
        write(printstring,*) "Attempt/accept/% moves for deletions"
       call prettyLine(ACCPT,printstring)
        itry = 0
        iacp = 0
      do igcs=1,NGCSolutes
        if (NDelTry(igcs).gt.0) then
        iprcnt = 100.0 * dble(NDelAcp(igcs))/dble(NDelTry(igcs))
        itry = itry + NDelTry(igcs)
        iacp = iacp + NDelAcp(igcs)
       write(printstring,50) igcs,GCSName(igcs),NDelAcp(igcs),NDelTry(igcs),iprcnt
50          format('Deletin ',i5,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')        
            call printLine(ACCPT,printstring)
        else
        endif
      enddo
        if (itry.gt.0) then
          iprcnt = 100.0 * dble(iacp)/dble(itry)
           write(printstring,51) iacp,itry,iprcnt
          call printLine(ACCPT,printstring)
51        format('Total ',i8,'/',i8,'(',f8.3,'%)')
        endif
      endif

      if (NThetaSolKind.gt.0) then
        write(printstring,*) "Attempt/accept/% moves for solute theta moves"
        call prettyLine(ACCPT,printstring)
        itry = 0
        iacp = 0
        do isol=1,NThetaSolutes
          if (NSolTheTry(isol).gt.0) then
            iprcnt = 100.0 * dble(NSolTheAcp(isol))/dble(NSolTheTry(isol))
            itry = itry + NSolTheTry(isol)
            iacp = iacp + NSolTheAcp(isol)
            write(printstring,52)'Theta ', isol,SoluteName(isol)(1:lineLength(SoluteName(isol))),NSolTheAcp(isol),NSolTheTry(isol),iprcnt  
            call printLine(ACCPT,printstring)
          endif
        enddo
        if (itry.gt.0) then
          iprcnt = 100.0 * dble(iacp)/dble(itry)
          write(printstring,51) iacp,itry,iprcnt
          call printLine(ACCPT,printstring)
        endif
      endif
52          format(a,i5,' ',a,' ',i8,'/',i8,'(',f8.3,'%)')   

      return
      end
