      double precision function sumSoluteCharge(isol)
      include 'dimensions.inc'
      include 'enums.inc'
      include 'solutes.inc'
      include 'constants.inc'
      include 'printer.inc'
      integer isol,i,j
      integer lgth,lineLength
      double precision ichg(MAXSOLUTERESIDUES),charge
      
      sumSoluteCharge = ZERO
      
      if (isol.le.0 .or. isol.gt.NSolutes) return
      if (NSolResidues(isol).le.0) return
      do i=1,NSolResidues(isol)
        ichg(i) = ZERO
        do j=SolResidues(isol,i,1),SolResidues(isol,i,2)
          ichg(i) = ichg(i) + charge(ZSolPar(isol,j))  
          write(printstring,10) SolAtomNam(isol,j),SolResNam(isol,i),charge(ZSolPar(isol,j))
10        format("Charge on ",a,"-",a," = ",f12.5)
          call printLine(DETAIL,printstring)
        enddo
        write(printstring,11) SolResNam(isol,i),ichg(i)
11      format("Total charge on residue ",a," = ",f12.5)        
        call printLine(DETAIL,printstring)
        
        sumSoluteCharge = sumSoluteCharge + ichg(i)
      enddo

      if (abs(sumSoluteCharge).lt.0.00001) sumSoluteCharge = ZERO

      lgth = lineLength(SoluteName(isol))
      write(printstring,12) SoluteName(isol)(1:lgth),sumSoluteCharge
12    format("Charge on solute ",a," = ",f12.5)        
      call prettyLine(INFO,printstring)
      
      return
      end
