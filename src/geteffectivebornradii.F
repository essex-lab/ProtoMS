      subroutine GetEffectiveBornRadii
c##################################################
c Calculate the effective Born radius of every atom
c in the system. Must be called once the pairwise 
c Contributions have been calculated
c
c     Julien Michel November 2004
c##################################################
      include 'dimensions.inc'
      include 'constants.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'gb.inc'
      integer i,j,k,nres,resstart,resend,ipar
      double precision borni,bornif,bornib
      double precision psi

      if (NProteins .gt. ZERO) then
        do 100 i=1,NProteins
          nres = NResidues(i)
          do 110 j=1,nres
            resstart = BBZid(i,j,1)
            resend = BBZid(i,j,4)
            do 120 k=resstart,resend
              ipar = ZPar(i,k)
              if (ipar .lt. 0) then
                borni = 0.0d0
                bornif = borni
                bornib = borni
                write(printstring,*) 'No support for perturbed ',
     &            'protein atoms for GB'
                call printLine(FATAL,printstring)
                call closeProgram(-1)
              else
                borni = IntrBornDB(CLJType(ipar))
                bornif = borni
                bornib = borni
              endif
c             write(printstring,*) 'prot ',i,' res ',j,' at ',k,
c             &        'cont tot',ZProtContBorn(i,k),
c             &        'new ',NewZProtContBorn(i,k),' old ',OldZProtContBorn(i,k)
c             call printLine(DEBUG,printstring)
              ZProtContBorn(i,k) = ZProtContBorn(i,k) + 
     &          NewZProtContBorn(i,k) - OldZProtContBorn(i,k)
              ZProtContBornF(i,k) = ZProtContBornF(i,k) + 
     &          NewZProtContBornF(i,k) - OldZProtContBornF(i,k)
              ZProtContBornB(i,k) = ZProtContBornB(i,k) + 
     &          NewZProtContBornB(i,k) - OldZProtContBornB(i,k)
c             Consider applying mapping from Onufriev 
              if (gbproteins) then
                psi = 0.5*borni*ZProtContBorn(i,k)
                ZProtEffBorn(i,k) = 1 / ( (1 / borni) - 
     &            tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &            / (borni + 0.09) )

                psi = 0.5*bornif*ZProtContBornF(i,k)
                ZProtEffBornF(i,k) = 1 / ( (1 / bornif) - 
     &            tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &            / (bornif + 0.09) )

                psi = 0.5*bornib*ZProtContBornB(i,k)
                ZProtEffBornB(i,k) = 1 / ( (1 / bornib) - 
     &            tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &            / (bornib + 0.09) )
              else
c               Default Method (avoids negative born radii)
                ZProtEffBorn(i,k) = DMAX1( 1 / ( 1 / borni - 
     &            0.5*ZProtContBorn(i,k) ), borni )

                ZProtEffBornF(i,k) = DMAX1( 1 / ( 1 / bornif - 
     &            0.5*ZProtContBornF(i,k) ), bornif )

                ZProtEffBornB(i,k) = DMAX1( 1 / ( 1 / bornib - 
     &            0.5*ZProtContBornB(i,k) ), bornib )
            endif
c             Comment out to inspect contents of arrays
c              write(printstring,101) i,j,k,borni,
c     &          ZProtEffBorn(i,k),0.5*ZProtContBorn(i,k)
c 101          format ('PROTEIN',i2,' RESIDUE ',i4,' ATOM',i5,
c     &          ' Born init ',f6.3,' Born eff ',f6.3,' Cont ',f12.8)
c              call printLine(DEBUG,printstring)
c              write(printstring,*) 'BEff Fwd ', ZProtEffBornF(i,k),
c     &          ' BEff Bwd ',ZProtEffBornB(i,k)
c              call printLine(DEBUG,printstring) 
 120        continue
 110      continue
 100    continue
      endif

      if (NSolutes .gt. ZERO) then
        do 200 i=1,NSOLUTES
          if (IgnoreSolute(i)) goto 200
          do 210 j=1,NSolAtoms(i)
c            borni = ZSolIBorn(i,j)
            ipar = ZSolPar(i,j)
            if (ipar .lt. 0) then
              borni = VariableIntrBorn(-ipar,2)
c             What to do with FEP ?
              if (DifferentLambdas) then
                bornif = VariableIntrBorn(-ipar,5)
                bornib = VariableIntrBorn(-ipar,4)
              else
                bornif = borni
                bornib = borni
              endif
            else
              borni = IntrBornDB(CLJType(ipar))
              bornif = borni
              bornib = borni
            endif
c           write(printstring,*) 'sol ',i,' at ',j,'cont tot',ZSolContBorn(i,j),
c           &      'new ',NewZSolContBorn(i,j),' old ',OldZSolContBorn(i,j)
c           call printLine(DEBUG,printstring)
c            write(printstring,*) 'At ',j,' Cont was ',ZSolContBorn(i,j),
c     &        ' New is ',NewZSolContBorn(i,j),' Old is ',OldZSolContBorn(i,j)
c            call printLine(DEBUG,printstring)
            ZSolContBorn(i,j) = ZSolContBorn(i,j) + 
     &        NewZSolContBorn(i,j) - OldZSolContBorn(i,j)
c            write(printstring,*) 'At ',j,' Cont is now ',ZSolContBorn(i,j)
c            call printLine(DEBUG,printstring)
c           Need a Forward contribution
            ZSolContBornF(i,j) = ZSolContBornF(i,j) + 
     &        NewZSolContBornF(i,j) - OldZSolContBornF(i,j)
c           Need a Backward contribution
c            write(printstring,*) 'At ',j,' ContB was ',ZSolContBornB(i,j),
c     &        ' NewB is ',NewZSolContBornB(i,j),' OldB is ',OldZSolContBornB(i,j)
c            call printLine(DEBUG,printstring)
            ZSolContBornB(i,j) = ZSolContBornB(i,j) + 
     &        NewZSolContBornB(i,j) - OldZSolContBornB(i,j)
c            write(printstring,*) 'At ',j,' ContB is now ',ZSolContBornB(i,j)
c            call printLine(DEBUG,printstring)
c           Consider applying mapping from Onufriev 
            if (gbproteins) then
c
c             Should I save the pre-mapping value of ZSolEffBorn to
c             apply more effectively FastGB ?
c
              psi = 0.5*borni*ZSolContBorn(i,j)
              ZSolEffBorn(i,j) = 1 / ( (1 / borni) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (borni + 0.09) )

              psi = 0.5*bornif*ZSolContBornF(i,j)
              ZSolEffBornF(i,j) = 1 / ( (1 / bornif) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (bornif + 0.09) )

              psi = 0.5*bornib*ZSolContBornB(i,j)
              ZSolEffBornB(i,j) = 1 / ( (1 / bornib) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (bornib + 0.09) )
            else
              ZSolEffBorn(i,j) = DMAX1( 1 / ( 1 / borni - 
     &          0.5*ZSolContBorn(i,j) ), borni )

              ZSolEffBornF(i,j) = DMAX1( 1 / ( 1 / bornif - 
     &          0.5*ZSolContBornF(i,j) ), bornif )

              ZSolEffBornB(i,j) = DMAX1( 1 / ( 1 / bornib - 
     &          0.5*ZSolContBornB(i,j) ), bornib )
            endif
c           Comment out to inspect contents of arrays
c            write(printstring,102) i,j,SolAtomNam(i,j),borni,
c     &      ZSolEffBorn(i,j),ZSolContBorn(i,j)
c 102        format ('SOLUTE',i2,' ATOM',i5,' Name ',a4,' Born init ',f6.3,
c     &        ' Born eff ',f6.3,' Cont ',f8.5)
c            call printLine(DEBUG,printstring)
c            write(printstring,103) bornif,ZSolEffBornF(i,j),bornib,ZSolEffBornB(i,j)
c 103        format ('----',' IBorn Fwd ',f6.3,' BornEff Fwd ',f6.3,
c     &        ' IBorn Bwd ',f6.3,' BornEff Bwd ',f6.3)
c            call printLine(DEBUG,printstring)
c            write(printstring,104) j,ZSolContBorn(i,j),ZSolContBornF(i,j),ZSolContBornB(i,j)
c 104        format('At ',i3,' Cont ',f8.5,' Cont fwd ',f8.5,' Cont Bwd ',f8.5)
c            call printLine(DEBUG,printstring)
 210      continue
 200    continue
      endif

      if (NSolvents .gt. ZERO) then
        do 300 i=1,NSolvents
          do 310 j=1,NSvnAtoms(i)
c           borni = ZSvnIBorn(i,j)
            ipar = ZSvnPar(i,j)
            if (ipar .lt. 0) then
              borni = VariableIntrBorn(-ipar,2)
c             What to do with FEP ?
              if (DifferentLambdas) then
                bornif = VariableIntrBorn(-ipar,5)
                bornib = VariableIntrBorn(-ipar,4)
              else
                bornif = borni
                bornib = borni
              endif
            else
              borni = IntrBornDB(CLJType(ipar))
              bornif = borni
              bornib = borni
            endif

            ZSvnContBorn(i,j) = ZSvnContBorn(i,j) + 
     &        NewZSvnContBorn(i,j) - OldZSvnContBorn(i,j)

            ZSvnContBornF(i,j) = ZSvnContBornF(i,j) + 
     &        NewZSvnContBornF(i,j) - OldZSvnContBornF(i,j)

            ZSvnContBornB(i,j) = ZSvnContBornB(i,j) + 
     &        NewZSvnContBornB(i,j) - OldZSvnContBornB(i,j)

c           Consider applying mapping from Onufriev 
            if (gbproteins) then
              psi = 0.5*borni*ZSvnContBorn(i,j)
              ZSvnEffBorn(i,j) = 1 / ( (1 / borni) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (borni + 0.09) )

              psi = 0.5*bornif*ZSvnContBornF(i,j)
              ZSvnEffBornF(i,j) = 1 / ( (1 / bornif) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (bornif + 0.09) )

              psi = 0.5*bornib*ZSvnContBornB(i,j)
              ZSvnEffBornB(i,j) = 1 / ( (1 / bornib) - 
     &          tanh ( (ONUA + ONUB * psi + ONUG * psi **2) * psi ) 
     &          / (bornib + 0.09) )
            else
              ZSvnEffBorn(i,j) = DMAX1( 1 / ( 1 / borni - 
     &          0.5*ZSvnContBorn(i,j) ), borni )

              ZSvnEffBornF(i,j) = DMAX1( 1 / ( 1 / bornif - 
     &          0.5*ZSvnContBornF(i,j) ), bornif )

              ZSvnEffBornB(i,j) = DMAX1( 1 / ( 1 / bornib - 
     &          0.5*ZSvnContBornB(i,j) ), bornib )
            endif
c           Comment out to inspect contents of arrays
c           write(printstring,113) i,j,SvnAtomName(i,j),borni,
c    &        ZSvnEffBorn(i,j),ZSvnContBorn(i,j)
c113        format ('SOLVENT',i2,' ATOM',i5,' Name ',a4,' Born init ',f6.3,
c    &        ' Born eff ',f6.3,' Cont ',f8.5)
c           call printLine(DEBUG,printstring)
c           write(printstring,114) bornif,ZSvnEffBornF(i,j),bornib,ZSvnEffBornB(i,j)
c114        format ('----',' IBorn Fwd ',f6.3,' BornEff Fwd ',f6.3,
c    &        ' IBorn Bwd ',f6.3,' BornEff Bwd ',f6.3)
c           call printLine(DEBUG,printstring)
c           write(printstring,115) j,ZSvnContBorn(i,j),ZSvnContBornF(i,j),ZSvnContBornB(i,j)
c115       format('At ',i3,' Cont ',f8.5,' Cont fwd ',f8.5,' Cont Bwd ',f8.5)
c           call printLine(DEBUG,printstring)
 310      continue
 300    continue
      endif

      return
      end
