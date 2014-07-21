      
      subroutine simulation
      include 'dimensions.inc'
      include 'printer.inc'
      include 'flags.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'simparams.inc'
      include 'averages.inc'
      include 'rotamer.inc'
      include 'enums.inc'
      include 'surface.inc'
      include 'keyset.inc'
c###########################################################
c
c     This routine reads in the simulation description
c     and acts accordingly. If check is true then the
c     input is only checked (dry run), and no actual
c     simulation is performed
c
c     (C) Christopher Woods, October 2004
c
c###########################################################

      logical getOption,printdummies,showhidden,printstd
      integer getInt,lineLength,lgth,openFILE
      integer isol,iprnt,iprot,ires
      double precision retideltau,retideltal,getDouble,printsolvent
      double precision delta,deltaf,deltab,lb,lf
      character*64 tempstring
      integer i,j,nstps,ilgth,istrt,iend,k
      integer ichunk,numb
      logical check,getLogical,append
      integer nsplit
      double precision rad2Deg,deg2Rad,readDouble
      double precision oldscale1,oldscale2,oldscale3
c     see if this should be a dry-run
      if (getOption("dryrun")) then
        check = getLogical(words(1),.false.)
      else
        check = .false.
      endif
            
      append = .FALSE.
      numb=0

      if (PdbParams) then
c     use the fixresidues and fixbackbone lines from the protein pdb
        do k=1,nkeys
          if (keys(k)(1:9).eq.'pdbparams') then
            call fixresidues(values(k),check)
          endif
        enddo
      endif

c     read in each simulation chunk via a do loop
      ichunk = 0
      do k=1,nkeys
        if (keys(k)(1:5) .eq. 'chunk') then
          call split(values(k))
          ichunk = ichunk+1
          write(printstring,*) "Running chunk ",ichunk
          call printLine(INFO,printstring)
          
c         now see what we will do based on the first word...
          call lowercase(words(1))
          
          if (words(1).eq.'equilibrate') then
c           the format is;
c           'equilibrate nsteps (newprob) (protein=1) (solute=2) (solvent=3) (printinfo=10)'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret equilibration line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif
            
c           get the move probabilities and print move from this line
            call getMoveInfo()
            
            nstps = getInt(words(2),0)
            write(printstring,*) "Requested ",nstps," steps of equilibration"
            call prettyLine(INFO,printstring)

c           actually run the equilibration
            if (.not.check) call equilibrate(nstps)
          
          else if (words(1).eq.'simulate') then
c           the format is;
c           'simulate nsteps (newprob) (protein=1) (solute=2) (solvent=3) (printinfo=10)'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not intepret simulation line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif
            
c           get the move probabilities etc.
            call getMoveInfo()
            
            nstps = getInt(words(2),0)
            write(printstring,*) "Requested ",nstps," steps of simulation"
            call prettyLine(INFO,printstring)
            
c           actually perform the simulation
            if (.not.check) call simulate(nstps)            
 
          else if (words(1) .eq. 'printparameters') then
            call printParameters()

          else if (words(1) .eq. 'splitdddsasimulate') then
c           JM format is splitsimulate nstps split (probs)
            if (nwords.lt.3) then
              call printLine(WARNING,"Could not intepret simulation line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif            
c           get the move probabilities etc.
            call getMoveInfo()
            
            nstps = getInt(words(2),0)
            nsplit = getInt(words(3),0)
            write(printstring,*) "Requested ",nstps,
     .        " steps of splitDDDGBSAsimulation split= ",nsplit
            call prettyLine(INFO,printstring) 
c           actually perform the simulation
            if (.not.check) call splitDDDSAsimulate(nstps,nsplit) 
          else if (words(1) .eq. 'splitgbsasimulate') then
c           JM format is splitsimulate nstps split (probs)
            if (nwords.lt.3) then
              call printLine(WARNING,"Could not intepret simulation line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif            
c           get the move probabilities etc.
            call getMoveInfo()
            
            nstps = getInt(words(2),0)
            nsplit = getInt(words(3),0)
            write(printstring,*) "Requested ",nstps,
     .        " steps of split-simulationGBSA split= ",nsplit
            call prettyLine(INFO,printstring) 
c           actually perform the simulation
            if (.not.check) call splitGBSAsimulate(nstps,nsplit)                   
          else if (words(1).eq.'fakesim') then
            if (.not.check) call fakeSimulate()
          
          else if (words(1).eq.'singlepoint') then
c           the format is 'singlepoint'
            call prettyLine(INFO,"Requested the single point energy of the current system")
            if (.not.check) call singlePointEnergy
          
          else if (words(1).eq.'soluteenergy') then
c           the format is 'soluteenergy isol'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret soluteenergy line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif

c           get the solute and calculate its energy in detail            
            isol = getInt(words(2),0)
            if (isol.le.0 .or. isol.gt.NSolutes) then
              write(printstring,*) "Cannot request detailed energy of invalid solute ",isol
              call printLine(WARNING,printstring)
              cycle
            endif
            
            lgth = lineLength(SoluteName(isol))
            
            write(printstring,*) "Requested detailed energy of solute ",isol," ",SoluteName(isol)(1:lgth)
            call prettyLine(INFO,printstring)
            if (.not.check) call detailedSoluteEnergy(isol)            
          
          else if (words(1).eq.'residueenergy') then
c           the format is 'residueenergy iprot ires'
            if (nwords.lt.3) then
              call printLine(WARNING,"Could not interpret residueenergy line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif
            
            iprot = getInt(words(2),0)
            ires = getInt(words(3),0)
            if (iprot.le.0 .or. iprot.gt.NProteins .or. ires.le.0 .or. ires.gt.NResidues(iprot)) then
              write(printstring,*) "Cannot request detailed energy of invalid residue ",iprot," ",ires
              call printLine(WARNING,printstring)
              cycle
            endif
            
            lgth = lineLength(ProteinName(iprot))
            write(printstring,*) "Requested detailed energy of residue ",ResNam(iprot,ires),
     .                           "-",ResNum(iprot,ires)," in protein ",ProteinName(iprot)(1:lgth)
            call printLine(INFO,printstring)
            if (.not.check) call detailedResidueEnergy(iprot,ires)       

c     SG REST Scaling, not in official release
c          else if (words(1).eq.'rest2') then
c           the format is 'rest2'
c            call prettyLine(INFO,"Requested the REST2 energy of the current system")
c            oldscale1 = restscale1
c            oldscale2 = restscale2
c            oldscale3 = restscale3
c            restscale1 = 1.0d0
c            restscale2 = 0.50d0
c            restscale3 = 0.0d0
c            if (.not.check) call singlePointEnergy
c            restscale1 = oldscale1
c            restscale2 = oldscale2
c            restscale3 = oldscale3        
     
          else if (words(1).eq.'retienergy') then
c           the format is 'retienergy retidelta'
c           JM mod format can also be retienergy retideltau retideltal
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not intepret retienergy line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif

            retideltau = getDouble(words(2),ZERO)
            if (nwords .gt. 2) then 
              retideltal = getDouble(words(3),ZERO)
            else
              retideltal = retideltau
            endif
            write(printstring,*) "Requested retienergy using retidelta=",retideltau,retideltal
            call prettyLine(INFO,printstring)
            
            if (.not.check) call retiEnergy(retideltau,retideltal)

          else if (words(1).eq.'freeenergy') then
            call join(1,nwords,tempstring)
            call dump(tempstring,append,numb)

          else if (words(1) .eq. 'resetgb') then
            if (.not. check) call resetgbenergy()

          else if (words(1).eq.'averages' .or. words(1).eq.'results') then
            call join(1,nwords,tempstring)
            call dump(tempstring,append,numb)
 
          else if (words(1).eq.'rotamer') then
c           call to read in a rotamer library, or change a rotamer library parameter

c           format is 'rotamer read filename'
c           or 'rotamer set param1=value1 param2=value2'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret rotamer line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif

c           temporarily convert back the current deltas to degrees
            BigMoveDelta = rad2Deg(BigMoveDelta)
            RotamerDelta = rad2Deg(RotamerDelta)

c           we are setting rotamer parameters...
            do j=2,nwords
              call lowercase(words(j))
              ilgth = lineLength(words(j))

              if (words(j)(1:7).eq.'bigprob') then
                BigMoveProb = readDouble( words(j)(9:ilgth) ) / 100.0D+00
                if (BigMoveProb .gt. ONE) BigMoveProb = ONE
                if (BigMoveProb .lt. ZERO) BigMoveProb = ZERO
              else if (words(j)(1:8).eq.'bigdelta') then
                BigMoveDelta = readDouble( words(j)(10:ilgth) )
              else if (words(j)(1:8).eq.'rotdelta') then
                RotamerDelta = readDouble( words(j)(10:ilgth) )
              endif
            enddo

c           now see if we should read in the rotamer library
            if (words(2).eq.'read') then
              if (nwords.lt.3) then
                call printLine(WARNING,"Need to specify a filename containing the rotamer file!")
                cycle
              endif

              ilgth = lineLength(words(3))
              write(printstring,*) "Reading rotamer library file ",words(3)(1:ilgth)
              call prettyLine(INFO,printstring)
              if (check) cycle
     
c             actually read in the rotamer library file
              call readRotamer(words(3)(1:ilgth))
            endif

            write(printstring,*) "BigMoveProb = ",BigMoveProb*100.0,"%, BigMoveDelta = ",
     .           BigMoveDelta," degrees, RotamerDelta = ",RotamerDelta," degrees"
            call printLine(INFO,printstring)

c           convert the deltas to radians
            BigMoveDelta = deg2Rad(BigMoveDelta)
            RotamerDelta = deg2Rad(RotamerDelta)
         
          else if (words(1).eq.'restart') then
            call join(1,nwords,tempstring)
            call dump(tempstring,append,numb)

          else if (words(1).eq.'mcpro') then
c           format is;
c           'mcpro solute=x|protein=x file=filename  #comment'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret mcpro line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif

            call prettyLine(INFO,"Writing MCPro z-matrix and parameter files...")

c           scan the words to see if the file has been specified
            if (nwords.gt.2) then
              do j=3,nwords
                if (words(j)(1:1).eq.'#') cycle

                lgth = lineLength(words(j))
                if (words(j)(1:4).eq.'file') then
c                 change the pdb file! - connect the PDB stream to a new filehandle
                  write(printstring,*) "Redirecting PDB to ",words(j)(6:lgth)
                  call printLine(INFO,printstring)
                  streamhandle(MCPRO) = openFILE(words(j)(6:lgth))
                endif
              enddo          
            endif

            lgth = lineLength(words(2))

            if (words(2)(1:7).eq.'protein') then
c             get which protein will be printed
              iprnt = getInt(words(2)(9:lgth),0)
              if (check) then
                write(printstring,*) "Printing protein ",iprnt
                call printLine(INFO,printstring)
              endif
              if (.not.check) call WriteMCPro(PROTEIN,iprnt)
            else if (words(2)(1:6).eq.'solute') then
c             get which solute will be printed
              iprnt = getInt(words(2)(8:lgth),0)
              if (check) then
                write(printstring,*) "Printing solute ",iprnt
                call printLine(INFO,printstring)
              endif
              if (.not.check) call WriteMCPro(SOLUTE,iprnt)
            else
              call printLine(WARNING,"Could not work out what to print?")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif

          else if (words(1).eq.'pdb') then
            call join(1,nwords,tempstring)
            call dump(tempstring,append,numb)

          else if (words(1).eq.'fixresidues'.or.words(1).eq.'fixbackbone') then
            call join(1,nwords,tempstring)
            call fixresidues(tempstring,check)

          else if (words(1).eq.'transrot') then
            if (nwords.ge.3) then
              i = getInt(words(2),1)
              SolADel(i) = getDouble(words(3),SolADel(i))
            endif
            if (nwords.ge.4) SolRDel(i) = deg2Rad(getDouble(words(4),Rad2deg(SolRDel(i))))
            if (nwords.ge.3) then 
              write(printstring,*)"New displacements for solute ",i," translation=",SolADel(i)," rotation=",SolRDel(i)
              call printLine(INFO,printstring)
            endif

          else if (words(1).eq.'geometry') then
c           change the geometry of a molecule
            call processGeometryLine(check)
      
          else if (words(1).eq.'solvent') then
c           commands to build solvent boxes / caps
            call processSolventLine(check)
                  
          else if (words(1).eq.'lambda') then
c           change lambda - calculate change in energy and print that as well
c           this could be used for slow/fast growth. Format is;
c           'lambda newlam (lamf lamb)'  or 'lambda delta newlam'          
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret lambda line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif
            
            call lowercase(words(2))
            
            if (words(2).eq.'delta') then
c             change lambda (lambdab,lambdaf) by delta
              if (nwords.lt.3) then
                call printLine(WARNING,"Must specify how much you want to change lambda!")
                call join(1,nwords,printstring)
                call printLine(WARNING,printstring)
                cycle
              endif
              
              delta = getDouble(words(3),ZERO)
              
              lambda = lambda + delta
              lambdaf = lambdaf + delta
              lambdab = lambdab + delta
            else
c             try to read the new value of lambda
              if (nwords.ge.4) then
c               get lambda, lambdaf and lambdab
                lambda = getDouble(words(2),ZERO)
                lambdaf = getDouble(words(3),ZERO)
                lambdab = getDouble(words(4),ZERO)
              else
                if (DifferentLambdas) then
c                 just get lambda - preserve the deltas to lambdaf and lambdab
                  deltaf = lambdaf - lambda
                  deltab = lambdab - lambda
                
                  lambda = getDouble(words(2),ZERO)
                  lambdaf = lambda + deltaf
                  lambdab = lambdab + deltab                            
                else
                  lambda = getDouble(words(2),ZERO)
                  lambdaf = lambda
                  lambdab = lambda
                endif
              endif
            endif
            
c           sanitise the values of lambda
            if (lambda.lt.ZERO) lambda = ZERO
            if (lambda.gt.ONE) lambda = ONE
            
            if (lambdab.lt.ZERO) lambdab = ZERO
            if (lambdab.gt.ONE) lambdab = ONE
            
            if (lambdaf.lt.ZERO) lambdaf = ZERO
            if (lambdaf.gt.ONE) lambdaf = ONE
                
            write(printstring,*) "Changing lambda to ",lambda,"(",lambdab,",",lambdaf,")"
            call printLine(INFO,printstring)
            
            if (.not.check) call changeLambda()
            
          else if (words(1).eq.'setstream') then
c           change the files associated with the streams! format is;
c           'setstream info='filename' fatal='filename' etc.'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret setstream line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              cycle
            endif
            
            do i=2,nwords
c             split the word by '=' sign
              lgth = lineLength(words(i))
c             scan the letters to find the '=' sign
              do j=1,lgth
                if (words(i)(j:j).eq.'=') then
                  if (j.eq.1 .or. j.eq.lgth) then
                    call printLine(WARNING,"Could not interpret setstream line")
                    call join(1,nwords,printstring)
                    call printLine(WARNING,printstring)
                    cycle
                  endif
                  
c                 set the stream file
                  if (check) then
                    write(printstring,*) "Directing stream ",words(i)(1:j-1)," to file ",
     .                                   words(i)(j+1:lgth)
                    call printLine(INFO,printstring)
                    goto 200
                  else
                    call setStream(words(i)(1:j-1),words(i)(j+1:lgth))
                    goto 200
                  endif
                endif
              enddo
200           continue
            enddo
          else if (words(1) .eq. 'surface') then
            if (doSA) then
              call calculateTotalSurface(.true.)
            else
              call printLine(WARNING,"Cannot calculate surface as surface area calculations were not specified !")
              call printLine(WARNING,printstring)
            endif
c           Deal with setid/clear id lines
          else if (words(1) .eq. 'id') then
            call processIdline(check)
c           Deal with the adding of restrains
          else if (words(1) .eq. 'restraint') then
            call processRestraintLine(check)
          else if (words(1) .eq. 'hardwall') then
           call processHardwallLine(check)
          else if (words(1) .eq. 'addclj') then
            call addSoluteCLJ(getInt(words(2),ZERO),getInt(words(3),ZERO),getDouble(words(4),ZERO))
          else
            call printLine(WARNING,"Cannot intepret simulation chunk!")
            call join(1,nwords,printstring)
            call printLine(WARNING,printstring)
          endif
          cycle
        endif
c     end of do-while loop
      enddo

      return
      
      end


