      subroutine dump(line,append,numb)

      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'averages.inc'
      include 'solutes.inc'
      include 'solvents.inc'
      include 'proteins.inc'
      include 'enums.inc'
      include 'flags.inc'
      include 'parallel.inc'

      logical getLogical,printdummies,showhidden,printstd,append,newfile,getoption
      double precision printsolvent,getDouble,lb,lf,lambda
      double precision nrg,nrgb,nrgf
      integer lgth,j,iprnt,ilgth,lineLength,openFILE,getint,i
      integer numb
      character*64 line,tempstring
      character*600 tempname,cmd
      character*255 fullfilename

      
      call split(line)

c     numb will be 0 is this is comming from a chunk
c     numb will be 1 if this is the first time we print in a file from a dump
c     numb will be more than 1 if we are appending to a file from dump

c PDB -------------------------------------------
          if (words(1).eq.'pdb') then
c           format is;
c           'pdb all|protein=x|solute=x (file=filename) (solvent=all|dist) (printdummies) (standard) #comment'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret pdb line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              goto 100
            endif

            printsolvent = ZERO
            printdummies = .false.
            showhidden = .false.
            printstd = .false.
            newfile = .false.

            if (numb.le.1) call prettyLine(INFO,"Printing PDB of the system")
            
c           scan the words to see if the file or solvent have been specified
            if (nwords.gt.2) then
              do j=3,nwords
                if (words(j)(1:1).eq.'#') goto 101
                
                lgth = lineLength(words(j))
                if (words(j)(1:4).eq.'file') then
c                 change the pdb file! - connect the PDB stream to a new filehandle
                  if (numb.le.1) then
                    write(printstring,*) "Redirecting PDB to ",words(j)(6:lgth)
                    call printLine(INFO,printstring)
                  endif
                  streamhandle(PDB) = openFILE(words(j)(6:lgth))
                  newfile = .TRUE.
                else if (words(j)(1:7).eq.'solvent' .and. NSolvents.gt.0) then
c                 change how much of the solvent is printed
                  if (words(j)(9:11).eq.'all') then
c                   print all of the solvent
                    printsolvent = -1.0
                  else
c                   print all solvent molecules within a certain distance from
c                   everything else that was printed
                    printsolvent = getDouble(words(j)(9:lgth),ZERO)
                  endif
                else if (words(j).eq.'showdummies') then
                  printdummies = .true.
                else if (words(j).eq.'showhidden') then
                  showhidden = .true.
                else if (words(j).eq.'standard') then
                  printstd = .true.
                endif
              enddo
101           continue
            endif
                
            lgth = lineLength(words(2)) 
            if (.not.newfile) then
              if (getOption('streampdb')) then
                streamhandle(PDB) = openFile(words(1))
              else
                streamhandle(PDB) = STDOUT
              endif  
              call split(line)
            endif     

            fileappend(streamhandle(PDB)-40) = (append.or.appendstreams)
c (append.or.appendstreams)

            if (numb.le.1) then
              call printLine(PDB,"HEADER Simulation output from ProtoMS")
            endif

            if (numb.ne.0) then
              if (numb.le.9999) then
                write(tempstring,"(A,I4)")"MODEL     ",numb
              else
                write(tempstring,"(A,I8)")"MODEL     ",numb
              endif
              call printline(PDB, tempstring)
            endif  
                    
            if (words(2).eq.'all') then
c             all of the proteins and solutes will be printed!            
              call printPDB(0,0,printsolvent,printdummies,showhidden,printstd,numb)
            else if (words(2)(1:7).eq.'protein') then
              if (words(2)(9:11).eq.'all') then
c               all of the proteins will be printed!
                iprnt = 0
              else
c               get which protein will be printed
                iprnt = getInt(words(2)(9:lgth),0)
              endif
              call printPDB(PROTEIN,iprnt,printsolvent,printdummies,showhidden,printstd,numb)
            else if (words(2)(1:6).eq.'solute') then
              if (words(2)(8:10).eq.'all') then
c               all of the solutes will be printed!
                iprnt = 0
              else
c               get which solute will be printed
                iprnt = getInt(words(2)(8:lgth),0)
              endif
              call printPDB(SOLUTE,iprnt,printsolvent,printdummies,showhidden,printstd,numb)
            else if (words(2)(1:7).eq.'solvent') then
              iprnt = 0
              printsolvent = -1
              call printPDB(SOLVENT,iprnt,printsolvent,printdummies,showhidden,printstd,numb)
            else
              call printLine(WARNING,"Could not work out what to print?")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif  
            if (numb.ne.0) then
              call printline(PDB,'ENDMDL')
            endif   

c FREEENERGY -------------------------------------------------
c This is done by printaverages now

c          else if (words(1).eq.'freeenergy') then
c           SG: format is 'freeenergy lambdab labmbdaf'
c            if (nwords.lt.3) then
c              call printLine(WARNING,"Could not intepret freeenergy line")
c              call join(1,nwords,printstring)
c              call printLine(WARNING,printstring)
c              goto 100
c            else
c              lb = getDouble(words(2),lambda)
c              lf = getDouble(words(3),lambda)
c            endif
c            write(printstring,*) "Requested free energy evaluation using lambdab and lambdaf=",lb,lf
c            call prettyLine(INFO,printstring) 
c            call dGEnergies(lb,lf)

c RESTART ------------------------------------------------------------

          else if (words(1).eq.'restart') then
c           the restart files will always be ovewriten
            append = .FALSE.
c           format is;
c           'restart read|write filename'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret restart line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              goto 100
            endif
            
            if (words(2).eq.'read') then
              if (nwords.lt.3) then
                call printLine(WARNING,"Need to specify filename containing restart file")
                goto 100
              endif
              
              ilgth = lineLength(words(3))
              if (numb.le.1) then
                write(printstring,*) "Reading restart file from file '",words(3)(1:ilgth),"'"
                call prettyLine(INFO,printstring)
              endif
c                    do i=1,NGCSolutes
c               PRINT*, "GCS,theta",i,thetaparam(i)
c                   enddo

              call getFullFilename(words(3),fullfilename)
              ilgth = lineLength(fullfilename)              
              call readRestart(fullfilename(1:ilgth))
            else if (words(2).eq.'write') then
              if (nwords.ge.3) then
c               redirect the stream to 'filename'
                if (numb.gt.1) then
                  do i=1,nOpenFiles
                    if (i.eq.(streamhandle(RESTART)-40)) then
                      tempname=openfiles(i)(1:lineLength(openfiles(i)))//".prev"
                      write(cmd,"(4A)")"mv ",openfiles(i)(1:lineLength(openfiles(i)))," ",tempname(1:lineLength(tempname))
                      call system(cmd)
                      fileopened(i)=.false. 
                    endif
                  enddo         
                endif
                ilgth = lineLength(words(3))
                streamhandle(RESTART) = openFile(words(3)(1:ilgth))
c               restart files will never be appended
                fileappend(streamhandle(RESTART)-40) = .FALSE.
              endif
              if (numb.le.1) call prettyLine(INFO,"Writing restart file")
c              if (numb.ne.0) then
c                write(tempstring,"(A,I10)")"RESTART FILE ",numb
c                call prettyline(RESTART, tempstring)
c              endif  
              call writeRestart()
            else
              call printLine(WARNING,"Cannot interpret restart line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif  

c RESULTS --------------------------------------------------------------------

          else if (words(1).eq.'averages' .or. words(1).eq.'results') then
            newfile = .false.
c           format is;
c           'averages reset|write (filename)'
            if (nwords.lt.2) then
              call printLine(WARNING,"Could not interpret averages line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
              goto 100
            endif

            if (words(2).eq.'reset') then
              if (numb.le.1) call printLine(INFO,"Resetting free energy average to zero.")
              nsteps = 0
            else if ((words(2).eq.'write').or.(words(2).eq.'writeinst')) then
              
              if (nwords.ge.3) then
c               the new file is the third word
                ilgth = lineLength(words(3))
                streamhandle(RESULTS) = openFile(words(3)(1:ilgth))
                fileappend(streamhandle(RESULTS)-40) = (append.or.appendstreams)
                newfile = .TRUE.
              endif
              if (.not.newfile) then
                if (getOption('streamresults')) then
                  streamhandle(RESULTS) = openFile(words(1))
                else
                  streamhandle(RESULTS) = STDOUT
                endif  
                call split(line)
              endif    
              if (numb.gt.1) then
                call printline(RESULTS, "  -  ")
              endif
              if (numb.ne.0) then
                write(tempstring,"(A,I10)")"RESULTS FILE ",numb
                call prettyline(RESULTS, tempstring)
              endif  
              if (words(2)(1:9).eq.'writeinst') then
                call printAllEnergies
              else 
                call printAverages
              endif

            else
              call printLine(WARNING,"Could not interpret averages line")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif 
            if (.not.newfile) then
              if (getOption('streampdb')) then
                streamhandle(PDB) = openFile(words(1))
              else
                streamhandle(PDB) = STDOUT
              endif  
              call split(line)
            endif  

          else if (words(1).eq.'pmeenergy') then
            call pmePerturbeEnergy(nrg,nrgb,nrgf)                   
          endif   

  100 return

      end

