      subroutine setupReplicas
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'printer.inc'
      include 'keyset.inc'
      include 'flags.inc'
      include 'parallel.inc'
      include 'sep.inc'
      include 'mpif.h'

      integer i,j,lineLength,key,k,ind
      integer mpistat(MPI_STATUS_SIZE),ierr
      logical getOption,closepms
      double precision getDouble
      character*255 cmd,tmpstr,thistfolder

      doLambdaRE = .false.
      doTemperatureRE = .false.
      doMultiGCMC = .false.
      doMultiJaws2 = .false.
      closepms = .false.

      ntemladders = 0
      nlamladders = 0
      ngcladders = 0
      njawsladders = 0

      if (ntasks.lt.2) return

      if (getoption('temperaturere')) then
c       Read in frequency and number of temperatures
        trefreq = nint(getDouble(words(1),1E5))
        ntemperatures = nwords - 1
        do i=2,nwords
          temval(i-1) = getDouble(words(i),25.0)
c          pt_beta_t(1,i-1) = dble(1.0) / (GASR*temval(i-1))
        enddo
        if (ntasks.lt.ntemperatures) then
          write(printstring,*)"Number of tasks is lower than the number of temperatures!",ntasks,ntemperatures
          call closeProgram(100,printstring)
        endif  
        doTemperatureRE = .true.
c       Storing number of temperature ladders needed and where they should be applied
        ntemladders = 1
        if (getoption('temperatureladder')) then
          ntemladders = nwords - 1
          if (words(1).eq.'lambda') then
            do i=2,nwords
              lamwithlad(i-1) = getDouble(words(i),-1.0)
              if (lamwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for lambda",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(1).eq.'gcmcre') then
            do i=2,nwords
              gcwithlad(i-1) = getDouble(words(i),-1.0)
              if (gcwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for adams ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(1).eq.'jbias') then
            do i=2,nwords
              jawswithlad(i-1) = getDouble(words(i),-1.0)
              if (jawswithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for jbias ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          else
            write(printstring,*)"Cannot generate temperature ladder for ",words(2),". Incorrect value."
            call closeProgram(100,printstring)
          endif
        elseif (getoption('lambdare')) then
          nlamladders = nwords - 1
          do i=2,nwords
            lamwithlad(i-1) = getDouble(words(i),-1.0)
            if (lamwithlad(i).eq.-1.0) then
              write(printstring,*)"Cannot generate temperature ladder for lambda ",words(i),". Incorrect value."
              call closeProgram(100,printstring)
            endif     
          enddo
        endif
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out"
        endif
c       Here I need to take care of generating the subgroups
c       I am currently creating the subgroups one at a time in the loop
        mytemid = -1
        mytemladder = -1
        do i=1,ntemladders
          if (((i-1)*ntemperatures).le.myid.and.myid.lt.(i*ntemperatures)) then
c           all these ids should belong to the same temperature ladder (same group of temperature replicas)
            mytemladder = i
          endif
        enddo


        call assigntemperaturere(0,.false.,.true.)


        do i=1,ntemperatures   
          task2temperature(i) = i
          temperature2task(i) = i-1
          ntretry = 0
          ntreswap = 0
        enddo

c        do i=1,nkeys
c          if (keys(i)(1:13).eq.'temperaturere') then
c            if (ntemtasks.lt.2) then
c              write(6,*)"Error: cannot run temperature-RE with this few tasks!"
c              call exit(0)
c            endif  
c            call split(values(i))
c            trefreq = nint(getDouble(words(1),1E5))
c            ntemperatures(1) = nwords - 2
c            do j=3,nwords
c              tempval(1,j-2) = getDouble(words(j),0.0)
c            enddo            
c          endif          
c        enddo
c        if (ntemperatures(1).gt.1) then
c          if (getOption('outfolder')) then
c            outfolder = words(1)
c          else 
c            outfolder = "out"
c          endif
c         Here, I should create folders and assign temperatures          
c        endif
      endif
      if (getOption('lambdare')) then      
c       Read in frequency and number of lambdas
        lrefreq = nint(getDouble(words(1),1E5))
        nlambdas = nwords - 1
        do i=2,nwords
          lambdaval(i-1) = getDouble(words(i),0.0)
        enddo
c       Storing number of lambda ladders needed and where they should be applied
        nlamladders = 1
c 		GR: Simplifying the syntax to generate ladders for GCMC and lambda RE.
        if (getoption('multigcmc')) then
          nlamladders = nwords - 1
          do i=2,nwords
            gcwithlad(i-1) = getDouble(words(i),-1.0)
            if (gcwithlad(i).eq.-1.0) then
              write(printstring,*)"Cannot generate lambda ladder for gcmc ",words(i),". Incorrect value."
              call closeProgram(100,printstring)
            endif
          enddo
        endif     
        if (getoption('lambdaladder')) then
          nlamladders = nwords - 1
          if (words(1).eq.'temperature') then
            do i=2,nwords
              temwithlad(i-1) = getDouble(words(i),-1.0)
              if (temwithlad(i-1).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for temperature ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(1).eq.'gcmcre') then
            do i=2,nwords
              gcwithlad(i-1) = getDouble(words(i),-1.0)
              if (gcwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for adams ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(1).eq.'jbias') then
            do i=2,nwords
              jawswithlad(i-1) = getDouble(words(i),-1.0)
              if (jawswithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for jbias ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          else
            write(printstring,*)"Cannot generate lambda ladder for ",words(2),". Incorrect value."
            call closeProgram(100,printstring)
          endif
        elseif (getoption('multigcmc')) then
          nlamladders = nwords - 1
          do i=2,nwords
            gcwithlad(i-1) = getDouble(words(i),-1.0)
            if (gcwithlad(i).eq.-1.0) then
              write(printstring,*)"Cannot generate lambda ladder for gcmc ",words(i),". Incorrect value."
              call closeProgram(100,printstring)
            endif     
          enddo
        elseif (getoption('temperature')) then
          nlamladders = nwords - 1
          do i=2,nwords
            temwithlad(i-1) = getDouble(words(i),-1.0)
            if (temwithlad(i).eq.-1.0) then
              write(printstring,*)"Cannot generate lambda ladder for temperature ",words(i),". Incorrect value."
              call closeProgram(100,printstring)
            endif     
          enddo
        endif
        if (ntasks.lt.nlambdas) then
          write(printstring,*)"Number of tasks is lower than the number of lambdas!",ntasks,nlambdas
          call closeProgram(100,printstring)
        elseif (doTemperatureRE.and.(ntasks.lt.(nlambdas*nlamladders+ntemperatures*ntemladders-nlamladders*ntemladders))) then 
          write(printstring,*)"Number of tasks is lower than the number required for the lambdare-temperaturere run!"
     .                        ,ntasks,(nlambdas*nlamladders+ntemperatures*ntemladders-nlamladders*ntemladders)
          call closeProgram(100,printstring)
        endif  
        doLambdaRE = .true.
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out"
        endif

c       The combination for temperature and lambda groups goes below
c       Here I need to take care of generating the subgroups
c       I am currently creating the subgroups one at a time in the loop
        mylamid = -1
        mylamladder = -1
        key = -1
        do i=1,nlamladders
          if (doTemperatureRE) then
c           Matching tasks wich should belong to both a lambda ladder and a temperature ladder
            if ((mytemid.ge.ZERO).and.(temval(mytemperature).eq.temwithlad(i))) then
              do j=1,ntemladders
                if (mytemladder.eq.j) then
                  do k=1,nlambdas
                    if (lambdaval(k).eq.lamwithlad(j)) then
                      key = k-1
                    endif
                  enddo
                endif
              enddo
              mylamladder = i
c           Assigning the rest of the tasks for lambda values
            elseif ((ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)).le.myid.and.
     .             myid.lt.(ntemladders*ntemperatures+(nlambdas-ntemladders)*i)) then
              key = myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)) + 1
              do j=1,ntemladders
                do k=1,nlambdas
                  if (lambdaval(k).eq.lamwithlad(j)) then 
                    ind = k-1
                  endif
                enddo
                if ((myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1))).gt.ind) then
                  key = myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)+j)
                endif
              enddo
              mylamladder = i
            endif
          else
            if (((i-1)*nlambdas).le.myid.and.myid.lt.(i*nlambdas)) then
c             all these ids should belong to the same lambda ladder
              mylamladder = i
            endif
          endif
        enddo

        call assignlambdare(key,.false.,.true.)


        do i=1,nlambdas  
c         Here is where I need to correlate a lambda to a group.
          task2lambda(i) = i
          lambda2task(i) = i-1
        enddo
        nlretry = 0
        nlreswap = 0
c     End of assigning ladders
      endif
      if (getoption('multigcmc')) then
        gcrefreq = nint(getDouble(words(1),1E5))
c       Read in all B-values
        nadams = nwords - 1
        do i=2,nwords
          adamparam(i-1) = getDouble(words(i),0.0)
        enddo     
        if (ntasks.ne.nadams) then
          call closeProgram(100,"Number of Adam parameters is not equal to the number of tasks!")
        endif             
        doMultiGCMC = .true.
c 		GR: Simplifying the syntax to generate ladders for GCMC and lambda RE.
        if (getoption('lambdare')) then
          ngcladders = nwords - 1
          do i=2,nwords
            lamwithlad(i-1) = getDouble(words(i),-1.0)
            if (lamwithlad(i).eq.-1.0) then
              write(printstring,*)"Cannot generate gcmc ladder for lambda",words(i),". Incorrect value."
              call closeProgram(100,printstring)
            endif
          enddo
        endif
c       Read in template folder    
c 		GR: Remove folder creating and modify assignadamre.F amd call assignadamre.F
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out_gcmc"
        endif
        if (myid.eq.0) then
          do i=1,nadams
            write(tmpstr,"(4A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
            call writedoublestring(tmpstr,adamparam(i),cmd)
            call system(cmd) 
            if (i.gt.1) call MPI_Send(i,1,MPI_INTEGER,i-1,11,MPI_COMM_WORLD,ierr)
          enddo
          myadam = 1
          write(tmpstr,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(1),gcmcfolder)
        else
          call MPI_Recv(myadam,1,MPI_INTEGER,0,11,MPI_COMM_WORLD,mpistat,ierr)
          write(tmpstr,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
        endif
c Assigning gc ladders with lambda ladders
        do i=1,ngcladders
          if (doLambdaRE) then
c           Matching tasks wich should belong to both a lambda ladder and a temperature ladder
            if ((mylamid.ge.ZERO).and.(lambdaval(mylambda).eq.lamwithlad(i))) then
              do j=1,nlamladders
                if (mylamladder.eq.j) then
                  do k=1,nadams
                    if (adamparam(k).eq.gcwithlad(j)) then
                      key = k-1
                    endif
                  enddo
                endif
              enddo
              mygcladder = i
            endif
          else
            if (((i-1)*nadams).le.myid.and.myid.lt.(i*nadams)) then
c             all these ids should belong to the same lambda ladder
              mygcladder = i
            endif
          endif
        enddo
c 		Finished assigning gcmc ladders.
        do i=1,nadams  
c         Here is where I need to correlate a lambda to a group.
          task2adams(i) = i
          adams2task(i) = i-1
        enddo
        ngcretry = 0
        ngcreswap = 0
      endif

      if (getoption('multijaws2')) then
c       Read in all jbiases
        njaws2 = nwords
        do i=1,nwords
          alljbias(i) = getDouble(words(i),0.0)
        enddo     
        if (ntasks.ne.njaws2) then
          call closeProgram(100,"Number of bias parameters is not equal to the number of tasks!")
        endif             
        doMultiJaws2 = .true.
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out_jaws2"
        endif
        if (myid.eq.0) then
          do i=1,njaws2      
            write(tmpstr,"(4A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),PATHSEP,"j_"
            call writedoublestring(tmpstr,alljbias(i),cmd)
            call system(cmd) 
            if (i.gt.1) call MPI_Send(i,1,MPI_INTEGER,i-1,11,MPI_COMM_WORLD,ierr)
          enddo
          myjbias = 1
          write(tmpstr,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,"j_"
          call writedoublestring(tmpstr,alljbias(1),jaws2folder)
        else
          call MPI_Recv(myjbias,1,MPI_INTEGER,0,11,MPI_COMM_WORLD,mpistat,ierr)
          write(tmpstr,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,"j_"
          call writedoublestring(tmpstr,alljbias(myjbias),jaws2folder)
        endif
      endif


      if (doLambdaRE.and.doTemperatureRE.and.(mytemid.eq.-1).and.(mylamid.eq.-1)) then
        closepms = .true.
      elseif (doLambdaRE.and..not.doTemperatureRE.and.(mylamid.eq.-1)) then
        closepms = .true.
      elseif (.not.doLambdaRE.and.doTemperatureRE.and.(mytemid.eq.-1)) then
        closepms = .true.
      endif

      call MPI_BCAST(closepms,1,MPI_LOGICAL,ntasks-1,MPI_COMM_WORLD,ierr)

      if (closepms) call closeProgram(100,"Number of tasks higher than required!")

      appendto = .False.

      end
