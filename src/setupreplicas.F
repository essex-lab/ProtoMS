      subroutine setupReplicas
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'printer.inc'
      include 'keyset.inc'
      include 'flags.inc'
      include 'parallel.inc'
      include 'mpif.h'

      integer i,j,lineLength,key,k,ind
      integer mpistat(MPI_STATUS_SIZE),ierr
      logical getOption
      double precision getDouble
      character*255 cmd,tmpstr,mylamcom

      doLambdaRE = .false.
      doTemperatureRE = .false.
      doMultiGCMC = .false.
      doMultiJaws2 = .false.

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
          temval(i-1) = getDouble(words(i),0.0)
        enddo
        if (ntasks.lt.ntemperatures) then
          write(printstring,*)"Number of temperatures is lower than the number of tasks!",ntasks,ntemperatures
          call closeProgram(100,printstring)
        endif  
        doTemperatureRE = .true.
c       Storing number of temperature ladders needed and where they should be applied
        ntemladders = 1
        if (getoption('temperatureladder')) then
          ntemladders = nwords - 1
          if (words(2).eq.'lambda') then
            do i=2,nwords
              lamwithlad(i) = getDouble(words(i),-1.0)
              if (lamwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for lambda ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(2).eq.'adams') then
            do i=2,nwords
              gcwithlad(i) = getDouble(words(i),-1.0)
              if (gcwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for adams ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(2).eq.'jbias') then
            do i=2,nwords
              jawswithlad(i) = getDouble(words(i),-1.0)
              if (jawswithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate temperature ladder for jbias ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          else
            write(printstring,*)"Cannot generate temperature ladder for ",words(2),". Incorrect value."
            call closeProgram(100,printstring)
          endif
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
        do i=1,ntemladders
          if (((i-1)*ntemperatures).lt.myid.le.(i*ntemperatures)) then
c           all these ids should belong to the same temperature ladder (same group of temperature replicas)
            write(mytemcom,*)"temladder",i
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,i,0,mytemcom,ierr)
            call MPI_COMM_RANK(group_com,mytemid,ierr)
          else
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,0,MPI_COMM_NULL,ierr)
          endif
        enddo
        if (mytemid.eq.0) then
          do i=1,ntemperatures  
            if (myid.eq.0) then         
              write(cmd,"(3A,F5.3,A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/t-",temval(i),'/'
              call system(cmd) 
            endif
            task2temperature(i) = i
            temperature2task(i) = i-1
            do j=1,ntemladders
              write(thistemcom,*)"temladder",j+1
              if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,thistemcom,ierr)
            enddo
          enddo
c          j = myid/ntemperatures
c          write(mytemcom,*)"temladder",j+1
          mytemperature = 1
          write(trefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(1)
          ntretry = 0
          ntreswap = 0
        elseif (mytemid.gt.0) 
          j = myid/ntemperatures
          if (j.lt.ntemladders) then
            write(mytemcom,*)"temladder",j+1
            call MPI_Recv(mytemperature,1,MPI_INT,0,11,mytemcom,mpistat,ierr)
            write(trefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(mytemperature)
          endif
        endif
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
      elseif (getOption('lambdare')) then      
       Read in frequency and number of lambdas
        lrefreq = nint(getDouble(words(1),1E5))
        nlambdas = nwords - 1
        do i=2,nwords
          lambdaval(i-1) = getDouble(words(i),0.0)
        enddo
        if (ntasks.lt.nlambdas) then
          write(printstring,*)"Number of lambdas is lower than the number of tasks!",nlamtasks,nlambdas
          call closeProgram(100,printstring)
        endif  
        doLambdaRE = .true.
c       Storing number of lambda ladders needed and where they should be applied
        nlamladders = 1
        if (getoption('lambdaladder')) then
          nlamladders = nwords - 1
          if (words(2).eq.'temperature') then
            do i=2,nwords
              temwithlad(i) = getDouble(words(i),-1.0)
              if (temwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for temperature ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(2).eq.'adams') then
            do i=2,nwords
              gcwithlad(i) = getDouble(words(i),-1.0)
              if (gcwithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for adams ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          elseif (words(2).eq.'jbias') then
            do i=2,nwords
              jawswithlad(i) = getDouble(words(i),-1.0)
              if (jawswithlad(i).eq.-1.0) then
                write(printstring,*)"Cannot generate lambda ladder for jbias ",words(i),". Incorrect value."
                call closeProgram(100,printstring)
              endif
            enddo
          else
            write(printstring,*)"Cannot generate lambda ladder for ",words(2),". Incorrect value."
            call closeProgram(100,printstring)
          endif
        endif
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out"
        endif
c       The combination fo temperature and lambda groups goes below
c       Here I need to take care of generating the subgroups
c       I am currently creating the subgroups one at a time in the loop
        mylamid = -1
        do i=1,nlamladders
          if (doTemperaturesRE.eq..true.) then
c           Matching tasks wich should belong to both a lambda ladder and a temperature ladder
            if (mytemid.ge.0).and.(mytemperature.eq.temwithlad(i)) then
              do j=1,ntemladders
                if (mytemcom.eq.("temladder"+j)) then
                  do k=1,nlambdas
                    if (lambdavals(k).eq.lamwithladder(j)) then
                      key = k
                    endif
                  enddo
                endif
              enddo
              write(mylamcom,*)"lamladder",i
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,i,key,mylamcom,ierr)
              call MPI_COMM_RANK(group_com,mylamid,ierr)
              run_lambda = .true.
c           Assigning the rest of the tasks for lambda values
            elseif (ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)).lt.myid.le.(ntemladders*ntemperatures+(nlambdas-ntemladders)*i) then
              key = myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)
              do j=1,ntemladders
                do k=1,nlambdas
                  if (lambdavals(k).eq.lamwithladder(j)) then 
                    ind = k
                  endif
                enddo
                if ((myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1))).gt.ind) then
                  key = myid-(ntemladders*ntemperatures+(nlambdas-ntemladders)*(i-1)+j
                endif
              enddo
              write(mylamcom,*)"lamladder",i
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,i,key,mylamcom,ierr)
              call MPI_COMM_RANK(group_com,mylamid,ierr)
              run_lambda = .true.
            else
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,0,MPI_COMM_NULL,ierr)
            endif
          else
            if (((i-1)*nlambdas).lt.myid.le.(i*nlambdas)) then
c             all these ids should belong to the same lambda ladder
              write(mylamcom,*)"lamladder",i
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,i,0,mylamcom,ierr)
              call MPI_COMM_RANK(mylamcom,mylamid,ierr)
              run_lambda = .true.
            else
              call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,0,MPI_COMM_NULL,ierr)
            endif
          endif 
        enddo
        if (mylamid.eq.0) then
          do i=1,nlambdas  
            if (doTemperatureRE.eq..true.) then
              if (mytemid.eq.0) then 
                write(cmd,"(3A,F5.3,A)")"mkdir -p ",trefolder(1:lineLength(trefolder)),"/lam-",lamval(i),'/')
                call system(cmd) 
              endif
            else 
              if (myid.eq.0) then      
                write(cmd,"(3A,F5.3,A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/lam-",lamval(i),'/')
                call system(cmd) 
              endif
            endif
c           Here is where I need to correlate a lambda to a group. Not done yet
            task2lambda(i) = i
            lambda2task(i) = i-1
            do j=1,nlamladders
              write(thislamcom,*)"lamladder",j+1
              if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,thislamcom,ierr)
            enddo
          enddo
          mylambda = 1
          write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",lamval(1)
          nlretry = 0
          nlreswap = 0
        elseif (mylamid.gt.0) 
          if (run_lambda.eq..true.) then
            call MPI_Recv(mylambda,1,MPI_INT,0,11,mylamcom,mpistat,ierr)
            if (doTemperatureRE.eq..true.) then
              write(lrefolder,"(2A,F5.3)")trefolder(1:lineLength(trefolder)),"/lam-",lamval(mylambda)
            else
              write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",lamval(mylambda)
            endif
          endif
        endif
c      elseif (getOption('lambdare')) then      
c       Read in frequency and number of lambdas
c        lrefreq = nint(getDouble(words(1),1E5))
c        nlambdas = nwords - 1
c        do i=2,nwords
c          lambdaval(i-1) = getDouble(words(i),0.0)
c        enddo
c        if (ntasks.lt.nlambdas) then
c          write(printstring,*)"Number of lambdas is lower than the number of tasks!",nlamtasks,nlambdas
c          call closeProgram(100,printstring)
c        endif  
c        doLambdaRE = .true.
cc       Read in template folder    
c        if (getOption('outfolder')) then
c          outfolder = words(1)
c        else 
c          outfolder = "out"
c        endif
c        if (myid.eq.0) then
c          do i=1,nlambdas           
c            write(cmd,"(3A,F5.3,A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(i),'/'
c            call system(cmd) 
c            task2lambda(i) = i
c            lambda2task(i) = i-1
c            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,'lambdare',ierr)
c          enddo
c          mylambda = 1
c          write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(1)
cc          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
c          nlretry = 0
c          nlreswap = 0
c        else
c          call MPI_Recv(mylambda,1,MPI_INT,0,11,'lambdare',mpistat,ierr)
c          write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(mylambda)
cc          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
c        endif
      elseif (getoption('multigcmc')) then
c       Read in all B-values
        nadams = nwords
        do i=1,nwords
          adamparam(i) = getDouble(words(i),0.0)
        enddo     
        if (ntasks.lt.nadams) then
          call closeProgram(100,"Number of Adam parameters is lower than the number of tasks!")
        endif             
        doMultiGCMC = .true.
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out_gcmc"
        endif
        if (myid.eq.0) then
          do i=1,nadams      
            write(tmpstr,"(3A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/b_"
            call writedoublestring(tmpstr,adamparam(i),cmd)
            call system(cmd) 
            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,'multigcmc',ierr)
          enddo
          myadam = 1
          write(tmpstr,"(2A)")outfolder(1:lineLength(outfolder)),"/b_"
          call writedoublestring(tmpstr,adamparam(1),gcmcfolder)
          write(6,*)"Task = ",myid," folder = ",gcmcfolder(1:lineLength(gcmcfolder))
        else
          call MPI_Recv(myadam,1,MPI_INT,0,11,'multigcmc',mpistat,ierr)
          write(tmpstr,"(2A)")outfolder(1:lineLength(outfolder)),"/b_"
          call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
          write(6,*)"Task = ",myid," folder = ",gcmcfolder(1:lineLength(gcmcfolder))
        endif
      elseif (getoption('multijaws2')) then
c       Read in all jbiases
        njaws2 = nwords
        do i=1,nwords
          alljbias(i) = getDouble(words(i),0.0)
        enddo     
        if (ntasks.lt.njaws2) then
          call closeProgram(100,"Number of bias parameters is lower than the number of tasks!")
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
            write(tmpstr,"(3A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/j_"
            call writedoublestring(tmpstr,alljbias(i),cmd)
            call system(cmd) 
            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,'multijaws2',ierr)
          enddo
          myjbias = 1
          write(tmpstr,"(2A)")outfolder(1:lineLength(outfolder)),"/j_"
          call writedoublestring(tmpstr,alljbias(1),jaws2folder)
          write(6,*)"Task = ",myid," folder = ",jaws2folder(1:lineLength(jaws2folder))
        else
          call MPI_Recv(myjbias,1,MPI_INT,0,11,'multijaws2',mpistat,ierr)
          write(tmpstr,"(2A)")outfolder(1:lineLength(outfolder)),"/j_"
          call writedoublestring(tmpstr,alljbias(myjbias),jaws2folder)
          write(6,*)"Task = ",myid," folder = ",jaws2folder(1:lineLength(jaws2folder))
        endif
      endif

      appendto = .False.

      end
