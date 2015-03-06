      subroutine setupReplicas
      include 'dimensions.inc'
      include 'splitter.inc'
      include 'printer.inc'
      include 'keyset.inc'
      include 'flags.inc'
      include 'parallel.inc'
      include 'mpif.h'

      integer i,j,lineLength
      integer mpistat(MPI_STATUS_SIZE),ierr
      logical getOption
      double precision getDouble
      character*255 cmd,tmpstr

      doLambdaRE = .false.
      doTemperatureRE = .false.
      doMultiGCMC = .false.
      doMultiJaws2 = .false.

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
            write(group_com,*)"temladder",i
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,i,0,group_com,ierr)
            call MPI_COMM_RANK(group_com,mytemid,ierr)
          else
            call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,0,MPI_COMM_NULL,ierr)
          endif
        if (mytemid.eq.0) then
          do i=1,ntemperatures  
            if (myid.eq.0) then         
              write(cmd,"(3A,F5.3,A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/t-",temval(i),'/'
              call system(cmd) 
            endif
            task2temperature(i) = i
            temperature2task(i) = i-1
            j = myid/ntemperatures
            if (j.lt.ntemladders) then
              write(mytemcom,*)"temladder",j
              if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,mytemcom,ierr)
          enddo
          mytemperature = 1
          write(trefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(1)
c          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
          ntretry = 0
          ntreswap = 0
        elseif (mytemid.gt.0) 
          mysubcom = GetSubComunicator()
          call MPI_Recv(mytemperature,1,MPI_INT,0,11,mysubcom,mpistat,ierr)
          write(trefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(mytemperature)
c          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
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
c       Read in frequency and number of lambdas
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
c       Read in template folder    
        if (getOption('outfolder')) then
          outfolder = words(1)
        else 
          outfolder = "out"
        endif
        if (myid.eq.0) then
          do i=1,nlambdas           
            write(cmd,"(3A,F5.3,A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(i),'/'
            call system(cmd) 
            task2lambda(i) = i
            lambda2task(i) = i-1
            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,'lambdare',ierr)
          enddo
          mylambda = 1
          write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(1)
c          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
          nlretry = 0
          nlreswap = 0
        else
          call MPI_Recv(mylambda,1,MPI_INT,0,11,'lambdare',mpistat,ierr)
          write(lrefolder,"(2A,F5.3)")outfolder(1:lineLength(outfolder)),"/lam-",lambdaval(mylambda)
c          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
        endif
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
