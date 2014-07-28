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

      if (ntasks.lt.2) return

      if (getOption('lambdare')) then      
c       Read in frequency and number of lambdas
        lrefreq = nint(getDouble(words(1),1E5))
        nlambdas = nwords - 1
        do i=2,nwords
          lambdaval(i-1) = getDouble(words(i),0.0)
        enddo
        if (ntasks.ne.nlambdas) then
          call closeProgram(100,"Number of lambdas is not equal to the number of tasks!")
        endif  
        doLambdaRE = .true.
c       Read in template folder    
        if (getOption('refolder')) then
          reoutfolder = words(1)
        else 
          reoutfolder = "out"
        endif
        if (myid.eq.0) then
          do i=1,nlambdas           
            write(cmd,"(3A,F5.3,A)")"mkdir -p ",reoutfolder(1:lineLength(reoutfolder)),"/lam-",lambdaval(i),'/'
            call system(cmd) 
            task2lambda(i) = i
            lambda2task(i) = i-1
            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,MPI_COMM_WORLD,ierr)
          enddo
          mylambda = 1
          write(lrefolder,"(2A,F5.3)")reoutfolder(1:lineLength(reoutfolder)),"/lam-",lambdaval(1)
          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
          nlretry = 0
          nlreswap = 0
        else
          call MPI_Recv(mylambda,1,MPI_INT,0,11,MPI_COMM_WORLD,mpistat,ierr)
          write(lrefolder,"(2A,F5.3)")reoutfolder(1:lineLength(reoutfolder)),"/lam-",lambdaval(mylambda)
          write(6,*)"Task = ",myid," folder = ",lrefolder(1:lineLength(lrefolder))
        endif
      elseif (getoption('multigcmc')) then
c       Read in all B-values
        nadams = nwords
        do i=1,nwords
          adamparam(i) = getDouble(words(i),0.0)
        enddo     
        if (ntasks.ne.nadams) then
          call closeProgram(100,"Number of Adam parameters is not equal to the number of tasks!")
        endif             
        doMultiGCMC = .true.
c       Read in template folder    
        if (getOption('refolder')) then
          reoutfolder = words(1)
        else 
          reoutfolder = "out_gcmc"
        endif
        if (myid.eq.0) then
          do i=1,nadams      
            write(tmpstr,"(3A)")"mkdir -p ",reoutfolder(1:lineLength(reoutfolder)),"/b_"
            call writedoublestring(tmpstr,adamparam(i),cmd)
            call system(cmd) 
            if (i.gt.1) call MPI_Send(i,1,MPI_INT,i-1,11,MPI_COMM_WORLD,ierr)
          enddo
          myadam = 1
          write(tmpstr,"(2A)")reoutfolder(1:lineLength(reoutfolder)),"/b_"
          call writedoublestring(tmpstr,adamparam(1),gcmcfolder)
          write(6,*)"Task = ",myid," folder = ",gcmcfolder(1:lineLength(gcmcfolder))
        else
          call MPI_Recv(myadam,1,MPI_INT,0,11,MPI_COMM_WORLD,mpistat,ierr)
          write(tmpstr,"(2A)")reoutfolder(1:lineLength(reoutfolder)),"/b_"
          call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
          write(6,*)"Task = ",myid," folder = ",gcmcfolder(1:lineLength(gcmcfolder))
        endif
      else
        do i=1,nkeys
          if (keys(i)(1:13).eq.'temperaturere') then
            if (ntasks.lt.2) then
              write(6,*)"Error: cannot run temperature-RE with this few tasks!"
              call exit(0)
            endif  
            call split(values(i))
            trefreq = nint(getDouble(words(1),1E5))
            ntemperatures(1) = nwords - 2
            do j=3,nwords
              tempval(1,j-2) = getDouble(words(j),0.0)
            enddo            
          endif          
        enddo
        if (ntemperatures(1).gt.1) then
          if (getOption('refolder')) then
            reoutfolder = words(1)
          else 
            reoutfolder = "out"
          endif
c         Here, I should create folders and assign temperatures          
        endif
      endif

      appendto = .False.

      end
