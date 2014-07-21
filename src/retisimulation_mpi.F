      subroutine retiSimulation_mpi(id,ntasks)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'mpif.h'
c#######################################################
c
c     Driver for doing replica-exchange thermodynamic
c     integration
c 
c     Samuel Genheden February 2013
c
c#######################################################
      integer id,ntasks,iter,iretiblock,iptblock,totiter,doequil,doreti,dopt,i,j,stat,ierr
      integer mpistat(MPI_STATUS_SIZE),lineLength
      logical fileExists,doquit
      character*300 cmd,tempstr,tempstr2

c     MASTER NODE WORK
      if (id.eq.0) then
        doquit = .false.

c       Return if no iterations where requested, or less than 2 lambda values where given
        if ((niter*nretiblock.lt.1).or.(nlambda.eq.0)) then
          doquit = .true.      
        endif

c       Also return if slave script do not exists
        if (.not.fileExists(slave)) then
          write(printstring,"(A,F5.1,A)")"Could not find slave script! Exit."
          call prettyLine(FATAL,printstring)   
          doquit = .true.            
        endif

c       Check if tasks is at least equal to the number of lambda values
        if (nptblock.eq.0) then
          nminslaves = nlambda
        else
          nminslaves = 0
          do i=1,nlambda
            nminslaves = nminslaves + ntemp(i)
          enddo
        endif
        if (ntasks.lt.nminslaves) then
          write(printstring,"(A,I5,A)") 
     .      "Tasks less than the number of lambda values and temperature replicas (",nminslaves,")! Exit."
          call prettyLine(FATAL,printstring) 
          doquit = .true.
        endif
 
        if (doquit.eqv..true.) then
          call killSlaves(ntasks)
          return
        endif      

c       Create the names of the restart, reti and energy files
        if (nptblock.eq.0) then
          do i=1,nlambda
            write(cmd,"(3A,F4.2,A)")"mkdir -p ",outdir(1:lineLength(outdir)),'/lam-',lambdaval(i),'/'
            call system(cmd)
            write(rstfiles(i,1),"(2A,F4.2,A)")outdir(1:lineLength(outdir)),'/lam-',lambdaval(i),'/restart'
            write(retifiles(i),"(2A,F4.2,A)")outdir(1:lineLength(outdir)),'/lam-',lambdaval(i),'/reti'
          enddo
        else 
          do i=1,nlambda           
            write(tempstr,"(2A)")outdir(1:lineLength(outdir)),'/t-'
            call writeNumberString(tempstr,int(retitemp-DEG2KEL),tempstr2)
            write(tempstr,"(2A,F4.2,A)")tempstr2(1:lineLength(tempstr2)),'/lam-',lambdaval(i),'/'
            write(cmd,"(3A)")"mkdir -p ",tempstr(1:lineLength(tempstr)),'/'
            call system(cmd) 

            write(rstfiles(i,1),"(2A)")tempstr(1:lineLength(tempstr)),'/restart'
            write(retifiles(i),"(2A)")tempstr(1:lineLength(tempstr)),'/reti'
            if (ntemp(i).gt.1) then
              write(pttrajfiles(i),"(2A)")tempstr(1:lineLength(tempstr)),'/temptraj'     
              write(enedistfiles(i),"(2A)")tempstr(1:lineLength(tempstr)),'/enedist'               
              pttrajstream(i) = 50+i
              enediststream(i) = 100+i
              open(unit=pttrajstream(i),file=pttrajfiles(i),status="unknown",form='formatted')              
              open(unit=enediststream(i),file=enedistfiles(i),status="unknown",form='formatted') 
              do j=1,ntemp(i)   
                write(tempstr,"(2A)")outdir(1:lineLength(outdir)),'/t-'               
                call writeNumberString(tempstr,int(pttemperature(i,j)-DEG2KEL),tempstr2)   
                write(tempstr,"(2A,F4.2,A)")tempstr2(1:lineLength(tempstr2)),'/lam-',lambdaval(i),'/'
                if (j.gt.1) then
                  write(cmd,"(3A)")"mkdir -p ",tempstr(1:lineLength(tempstr)),'/'
                  call system(cmd) 
                endif
                write(energyfiles(i,j),"(2A)")tempstr(1:lineLength(tempstr)),'/spenergy'     
                write(rstfiles(i,j),"(2A)")tempstr(1:lineLength(tempstr)),'/restart'                
              enddo
            endif
          enddo
        endif

c       Write out current lambda values and initialize lambdacurr
        do i=1,nlambda
          lambdacurr(i) = lambdaval(i)          
          if ((nptblock.gt.0).and.(ntemp(i).gt.1)) then
            do j=1,ntemp(i)
              pttempcurr(i,j) = pttemperature(i,j)
            enddo
            write(pttrajstream(i),"(I6,50I7)")0,(int(pttemperature(i,j)),j=1,ntemp(i))
          endif
        enddo
        write(printstring,"(I6,50F5.2)")0,(lambdacurr(i),i=1,nlambda)
        call prettyline(INFO,printstring)

c       Loop over all requested iterations
        doequil = 1
        totiter = 0
        nretitry = 0
        nretiswap = 0
        if (nptblock.gt.0) then
          do i=1,nlambda
            npttry(i) = 0
            nptswap(i) = 0
          enddo
        endif
        do iter=1,niter
          do iretiblock=1,nretiblock
            doreti = 0
c           Check if we are doing temperature swapping            
            if (nptblock.gt.0) then
              do iptblock=1,nptblock
                totiter = totiter + 1
                dopt = 0
                if (iretiblock.eq.nretiblock) doreti = 1
                if (iptblock.eq.nptblock) dopt = 1
                call spawnProcesses(totiter,doequil,doreti,dopt)
                call checkOutput(stat)
                if (stat.ne.0) then 
                  call killSlaves(ntasks)
                  return
                endif
                doequil = 0  
              enddo
c             Attempt temperature swapping
              call tempSwap(totiter,stat)
              if (stat.ne.0) then 
                call killSlaves(ntasks)
                return
              endif
              do i=1,nlambda
                if (ntemp(i).gt.1) then
                  write(pttrajstream(i),"(I6,50I7)")totiter,(int(pttempcurr(i,j)),j=1,ntemp(i))
                endif
              enddo
c           If we are not doing temperature swapping           
            else
              totiter = totiter + 1 
              if (iretiblock.eq.nretiblock) doreti = 1
              call spawnProcesses(totiter,doequil,doreti,0)
              call checkOutput(stat)
              if (stat.ne.0) then 
                call killSlaves(ntasks)
                return
              endif
              doequil = 0  
            endif  
          enddo   
c         Attempt lambda swaping          
          call retiSwap(stat)    
          if (stat.ne.0) then 
            call killSlaves(ntasks)
            return
          endif
          write(printstring,"(I6,50F5.2)")totiter,(lambdacurr(i),i=1,nlambda)
          call prettyline(INFO,printstring)
        enddo
    
        write(printstring,"(A,2I6,A,F8.3,A)")
     .  "Attempted/succesful lambda swaps: ",nretitry,nretiswap," (",dble(nretiswap)/dble(nretitry)*100,"%)" 
        call prettyline(INFO,printstring)
        if (nptblock.gt.0) then
          do i=1,nlambda
            if (ntemp(i).gt.1) then
              write(pttrajstream(i),"(A,2I6,A,F8.3,A)")
     .  "Attempted/succesful temperature swaps: ",npttry(i),nptswap(i)," (",dble(nptswap(i))/dble(npttry(i))*100,"%)" 
              close(pttrajstream(i))
              close(enediststream(i))
            endif
          enddo
        endif

        call killSlaves(ntasks)
c     SLAVE WORK
      else
        do
          call MPI_Recv(cmd,300,MPI_CHARACTER,0,11,MPI_COMM_WORLD,mpistat,ierr)
c          write(6,*)"Slave Recv ",cmd(1:lineLength(cmd))
          if (cmd(1:4).eq."Exit") then
            return
          else
            call system(cmd(1:lineLength(cmd)))              
            write(cmd,"(A,I2)")"Done ",id
c            write(6,*)cmd(1:lineLength(cmd))
            call MPI_Send(cmd,300,MPI_CHARACTER,0,11,MPI_COMM_WORLD,ierr)            
          endif
        enddo       
      endif

      end subroutine
