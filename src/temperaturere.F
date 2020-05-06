      subroutine temperatureRE
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'gb.inc'
      include 'averages.inc'
      include 'printer.inc'
      include 'parallel.inc'
      include 'splitter.inc'
      include 'sep.inc'
      include 'mpif.h' 
      include 'metropolis.inc'
c###################################################################
c
c     This routine peforms temperature-replica exchange
c
c     Ana I. Cabedo Martinez, 2015
c     from lambdare.F (S. Genheden 2014)
c
c###################################################################

      integer i,j,ranInt,task1,task2,choice,aidx,bidx,linelength,k
      integer mpistat(MPI_STATUS_SIZE),ierr,nrecieved
      integer lamsbuffer(4),lamrbuffer(4),tsbuffer(nlambdas),lmbuffer(nlambdas)
      logical updated(ntemperatures),reTest
      double precision nrg,nrgb,nrgf,inrgf,inrgb,teml,temh,temp,oldtem
      double precision totalEnergy,extraEnergy,testbeta,resttestEnergy
      double precision allteml(ntemperatures),alltemh(ntemperatures),allnrg(ntemperatures),allbeta(ntemperatures)
      double precision dell,delh
      double precision buffer(4)
      character*255 tempstr,tempstr2

c     check that temperature swaps are only attempted between tasks that are running temperature values
      if (mytemid.lt.0) return

      do i = 1,ntemperatures
        updated(i) = .false.
      enddo

      call closeFiles

      if (mytemperature.eq.1) then
        teml = temval(mytemperature)
      else
        teml = temval(mytemperature-1)
      endif

      if (mytemperature.eq.ntemperatures) then
        temh = temval(mytemperature)
      else
        temh = temval(mytemperature+1)
      endif

      call setAllMoved
      if (.not.doREST) then
        nrg = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
        testbeta = mt_beta_t
      else
        nrg = restTestEnergy()
        testbeta = rest_beta      
      endif
      
      if (mytemid.eq.0) then
        allteml(mytemperature) = teml
        alltemh(mytemperature) = temh
        allnrg(mytemperature) = nrg
        allbeta(mytemperature) = testbeta
        updated(1) = .false.
        nrecieved = 1
        do i=2,ntemperatures
          call MPI_Recv(buffer,4,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11,mytemcom,mpistat,ierr)
          task1 = mpistat(MPI_SOURCE)
          allteml(task2temperature(task1+1)) = buffer(1)
          alltemh(task2temperature(task1+1)) = buffer(2)
          allnrg(task2temperature(task1+1)) = buffer(3)
          allbeta(task2temperature(task1+1)) = buffer(4)
          updated(task1+1) = .false.
        enddo
      else
        buffer(1) = teml
        buffer(2) = temh
        buffer(3) = nrg
        buffer(4) = testbeta
        call MPI_SSend(buffer,4,MPI_DOUBLE_PRECISION,0,11,mytemcom,ierr)
      endif


      if (mytemid.eq.0) then
c       Randomly selects odd/even pairs
        choice = ranInt(0,1)
        do i=1,ntemperatures-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1

          ntretry = ntretry + 1
c          write(6,*)myid," pt_beta_t ",pt_beta_t
          if (reTest(allbeta(aidx),allbeta(bidx),allnrg(aidx),allnrg(bidx))) then

            task1 = temperature2task(aidx)
            task2 = temperature2task(bidx)

            task2temperature(task1+1) = bidx
            task2temperature(task2+1) = aidx
            temperature2task(aidx) = task2
            temperature2task(bidx) = task1

            ntreswap = ntreswap + 1
c           The main task already knows bidx and aidx, so no need to send that information in that case
            if (task1.ne.0) then 
              call MPI_SSend(bidx,1,MPI_INTEGER,task1,11,mytemcom,ierr)
            else
              mytemperature = bidx
            endif
            if (task2.ne.0) then
              call MPI_SSend(aidx,1,MPI_INTEGER,task2,11,mytemcom,ierr)
            else
              mytemperature = aidx
            endif
            updated(task1+1) = .true.
            updated(task2+1) = .true.
          else
c            write(6,*)"Failing ",temval(aidx),temval(bidx),allnrg(aidx),allnrgh(aidx),allnrg(bidx),allnrgl(bidx)
          endif
        enddo

        do i=2,ntemperatures
          if (.not.updated(i)) call MPI_SSend(task2temperature(i),1,MPI_INTEGER,i-1,11,mytemcom,ierr)
        enddo

      else    
         call MPI_Recv(mytemperature,1,MPI_INTEGER,0,11,mytemcom,mpistat,ierr)
      endif

      call MPI_BCAST (ntretry, 1, MPI_INTEGER, 0, mytemcom, ierr)
      call MPI_BCAST (ntreswap, 1, MPI_INTEGER, 0, mytemcom, ierr)
      
      call MPI_BCAST (updated,ntemperatures,MPI_LOGICAL,0,mytemcom,ierr)
      call MPI_BCAST (temperature2task,ntemperatures,MPI_INTEGER,0,mytemcom,ierr)
      call MPI_BCAST (task2temperature,ntemperatures,MPI_INTEGER,0,mytemcom,ierr)
      
      if (.not.doREST) then
        oldtem = temperatureC
        temperatureC = temval(mytemperature)
        temperature = temperatureC + DEG2KEL
        mt_beta_t = dble(1.0) / (GASR*temperature) 
      else 
        oldtem = resttempC
        resttempC = temval(mytemperature)
        resttemp = resttempC + DEG2KEL
        rest_beta = dble(1.0) / (GASR*resttemp)
        restfac1 = rest_beta / mt_beta_t
        restfac2 = (mt_beta_t+rest_beta)/(2.0d0*mt_beta_t)
        do i=1,7
          if (restlevel(i).eq.1) then
            restscaling(i) = restfac1
          else if (restlevel(i).eq.2) then
            restscaling(i) = restfac2
          else
            restscaling(i) = 1.0d0
          endif
        enddo
        call setAllMoved
        call initialiseExtraEnergies(.true.)      
        nrg = totalEnergy(OLD,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)    
        call clearMoveList
        call updateExtraEnergies
      endif

      if (updated(mytemid+1).and.doLambdaRE) then
        lamsbuffer(1) = mylamladder
        lamsbuffer(2) = mylamid
        lamsbuffer(3) = ntretry
        lamsbuffer(4) = ntreswap
        do i=1,ntemperatures
          if (oldtem.eq.temval(i)) then 
            call MPI_Sendrecv(lamsbuffer,4,MPI_INTEGER,temperature2task(i),11,lamrbuffer,4,
     .           MPI_INTEGER,temperature2task(i),11,mytemcom,mpistat,ierr)
            call MPI_Sendrecv(task2lambda,nlambdas,MPI_INTEGER,temperature2task(i),11,tsbuffer,nlambdas,
     .           MPI_INTEGER,temperature2task(i),11,mytemcom,mpistat,ierr)
            call MPI_Sendrecv(lambda2task,nlambdas,MPI_INTEGER,temperature2task(i),11,lmbuffer,nlambdas,
     .           MPI_INTEGER,temperature2task(i),11,mytemcom,mpistat,ierr)
          endif
        enddo
        task2lambda = tsbuffer
        lambda2task = lmbuffer
        mylamladder = lamrbuffer(1)
        mylamid = lamrbuffer(2)
        ntretry = lamrbuffer(3)
        ntreswap = lamrbuffer(4)
      endif

c     Re-setting all running averages    
      nsteps = 0
      
      write(tempstr,"(1A,F6.3)")"t-",oldtem
      write(trefolder,"(3A,F6.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"t-",temval(mytemperature)
      if (oldtem.ge.1000.0) then
        write(tempstr,"(1A,F8.3)")"t-",oldtem
      else if (oldtem.ge.100.0) then
        write (tempstr,"(1A,F7.3)")"t-",oldtem
      endif
      if (temval(mytemperature).ge.1000.0) then
        write(trefolder,"(3A,F8.3)")outfolder(1:lineLength(outfolder)),
     .                                    PATHSEP,"t-",temval(mytemperature)
      else if (temval(mytemperature).ge.100.0) then
        write(trefolder,"(3A,F7.3)")outfolder(1:lineLength(outfolder)),
     .                                    PATHSEP,"t-",temval(mytemperature)
      endif
      if (doLambdaRE) write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lambdaval(mylambda)


c     Replacing the old temperature for the new one in the file name.
      k = 1
      if ((dolambdaRE).or.(doMultiGCMC).or.(doMultiJaws2)) k = 2
      do i=1,nOpenFiles    
        if (wasopened(i)) fileappend(i) = .true.
        if (i.eq.streamhandle(RESTART)-40) then
          fileappend(i) = .False.
        endif
        fileopened(i) = .False.
        call split2(openfiles(i),PATHSEP)
        tempstr2 = words2(nwords2-k)
        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
          do j=1,nwords2-(k+1)
            words(j) = words2(j)
          enddo
           write(words(nwords2-k),"(1A,F6.3)")"t-",temval(mytemperature)
           if (temval(mytemperature).ge.1000.0) then
             write(words(nwords2-k),"(1A,F8.3)")"t-",temval(mytemperature)
           else if (temval(mytemperature).ge.100.0) then 
             write(words(nwords2-k),"(1A,F7.3)")"t-",temval(mytemperature)
           endif
           do j=1,k
             words(nwords2-(j-1)) = words2(nwords2-(j-1))
           enddo
           nwords = nwords2
           call join2(1,nwords,PATHSEP,openfiles(i))
        endif
      enddo

      end subroutine
