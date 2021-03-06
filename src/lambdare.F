      subroutine lambdaRE
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
c###################################################################
c
c     This routine peforms lambda-replica exchange
c
c     Samuel Genheden, 2014
c     modified - A. Cabedo Martinez, 2015
c
c###################################################################

      integer i,j,ranInt,task1,task2,choice,aidx,bidx,linelength,k
      integer mpistat(MPI_STATUS_SIZE),ierr,nrecieved
      integer temsbuffer(4),temrbuffer(4),tsbuffer(ntemperatures),tmbuffer(ntemperatures)
      integer gcsbuffer(4),gcrbuffer(4),asbuffer(nadams),ambuffer(nadams)
      logical updated(nlambdas),retiTest
      double precision nrg,nrgb,nrgf,lamb,lamf,temp,inrgf,inrgb,oldlam
      double precision totalEnergy,extraEnergy
      double precision allnrg(nlambdas),allnrgb(nlambdas),allnrgf(nlambdas)
      double precision alllamb(nlambdas),alllamf(nlambdas)
      double precision delf,delb
      double precision buffer(5)
      character*255 tempstr,tempstr2

c     check that lambda swaps are only attempted between tasks that are running lambda values
      if (mylamid.lt.0) return

c     The loop below should not be needed
      do i = 1,nlambdas
        updated(i) = .false.
      enddo
   
      call closeFiles

      if (mylambda.eq.1) then
        lamb = lambdaval(mylambda)
      else
        lamb = lambdaval(mylambda-1)
      endif

      if (mylambda.eq.nlambdas) then
        lamf = lambdaval(mylambda)
      else
        lamf = lambdaval(mylambda+1)
      endif

      call LREEnergy(lamb,lamf,nrg,nrgb,nrgf)


      if (mylamid.eq.0) then
        allnrg(mylambda) = nrg
        allnrgb(mylambda) = nrgb
        allnrgf(mylambda) = nrgf
        alllamf(mylambda) = lamf
        alllamb(mylambda) = lamb
        updated(1) = .false.
        nrecieved = 1
        do i=2,nlambdas
          call MPI_Recv(buffer,5,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11,mylamcom,mpistat,ierr)
          task1 = mpistat(MPI_SOURCE)
          allnrg(task2lambda(task1+1)) = buffer(1)
          allnrgb(task2lambda(task1+1)) = buffer(2)
          allnrgf(task2lambda(task1+1)) = buffer(3)
          alllamb(task2lambda(task1+1)) = buffer(4)
          alllamf(task2lambda(task1+1)) = buffer(5)
          updated(task1+1) = .false.
        enddo
      else
        buffer(1) = nrg
        buffer(2) = nrgb
        buffer(3) = nrgf
        buffer(4) = lamb
        buffer(5) = lamf
        call MPI_SSend(buffer,5,MPI_DOUBLE_PRECISION,0,11,mylamcom,ierr)
      endif


      if (mylamid.eq.0) then
c       Randomly selects odd/even pairs
        choice = ranInt(0,1)
        do i=1,nlambdas-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1

          nlretry = nlretry + 1
          if (retiTest(lambdaval(aidx),alllamf(aidx),alllamb(aidx),allnrg(aidx),allnrgf(aidx),allnrgb(aidx),

     .                 lambdaval(bidx),alllamf(bidx),alllamb(bidx),allnrg(bidx),allnrgf(bidx),allnrgb(bidx))) then
c            write(6,*)"Swapping ",lambdaval(aidx),lambdaval(bidx),allnrgb(bidx)-allnrg(aidx),allnrgf(aidx)-allnrg(bidx)
c            write(6,*)allnrg(aidx),allnrgf(aidx),allnrg(bidx),allnrgb(bidx)
c            write(6,*)-(allnrg(bidx)-allnrgf(aidx)-allnrgb(bidx)+allnrg(aidx))
            task1 = lambda2task(aidx)
            task2 = lambda2task(bidx)

            task2lambda(task1+1) = bidx
            task2lambda(task2+1) = aidx
            lambda2task(aidx) = task2
            lambda2task(bidx) = task1
            nlreswap = nlreswap + 1
            if (task1.ne.0) then
              call MPI_SSend(bidx,1,MPI_INTEGER,task1,11,mylamcom,ierr)
            else
              mylambda = bidx
            endif
            if (task2.ne.0) then
              call MPI_SSend(aidx,1,MPI_INTEGER,task2,11,mylamcom,ierr)
            else
              mylambda = aidx
            endif
            updated(task1+1) = .true.
            updated(task2+1) = .true.
          endif
        enddo

        do i=2,nlambdas
          if (.not.updated(i)) call MPI_SSend(task2lambda(i),1,MPI_INTEGER,i-1,11,mylamcom,ierr)
        enddo

      else    
         call MPI_Recv(mylambda,1,MPI_INTEGER,0,11,mylamcom,mpistat,ierr)
      endif

c     The broadcast of nlretry needs to be dealt with in here,
c     since it is ignoring the fact that I have different groups per communicator

      call MPI_BCAST (nlretry, 1, MPI_INTEGER, 0, mylamcom, ierr)
      call MPI_BCAST (nlreswap, 1, MPI_INTEGER, 0, mylamcom, ierr)
      call MPI_BCAST (updated,nlambdas,MPI_LOGICAL,0,mylamcom,ierr)
      call MPI_BCAST (lambda2task,nlambdas,MPI_INTEGER,0,mylamcom,ierr)
      call MPI_BCAST (task2lambda,nlambdas,MPI_INTEGER,0,mylamcom,ierr)


      oldlam = lambda
      lambda = lambdaval(mylambda)

      if (DifferentLambdas) then
        lambdaf = lambda+deltalam
        lambdab = lambda-deltalam
        if (lambdaf.gt.1.0) lambdaf = 1.0
        if (lambdaf.lt.0.0) lambdaf = 0.0
        if (lambdab.gt.1.0) lambdab = 1.0
        if (lambdab.lt.0.0) lambdab = 0.0
      endif


c     The communication in this code is not optimized. Task2temperature, temperature2task,
c     nlretry & nlreswap probably only need to be swaped if a process with temid 0 is involved
c     in the swap, since the reason to swap these is problems associated with the next temperatureswap
      if (updated(mylamid+1).and.doTemperatureRE) then
        temsbuffer(1) = mytemladder
        temsbuffer(2) = mytemid
        temsbuffer(3) = nlretry
        temsbuffer(4) = nlreswap
        do i=1,nlambdas
          if (oldlam.eq.lambdaval(i)) then
            call MPI_Sendrecv(temsbuffer,4,MPI_INTEGER,lambda2task(i),11,temrbuffer,4,
     .           MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
            call MPI_Sendrecv(task2temperature,ntemperatures,MPI_INTEGER,lambda2task(i),11,
     .           tsbuffer,ntemperatures,MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
            call MPI_Sendrecv(temperature2task,ntemperatures,MPI_INTEGER,lambda2task(i),11,
     .           tmbuffer,ntemperatures,MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
          endif
        enddo
        task2temperature = tsbuffer
        temperature2task = tmbuffer
        mytemladder = temrbuffer(1)
        mytemid = temrbuffer(2)
        nlretry = temrbuffer(3)
        nlreswap = temrbuffer(4)
      elseif (updated(mylamid+1).and.doMultiGCMC) then
        gcsbuffer(1) = mygcladder
        gcsbuffer(2) = mygcid
        gcsbuffer(3) = nlretry
        gcsbuffer(4) = nlreswap
        do i=1,nlambdas
          if (oldlam.eq.lambdaval(i)) then
            call MPI_Sendrecv(gcsbuffer,4,MPI_INTEGER,lambda2task(i),11,gcrbuffer,4,
     .           MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
            call MPI_Sendrecv(task2adams,nadams,MPI_INTEGER,lambda2task(i),11,
     .           asbuffer,nadams,MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
            call MPI_Sendrecv(adams2task,nadams,MPI_INTEGER,lambda2task(i),11,
     .           ambuffer,nadams,MPI_INTEGER,lambda2task(i),11,mylamcom,mpistat,ierr)
          endif
        enddo
        task2adams= asbuffer
        adams2task = ambuffer
        mygcladder = gcrbuffer(1)
        mygcid = gcrbuffer(2)
        nlretry = gcrbuffer(3)
        nlreswap = gcrbuffer(4)
      endif

      call recalculateParameters
      call recalculateGeometries    
      call setAllMoved
      call initialiseExtraEnergies(.true.) 
      if (doGB) then 
        call ZeroAllPDACont
        call ZeroPDACont
        call gbsavesystem
      endif     
      nrg = totalEnergy(OLD,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)      
      call clearMoveList
      call updateExtraEnergies

c     Re-setting all running averages    
      nsteps = 0
c		NEW GR: every where there is a k used to be a 1.
      k = 1
      if (doMultiGCMC)then
        k = 2
      endif
      write(tempstr,"(1A,F5.3)")"lam-",oldlam
      if (doTemperatureRE)then
        write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lambdaval(mylambda)
      else
        write(lrefolder,"(3A,F5.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"lam-",lambdaval(mylambda)
      endif
      do i=1,nOpenFiles    
        if (wasopened(i)) fileappend(i) = .true.
        if (i.eq.streamhandle(RESTART)-40) then
          fileappend(i) = .False.
        endif
        fileopened(i) = .False.
        call split2(openfiles(i),PATHSEP)
        tempstr2 = words2(nwords2-k) 
c        write(6,*)myid, "GR lambdare.F tempstr2, tempstr ", tempstr2,tempstr   
        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
c          write(6,*)myid, "GR lambdare.F MATCH! tempstr2 = tempstr ", tempstr2,tempstr 
          do j=1,nwords2-(k+1)
            words(j) = words2(j)
          enddo
          write(words(nwords2-k),"(1A,F5.3)")"lam-",lambda
          do j=1,k
            words(nwords2-(j-1)) = words2(nwords2-(j-1))
          enddo
          nwords = nwords2
          call join2(1,nwords,PATHSEP,openfiles(i))
        endif
      enddo

c      write(6,*)myid," ----In LambdaRE-----"
c      write(6,*)myid, " mytemid ",mytemid," mylamid ",mylamid," mytemladder ",mytemladder," mylamladder ",mylamladder

c      call MPI_FINALIZE(ierr)
c      call exit(0)

      end subroutine
