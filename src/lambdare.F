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
      include 'mpif.h'
c###################################################################
c
c     This routine peforms lambda-replica exchange
c
c     Samuel Genheden, 2014
c
c###################################################################

      integer i,j,ranInt,task1,task2,choice,aidx,bidx,linelength
      integer mpistat(MPI_STATUS_SIZE),ierr,nrecieved
      logical updated(nlambdas),retiTest
      double precision nrg,nrgb,nrgf,lamb,lamf,temp,inrgf,inrgb,oldlam
      double precision totalEnergy,extraEnergy
      double precision allnrg(nlambdas),allnrgb(nlambdas),allnrgf(nlambdas)
      double precision alllamb(nlambdas),alllamf(nlambdas)
      double precision delf,delb
      double precision buffer(5)
      character*255 tempstr,tempstr2

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

      if (myid.eq.0) then
        allnrg(mylambda) = nrg
        allnrgb(mylambda) = nrgb
        allnrgf(mylambda) = nrgf
        alllamf(mylambda) = lamf
        alllamb(mylambda) = lamb
        updated(1) = .false.
        nrecieved = 1
        do i=2,nlambdas
          call MPI_Recv(buffer,5,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11,MPI_COMM_WORLD,mpistat,ierr)
          task1 = mpistat(MPI_SOURCE)
          allnrg(task2lambda(task1+1)) = buffer(1)
          allnrgb(task2lambda(task1+1)) = buffer(2)
          allnrgf(task2lambda(task1+1)) = buffer(3)
          alllamb(task2lambda(task1+1)) = buffer(4)
          alllamf(task2lambda(task1+1)) = buffer(5)
          updated(task1+1) = .false.
        enddo
c        do i=1,nlambdas
c          write(6,*)allnrg(i),allnrgb(i),allnrgf(i),alllamb(i),alllamf(i)
c        enddo
      else
        buffer(1) = nrg
        buffer(2) = nrgb
        buffer(3) = nrgf
        buffer(4) = lamb
        buffer(5) = lamf
        call MPI_Send(buffer,5,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,ierr)
      endif

      if (myid.eq.0) then

c       Randomly selects odd/even pairs
        choice = ranInt(0,1)
        do i=1,nlambdas-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1

          nlretry = nlretry + 1
c          write(6,*)nlretry
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
              call MPI_Send(bidx,1,MPI_INT,task1,11,MPI_COMM_WORLD,ierr)
            else
              mylambda = bidx
            endif
            if (task2.ne.0) then
              call MPI_Send(aidx,1,MPI_INT,task2,11,MPI_COMM_WORLD,ierr)
            else
              mylambda = aidx
            endif
            updated(task1+1) = .true.
            updated(task2+1) = .true.
          else
c            write(6,*)"Failing ",lambdaval(aidx),lambdaval(bidx),allnrg(aidx),allnrgf(aidx),allnrg(bidx),allnrgb(bidx)
          endif
        enddo

        do i=2,nlambdas
          if (.not.updated(i)) call MPI_Send(task2lambda(i),1,MPI_INT,i-1,11,MPI_COMM_WORLD,ierr)
        enddo

      else    
         call MPI_Recv(mylambda,1,MPI_INT,0,11,MPI_COMM_WORLD,mpistat,ierr)
      endif

      call MPI_BCAST (nlretry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (nlreswap, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (DifferentLambdas) then
        delf = lambda-lambdaf
        delb = lambda-lambdab
        lambdaf = lambdaval(mylambda)-delf
        lambdab = lambdaval(mylambda)-delb
        if (lambdaf.gt.1.0) lambdaf = 1.0
        if (lambdaf.lt.0.0) lambdaf = 0.0
        if (lambdab.gt.1.0) lambdab = 1.0
        if (lambdab.lt.0.0) lambdab = 0.0
      endif
      oldlam = lambda
      lambda = lambdaval(mylambda)

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
      nrgb = nrgb + inrgb
      nrgf = nrgf + inrgf      
      call clearMoveList
      call updateExtraEnergies

c     Re-setting all running averages    
      nsteps = 0

      write(tempstr,"(1A,F5.3)")"lam-",oldlam
      write(lrefolder,"(2A,F5.3)")reoutfolder(1:lineLength(reoutfolder)),"/lam-",lambdaval(mylambda)
      do i=1,nOpenFiles    
        fileappend(i) = .True.
        fileopened(i) = .False.
        call split2(openfiles(i),"/")
        tempstr2 = words2(nwords2-1)
        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
          do j=1,nwords2-2
            words(j) = words2(j)
          enddo
           write(words(nwords2-1),"(1A,F5.3)")"lam-",lambda
           words(nwords2) = words2(nwords2)
           nwords = nwords2
           call join2(1,nwords,"/",openfiles(i))
        endif
      enddo

c      call MPI_FINALIZE(ierr)
c      call exit(0)

      end subroutine
