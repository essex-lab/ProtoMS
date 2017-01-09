      subroutine gcmcRE
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'constants.inc'
      include 'enums.inc'
      include 'averages.inc'
      include 'printer.inc'
      include 'parallel.inc'
      include 'splitter.inc'
      include 'sep.inc'
      include 'mpif.h' 
      include 'metropolis.inc'
      include 'stacks.inc'
c###################################################################
c
c     This routine peforms GCMC replica exchange
c
c     Gregory Ross, 2015
c     Based on lamdbare.F and temperature.F
c
c###################################################################

      integer i,j,ranInt,task1,task2,choice,aidx,bidx,linelength,k
      integer mpistat(MPI_STATUS_SIZE),ierr,nrecieved
      logical updated(nadams), adamsretest
      integer lamsbuffer(4),lamrbuffer(4),tsbuffer(nlambdas),lmbuffer(nlambdas)
      double precision alladam(nadams),allongc(nadams),oldB
      double precision buffer(2)
      character*255 adamstr,adamstr2



c		GR: Saving the old b value, which is needed to keep track of the folders that change.
      oldB = B
c     check that temperature swaps are only attempted between tasks that are running temperature values
      if (mygcid.lt.0) return

      do i = 1,nadams
        updated(i) = .false.
      enddo

      call closeFiles


c		GR: This part defines how the information is communicated. mygcid=0 is the by convention the command process, so does not SEND information; instead, it RECIEVES information.
      if (mygcid.eq.0) then
        alladam(myadam) = adamparam(myadam)
        allongc(myadam) = GCSCavity
        updated(1) = .false.
        nrecieved = 1
        do i=2,nadams
          call MPI_Recv(buffer,2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11,mygccom,mpistat,ierr)
          task1 = mpistat(MPI_SOURCE)
          alladam(task2adams(task1+1)) = buffer(1)
          allongc(task2adams(task1+1)) = buffer(2)
          updated(task1+1) = .false.
        enddo
      else
        buffer(1) = adamparam(myadam)
        buffer(2) = GCSCavity
        call MPI_SSend(buffer,2,MPI_DOUBLE_PRECISION,0,11,mygccom,ierr)
      endif



      if (mygcid.eq.0) then
c        write(6,*)myid,mygcid, "GR gcmcre.F Trying adam swap" 
c       Randomly selects odd/even pairs from 1 to the max number of Adams values and attempts to swap those pairs.
        choice = ranInt(0,1)
        do i=1,nadams-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1

          ngcretry = ngcretry + 1
c			GR: The Metropolis test is below. THIS IS THE IMPORTANT BIT!
          if (adamsretest(alladam(aidx),alladam(bidx),allongc(aidx),allongc(bidx))) then
            task1 = adams2task(aidx)
            task2 = adams2task(bidx)

            task2adams(task1+1) = bidx
            task2adams(task2+1) = aidx
            adams2task(aidx) = task2
            adams2task(bidx) = task1
            ngcreswap = ngcreswap + 1
c           The main task already knows bidx and aidx, so no need to send that information in that case
c            write(6,*)myid, mygcid,"GR gcmcre.F Sending GC swap info to other gcmc nodes"
            if (task1.ne.0) then 
              call MPI_SSend(bidx,1,MPI_INTEGER,task1,11,mygccom,ierr)
            else
              myadam = bidx
            endif
            if (task2.ne.0) then
              call MPI_SSend(aidx,1,MPI_INTEGER,task2,11,mygccom,ierr)
            else
              myadam = aidx
            endif
            updated(task1+1) = .true.
            updated(task2+1) = .true.
          else
          endif
        enddo
        do i=2,nadams
          if (.not.updated(i)) call MPI_SSend(task2adams(i),1,MPI_INTEGER,i-1,11,mygccom,ierr)
        enddo
      else    
         call MPI_Recv(myadam,1,MPI_INTEGER,0,11,mygccom,mpistat,ierr)
      endif

      call MPI_BCAST (ngcretry, 1, MPI_INTEGER, 0, mygccom, ierr)
      call MPI_BCAST (ngcreswap, 1, MPI_INTEGER, 0, mygccom, ierr)

      call MPI_BCAST (updated,nadams,MPI_LOGICAL,0,mygccom,ierr)
      call MPI_BCAST (adams2task,nadams,MPI_INTEGER,0,mygccom,ierr)
      call MPI_BCAST (task2adams,nadams,MPI_INTEGER,0,mygccom,ierr)

c		GR: Set the new adams value as B, which is used in Metropolis2.F
      B = adamparam(myadam)

c     The communication in this code is not optimized. Task2temperature, temperature2task, nlretry & nlreswap probably only need to be swaped
c     if a process with temid 0 is involved in the swap, since the reason to swap these is problems associated with the next temperatureswap

      if (updated(mygcid+1).and.doLambdaRE) then
        lamsbuffer(1) = mylamladder
        lamsbuffer(2) = mylamid
        lamsbuffer(3) = ngcretry
        lamsbuffer(4) = ngcreswap
        do i=1,nadams
          if (oldb.eq.adamparam(i)) then 
            call MPI_Sendrecv(lamsbuffer,4,MPI_INTEGER,adams2task(i),11,lamrbuffer,4,
     .           MPI_INTEGER,adams2task(i),11,mygccom,mpistat,ierr)
            call MPI_Sendrecv(task2lambda,nlambdas,MPI_INTEGER,adams2task(i),11,tsbuffer,nlambdas,
     .           MPI_INTEGER,adams2task(i),11,mygccom,mpistat,ierr)
            call MPI_Sendrecv(lambda2task,nlambdas,MPI_INTEGER,adams2task(i),11,lmbuffer,nlambdas,
     .           MPI_INTEGER,adams2task(i),11,mygccom,mpistat,ierr)
          endif
        enddo
        task2lambda = tsbuffer
        lambda2task = lmbuffer
        mylamladder = lamrbuffer(1)
        mylamid = lamrbuffer(2)
        ngcretry = lamrbuffer(3)
        ngcreswap = lamrbuffer(4)
      endif

c     Resetting running averages
      nsteps = 0

      if (doLambdaRE) then
c     Defining gcmcfolder as the the folder that correponds to the NEW b value, using adamstr as a dummy variable.
        write(adamstr,"(4A)")lrefolder(1:lineLength(lrefolder)),PATHSEP,"b_"
        call writedoublestring(adamstr,adamparam(myadam),gcmcfolder)
      else 
c     Defining gcmcfolder as the the folder that correponds to the NEW b value, using adamstr as a dummy variable.
        write(adamstr,"(4A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
        call writedoublestring(adamstr,adamparam(myadam),gcmcfolder)
c     Defining adamstr as the the folder that correponds to the OLD b value.
      endif
      call writedoublestring("b_",oldB,adamstr)

c      if (doLambdaRE) write(lrefolder,"(3A,F5.3)")gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP,"lam-",lambdaval(mylambda)

c		GR: Must account for the fact that there must be nested folders if replica exchange is attempted along another dimention.
      k = 1
c      if ((dolambdaRE).or.(doREST).or.(doMultiJaws2).or.(doTemperatureRE)) k = 1
      do i=1,nOpenFiles    
        if (wasopened(i)) fileappend(i) = .true.
        if (i.eq.streamhandle(RESTART)-40) then
          fileappend(i) = .False.
        endif
        fileopened(i) = .False.
        call split2(openfiles(i),PATHSEP)
c			GR: By default, the temperature is the first folder, so the -k bit always picks that out. 
        adamstr2 = words2(nwords2-k)
c			GR: This next bit sees which of the currently opened folders corresponds to the old adams value.
        if (adamstr2(1:lineLength(adamstr2)).eq.adamstr(1:lineLength(adamstr))) then
          do j=1,nwords2-(k+1)
            words(j) = words2(j)
          enddo
          call writedoublestring("b_",adamparam(myadam),words(nwords2-k))
          do j=1,k
            words(nwords2-(j-1)) = words2(nwords2-(j-1))
          enddo
          nwords = nwords2
          call join2(1,nwords,PATHSEP,openfiles(i))
        endif
      enddo



      end subroutine
