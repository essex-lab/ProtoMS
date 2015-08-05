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
      double precision alladam(nadams),allongc(nadams),oldB
      double precision buffer(2)
      character*255 adamstr,adamstr2

c 		GR: Assign to double pression a vector alladams(nadams), where nadams is defined elsewhere.
c		GR: Also, I need to assign allongc(nadams), which is the vector of all on GC solutes

c		GR: Not sure why I need this, maybe something to do with the folders?
      oldB = B
c     check that temperature swaps are only attempted between tasks that are running temperature values
      if (myid.lt.0) return

      do i = 1,nadams
        updated(i) = .false.
      enddo

      call closeFiles

c		GR: Need to replace temperature vectors with vector of adams values. teml --> adaml, temh --> adamh etc. DONE
c		GR: myadam is the temperature of the process. Where is it defined? 
c		GR: NGCSolutes is the number of "on" GC solutes.
c		GR: Make substutions: adaml = adamparam(myadam-1), adamh = adamparam(myadam+1), etc. DONE

c		GR: Don't think I need this bit. See the exchange test below.
c      if (myadam.eq.1) then
c        adaml = adamparam(myadam)
c      else
c        adaml = adamparam(myadam-1)
c      endif

c      if (myadam.eq.nadams) then
c        adamh = adamparam(myadam)
c      else
c        adamh = adamparam(myadam+1)
c      endif

c		This sets everything (residues, solvents, gcsolute etc) as having been moved. It is so the extra energies are calculated correctly. Not necessary.
c      call setAllMoved

c		GR: Delete the below when certain.
c      if (.not.doREST) then
c       nrg = totalEnergy(DONTCARE,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)
c        testbeta = mt_beta_t
c      else
c        nrg = restTestEnergy()
c        testbeta = rest_beta      
c      endif


c		GR: This part defines how the information is communicated. myid=0 is the by convention the command process, so does not SEND information; instead, it RECIEVES information.
c		GR: I require the equivalent of mytemcom
c		GR: Check that GCSCavity is properly read by this file.
      if (myid.eq.0) then
c        alladaml(myadam) = adaml
c        alladamh(myadam) = adamh
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

c		GR: I've started modifying the below.
      if (myid.eq.0) then
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
            if (task1.ne.0) then 
              call MPI_SSend(bidx,1,MPI_INT,task1,11,mygccom,ierr)
            else
              myadam = bidx
            endif
            if (task2.ne.0) then
              call MPI_SSend(aidx,1,MPI_INT,task2,11,mygccom,ierr)
            else
              myadam = aidx
            endif
            updated(task1+1) = .true.
            updated(task2+1) = .true.
          else
c            write(6,*)"Failing ",adamparam(aidx),adamparam(bidx),allnrg(aidx),allnrgh(aidx),allnrg(bidx),allnrgl(bidx)
          endif
        enddo

        do i=2,nadams
          if (.not.updated(i)) call MPI_SSend(task2adams(i),1,MPI_INT,i-1,11,mygccom,ierr)
        enddo

      else    
         call MPI_Recv(myadam,1,MPI_INT,0,11,mygccom,mpistat,ierr)
      endif

      call MPI_BCAST (ngcretry, 1, MPI_INTEGER, 0, mygccom, ierr)
      call MPI_BCAST (ngcreswap, 1, MPI_INTEGER, 0, mygccom, ierr)
      
      call MPI_BCAST (updated,nadams,MPI_LOGICAL,0,mygccom,ierr)
      call MPI_BCAST (adams2task,nadams,MPI_INTEGER,0,mygccom,ierr)
      call MPI_BCAST (task2adams,nadams,MPI_INTEGER,0,mygccom,ierr)

c		GR: Set the new adams value as B, which is used in Metropolis2.F
      B = adamparam(myadam)
      write(6,*)"### GR: old B, new B = ",oldB, B

c		GR: Do I need the below? What does it do? If it pertains to energies then it may not be needed.      
c      if (.not.doREST) then
c        oldtem = temperatureC
c        temperatureC = adamparam(myadam)
c        temperature = temperatureC + DEG2KEL
c        mt_beta_t = dble(1.0) / (GASR*temperature) 
c      else 
c        oldtem = resttempC
c        resttempC = adamparam(myadam)
c        resttemp = resttempC + DEG2KEL
c        rest_beta = dble(1.0) / (GASR*resttemp)
c        restfac1 = rest_beta / mt_beta_t
c        restfac2 = (mt_beta_t+rest_beta)/(2.0d0*mt_beta_t)
c        do i=1,7
c          if (restlevel(i).eq.1) then
c            restscaling(i) = restfac1
c          else if (restlevel(i).eq.2) then
c            restscaling(i) = restfac2
c          else
c            restscaling(i) = 1.0d0
c          endif
c        enddo
c        call setAllMoved
c        call initialiseExtraEnergies(.true.)      
c        nrg = totalEnergy(OLD,nrgb,nrgf) + extraEnergy(NEW,inrgb,inrgf)    
c        call clearMoveList
c        call updateExtraEnergies
c      endif

c		GR: May only be needed for combining GCMC with Lambde replica exchange. I'll come to this later.
c      if (updated(myid+1).and.doLambdaRE) then
c        lamsbuffer(1) = mylamladder
c        lamsbuffer(2) = mylamid
c        lamsbuffer(3) = ntretry
c        lamsbuffer(4) = ntreswap
c        do i=1,nadams
c          if (oldtem.eq.adamparam(i)) then 
c            call MPI_Sendrecv(lamsbuffer,4,MPI_INT,temperature2task(i),11,lamrbuffer,4,
c     .           MPI_INT,temperature2task(i),11,mytemcom,mpistat,ierr)
c            call MPI_Sendrecv(task2lambda,nlambdas,MPI_INT,temperature2task(i),11,tsbuffer,nlambdas,
c     .           MPI_INT,temperature2task(i),11,mytemcom,mpistat,ierr)
c            call MPI_Sendrecv(lambda2task,nlambdas,MPI_INT,temperature2task(i),11,lmbuffer,nlambdas,
c     .           MPI_INT,temperature2task(i),11,mytemcom,mpistat,ierr)
c          endif
c        enddo
c        task2lambda = tsbuffer
c       lambda2task = lmbuffer
c        mylamladder = lamrbuffer(1)
c        mylamid = lamrbuffer(2)
c        ntretry = lamrbuffer(3)
c        ntreswap = lamrbuffer(4)
c      endif


c     Re-setting all running averages    
      nsteps = 0

c		GR: This first line assigns adamstr as a string with 6 characters: 2 before the decimal place, the decimal place, and 3 decimals.
    
c     Defining gcmcfolder as the the folder that correponds to the NEW b value, using adamstr as a dummy variable.
c      write(adamstr,"(4A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
c      call writedoublestring(adamstr,adamparam(myadam),gcmcfolder)

c     Defining adamstr as the the folder that correponds to the OLD b value, using adamstr2 as a dummy variable.
c      write(adamstr2,"(4A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
c      call writedoublestring(adamstr,oldB,adamstr2)

c     Here, files that are open by current processor must have their output directory changed to the directory corresponding to the new B value.
c        do i=1,nOpenFiles    
c        if (wasopened(i)) fileappend(i) = .true.
c        if (i.eq.streamhandle(RESTART)-40) then
c          fileappend(i) = .False.
c        endif
c        fileopened(i) = .False.
c        call split2(openfiles(i),PATHSEP)
c        tempstr2 = words2(nwords2-1)
c        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
c          do j=1,nwords2-2
c            words(j) = words2(j)
c          enddo
c           write(words(nwords2-1),"(1A,F5.3)")"lam-",lambda
c           words(nwords2) = words2(nwords2)
c           nwords = nwords2
c           call join2(1,nwords,PATHSEP,openfiles(i))
c        endif
c      enddo
c     Here, files that are open by current processor must have their output directory changed to the directory corresponding to the new B value.
c      write(tempstr,"(1A,F5.3)")"lam-",oldlam
c      if (doTemperatureRE)then
c        write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lambdaval(mylambda)
c      else
c        write(lrefolder,"(3A,F5.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"lam-",lambdaval(mylambda)
c      endif
c      do i=1,nOpenFiles    
c        if (wasopened(i)) fileappend(i) = .true.
c        if (i.eq.streamhandle(RESTART)-40) then
c          fileappend(i) = .False.
c        endif
c        fileopened(i) = .False.
c        call split2(openfiles(i),PATHSEP)
c        tempstr2 = words2(nwords2-1)
c        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
c          do j=1,nwords2-2
c            words(j) = words2(j)
c          enddo
c           write(words(nwords2-1),"(1A,F5.3)")"lam-",lambda
c           words(nwords2) = words2(nwords2)
c           nwords = nwords2
c           call join2(1,nwords,PATHSEP,openfiles(i))
c        endif
c      enddo



c     ###### GR: From this point, the old code is used.
c      write(adamstr,"(1A,F6.3)")"b_",oldB
c      write(gcmcfolder,"(3A,F6.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_",adamparam(myadam)
c      if (oldB.ge.100.0) write(adamstr,"(1A,F7.3)")"b_",oldB
c      if (adamparam(myadam).ge.100.0) write(gcmcfolder,"(3A,F7.3)")outfolder(1:lineLength(outfolder)),
c     .                                    PATHSEP,"b_",adamparam(myadam)

c     Defining gcmcfolder as the the folder that correponds to the NEW b value, using adamstr as a dummy variable.
      write(adamstr,"(4A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
      call writedoublestring(adamstr,adamparam(myadam),gcmcfolder)
      write(6,*)"### GR gcmcre.F: gcmcfolder=",gcmcfolder

c     Defining adamstr as the the folder that correponds to the OLD b value.
      call writedoublestring("b_",oldB,adamstr)
      write(6,*)"### GR gcmcre.F: adamstr=",adamstr


      if (doLambdaRE) write(lrefolder,"(3A,F5.3)")gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP,"lam-",lambdaval(mylambda)


c		GR: Must account for the fact that there must be nested folders if replica exchange is attempted along another dimention.
      k = 1
      if ((dolambdaRE).or.(doREST).or.(doMultiJaws2).or.(doTemperatureRE)) k = 2
      do i=1,nOpenFiles    
        if (wasopened(i)) fileappend(i) = .true.
        if (i.eq.streamhandle(RESTART)-40) then
          fileappend(i) = .False.
        endif
        fileopened(i) = .False.
        call split2(openfiles(i),PATHSEP)
c			GR: By default, the temperature is the first folder, so the -k bit always picks that out. 
        adamstr2 = words2(nwords2-k)
        write(6,*)"### GR gcmcre.F: adamsrt2=",adamstr2
c			GR: This next bit sees which of the currently opened folders corresponds to the old adams value.
        if (adamstr2(1:lineLength(adamstr2)).eq.adamstr(1:lineLength(adamstr))) then
          write(6,*)"### GR adamstr2 is the same as adamstr"     
          do j=1,nwords2-(k+1)
            words(j) = words2(j)
          enddo
          write(6,*)"### GR gcmcre.F: words(nwords2-k)=",words(nwords2-k)
c           write(words(nwords2-k),"(1A,F6.3)")"b_",adamparam(myadam)
c           if (adamparam(myadam).ge.100.0) write(words(nwords2-k),"(1A,F7.3)")"b_",adamparam(myadam)
          call writedoublestring("b_",adamparam(myadam),words(nwords2-k))
          write(6,*)"### GR gcmcre.F: NEW words(nwords2-k)=",words(nwords2-k)
          do j=1,k
            words(nwords2-(j-1)) = words2(nwords2-(j-1))
          enddo
          nwords = nwords2
          call join2(1,nwords,PATHSEP,openfiles(i))
        endif
      enddo

      end subroutine
