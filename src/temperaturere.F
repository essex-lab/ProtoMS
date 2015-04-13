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
      logical updated(ntemperatures),reTest
      double precision nrg,nrgl,nrgh,teml,temh,temp,oldtem
      double precision totalEnergy,extraEnergy
      double precision allteml(ntemperatures),alltemh(ntemperatures),allnrg(ntemperatures)
      double precision dell,delh
      double precision buffer(3)
      character*255 tempstr,tempstr2


c     check that temperature swaps are only attempted between tasks that are running temperature values
      if (mytemid.lt.0) return

      write(6,*)"\n",myid," !!!mytrefolder before tswap ",trefolder
      write(6,*)"\n",myid," !!!mylrefolder before tswap ",lrefolder

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

      if (mytemid.eq.0) then
        allnrg(mytemperature) = nrg
        updated(1) = .false.
        nrecieved = 1
        do i=2,ntemperatures
          call MPI_Recv(buffer,3,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11,mytemcom,mpistat,ierr)
          task1 = mpistat(MPI_SOURCE)
          allteml(task2temperature(task1+1)) = buffer(1)
          alltemh(task2temperature(task1+1)) = buffer(2)
          allnrg(task2temperature(task1+1)) = buffer(3)
          updated(task1+1) = .false.
        enddo
      else
        buffer(1) = teml
        buffer(2) = temh
        buffer(3) = nrg
        call MPI_SSend(buffer,3,MPI_DOUBLE_PRECISION,0,11,mytemcom,ierr)
      endif


      if (mytemid.eq.0) then

c       Randomly selects odd/even pairs
        choice = ranInt(0,1)
        do i=1,ntemperatures-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1

          nlretry = nlretry + 1
c          write(6,*)nlretry
          if (reTest(pt_beta_t(mylambda,aidx),pt_beta_t(mylambda,bidx),allnrg(aidx),allnrg(bidx))) then

c            write(6,*)"Swapping ",temval(aidx),temval(bidx),allnrgl(bidx)-allnrg(aidx),allnrgh(aidx)-allnrg(bidx)
c            write(6,*)allnrg(aidx),allnrgh(aidx),allnrg(bidx),allnrgl(bidx)
c            write(6,*)-(allnrg(bidx)-allnrgh(aidx)-allnrgl(bidx)+allnrg(aidx))
            task1 = temperature2task(aidx)
            task2 = temperature2task(bidx)

            task2temperature(task1+1) = bidx
            task2temperature(task2+1) = aidx
            temperature2task(aidx) = task2
            temperature2task(bidx) = task1

            ntreswap = ntreswap + 1
            if (task1.ne.0) then
              call MPI_SSend(bidx,1,MPI_INT,task1,11,mytemcom,ierr)
            else
              mytemperature = bidx
            endif
            if (task2.ne.0) then
              call MPI_SSend(aidx,1,MPI_INT,task2,11,mytemcom,ierr)
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
          if (.not.updated(i)) call MPI_SSend(task2temperature(i),1,MPI_INT,i-1,11,mytemcom,ierr)
        enddo

      else    
         call MPI_Recv(mytemperature,1,MPI_INT,0,11,mytemcom,mpistat,ierr)
      endif

      call MPI_BCAST (ntretry, 1, MPI_INTEGER, 0, mytemcom, ierr)
      call MPI_BCAST (ntreswap, 1, MPI_INTEGER, 0, mytemcom, ierr)

      oldtem = temperatureC
      temperatureC = temval(mytemperature)
      temperature = temperatureC + DEG2KEL
      mt_beta_t = dble(1.0) / (GASR*temperature)

c     Re-setting all running averages    
      nsteps = 0
      
      write(tempstr,"(1A,F6.3)")"t-",oldtem
      write(6,*)myid," tempstr ",tempstr
      write(trefolder,"(2A,F6.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(mytemperature)
      if (oldtem.ge.100.0) write(tempstr,"(1A,F7.3)")"t-",oldtem
      if (temval(mytemperature).ge.100.0) write(trefolder,"(2A,F7.3)")outfolder(1:lineLength(outfolder)),
     .                                    "/t-",temval(mytemperature)
      if (doLambdaRE) write(lrefolder,"(2A,F5.3)")trefolder(1:lineLength(trefolder)),"/lam-",lambdaval(mylambda)
      write(6,*)myid," myLAMBDA ",mylambda," mylambdaval ",lambdaval(mylambda)

c     Replacing the old temperature for the new one in the file name.
      k = 1
      if ((dolambdaRE).or.(doMultiGCMC).or.(doMultiJaws2)) k = 2
      do i=1,nOpenFiles    
        if (wasopened(i)) fileappend(i) = .true.
        if (i.eq.streamhandle(RESTART)-40) then
          fileappend(i) = .False.
        endif
        fileopened(i) = .False.
        call split2(openfiles(i),"/")
        tempstr2 = words2(nwords2-k)
        write(6,*)myid," tempstr2 ",tempstr2
        if (tempstr2(1:lineLength(tempstr2)).eq.tempstr(1:lineLength(tempstr))) then
          do j=1,nwords2-(k+1)
            words(j) = words2(j)
          enddo
           write(words(nwords2-k),"(1A,F6.3)")"t-",temperatureC
           if (temperatureC.ge.100.0) write(words(nwords2-k),"(1A,F7.3)")"t-",temperatureC
           do j=1,k
             words(nwords2-(j-1)) = words2(nwords2-(j-1))
           enddo
           nwords = nwords2
           call join2(1,nwords,"/",openfiles(i))
           write(6,*)myid," openfiles(i) ",openfiles(i)
        endif
      enddo

      write(6,*)"\n",myid," !!!mytrefolder after tswap ",trefolder
      write(6,*)"\n",myid," !!!mylrefolder after tswap ",lrefolder

c      call MPI_FINALIZE(ierr)
c      call exit(0)

      end subroutine
