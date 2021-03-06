      subroutine assignlambdare(key,freeCom,mkdirs)
      include 'dimensions.inc'
      include 'parallel.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'sep.inc'
      include 'mpif.h'

c##########################################################################
c
c    This routine sets the correct groups and communicators for lambdare
c    either initially or after a temperaturere swap has been performed.
c
c    by Ana I. Cabedo Martinez (2015)
c
c##########################################################################


      integer key,colour,linelength
      integer mpistat(MPI_STATUS_SIZE),ierr,mygroup
      logical freeCom,mkdirs
      character*255 cmd

      if (freeCom.and.(mylamcom.ne.MPI_COMM_NULL)) call MPI_Comm_free(mylamcom,ierr)


      if (mylamladder.eq.-1) then
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,key,mylamcom,ierr)
        if (doTemperatureRE.and.mytemid.ne.-1.and.mkdirs) then
          call findDoubleInList(lamwithlad(mytemladder),lambdaval,nlambdas,mylambda)
          write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lamwithlad(mytemladder)
          write(cmd,"(3A)")"mkdir -p ",lrefolder(1:lineLength(lrefolder)),PATHSEP
          call system(cmd)
        endif
      else
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,mylamladder,key,mylamcom,ierr)
        call MPI_COMM_RANK(mylamcom,mylamid,ierr)
        if (mkdirs) then
          mylambda = mylamid + 1
          if (doTemperatureRE) then
            if (mytemid.lt.0) then
              call findDoubleInList(temwithlad(mylamladder),temval,ntemperatures,mytemperature)
              write(trefolder,"(3A,F6.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"t-",temval(mytemperature)
              if (temval(mytemperature).ge.1000.0) then
                write(trefolder,"(3A,F8.3)")outfolder(1:lineLength(outfolder)),
     .                                             PATHSEP,"t-",temval(mytemperature)
              else if (temval(mytemperature).ge.100.0) then
                write(trefolder,"(3A,F7.3)")outfolder(1:lineLength(outfolder)),
     .                                             PATHSEP,"t-",temval(mytemperature)
              endif
              write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lambdaval(mylambda)
              write(cmd,"(3A)")"mkdir -p ",lrefolder(1:lineLength(lrefolder)),PATHSEP
              call system(cmd) 
            else
              write(lrefolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"lam-",lambdaval(mylambda)
              write(cmd,"(3A)")"mkdir -p ",lrefolder(1:lineLength(lrefolder)),PATHSEP
              call system(cmd)
            endif
          else
            write(lrefolder,"(3A,F5.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"lam-",lambdaval(mylambda)
            write(cmd,"(3A)")"mkdir -p ",lrefolder(1:lineLength(lrefolder)),PATHSEP
            call system(cmd)
          endif 
        endif
      endif

      end subroutine


