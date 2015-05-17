      subroutine assigntemperaturere(key,freeCom,mkdirs)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'parallel.inc'
      include 'flags.inc'
      include 'printer.inc'
      include 'mpif.h'

c##########################################################################
c
c    This routine re-sets the correct groups and communicators for lambdare
c    after a temperaturere swap has been performed
c
c    by Ana I. Cabedo Martineez (2015)
c
c##########################################################################


      integer key,linelength
      integer mpistat(MPI_STATUS_SIZE),ierr,mygroup
      logical freeCom,mkdirs
      character*255 cmd

      if (freeCom.and.(mytemcom.ne.MPI_COMM_NULL)) call MPI_Comm_free(mytemcom,ierr)

      if (mytemladder.eq.-1) then
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,MPI_UNDEFINED,key,mytemcom,ierr)
c          call MPI_COMM_RANK(mytemcom,mytemid,ierr)
      else
        call MPI_COMM_SPLIT(MPI_COMM_WORLD,mytemladder,key,mytemcom,ierr)
        call MPI_COMM_RANK(mytemcom,mytemid,ierr)
        if (mkdirs) then
          mytemperature = mytemid + 1
          write(trefolder,"(2A,F6.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(mytemperature)
          if (temval(mytemperature).ge.100.0) write(trefolder,"(2A,F7.3)")outfolder(1:lineLength(outfolder)),"/t-",temval(mytemperature)
          if (mytemladder.eq.1) then
            write(cmd,"(3A)")"mkdir -p ",trefolder(1:lineLength(trefolder)),'/'
            call system(cmd) 
          endif
        endif
      endif

      end subroutine


