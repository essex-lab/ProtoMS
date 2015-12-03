      subroutine assignadamre(key,freeCom,mkdirs)
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
c    by Gregory A. Ross and Ana I. Cabedo Martinez (2015)
c
c##########################################################################


      integer key,colour,linelength
      integer mpistat(MPI_STATUS_SIZE),ierr,mygroup
      logical freeCom,mkdirs
      character*255 cmd, tmpstr

      write(6,*)myid,"GR assignadamre:",mylamid,mylamladder,mygcid,mygccom,mygcladder
      if (freeCom.and.(mygccom.ne.MPI_COMM_NULL)) call MPI_Comm_free(mygccom,ierr)

      call MPI_COMM_SPLIT(MPI_COMM_WORLD,mygcladder,key,mygccom,ierr)
      write(6,*)"GR assignadamre: MPI_COMM_SPLIT called"
      call MPI_COMM_RANK(mygccom,mygcid,ierr)
      write(6,*)"GR assignadamre: MPI_COMM_RANK called"
      if (mkdirs) then
        myadam = mygcid + 1
        if (doTemperatureRE) then
          if (mytemid.lt.0) then
            call findDoubleInList(temwithlad(mygcladder),temval,ntemperatures,mytemperature)
            write(trefolder,"(3A,F6.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"t-",temval(mytemperature)
            if (temval(mytemperature).ge.100.0) write(trefolder,"(3A,F7.3)")outfolder(1:lineLength(outfolder)),
     .                                           PATHSEP,"t-",temval(mytemperature)
c			GCMC folder needs to be formatted correctly
            write(tmpstr,"(4A)")"mkdir -p ",trefolder(1:lineLength(trefolder)),PATHSEP,"b_"
            call writedoublestring(tmpstr,adamparam(myadam),cmd)
            write(tmpstr,"(3A)")trefolder(1:lineLength(trefolder)),PATHSEP,"b_"
            call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
c            write(gcmcfolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"b_",adamparam(myadam)
c            write(cmd,"(3A)")"mkdir -p ",gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP
            call system(cmd) 
          else
c			GCMC folder needs to be formatted correctly
            write(tmpstr,"(4A)")"mkdir -p ",trefolder(1:lineLength(trefolder)),PATHSEP,"b_"
            call writedoublestring(tmpstr,adamparam(myadam),cmd)
            write(tmpstr,"(3A)")trefolder(1:lineLength(trefolder)),PATHSEP,"b_"
            call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
c            write(gcmcfolder,"(3A,F5.3)")trefolder(1:lineLength(trefolder)),PATHSEP,"b_",adamparam(myadam)
c            write(cmd,"(3A)")"mkdir -p ",gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP
            call system(cmd)
          endif
        else if (doLambdaRE) then
c			GCMC folder needs to be formatted correctly
          write(tmpstr,"(4A)")"mkdir -p ",lrefolder(1:lineLength(lrefolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(myadam),cmd)
          write(tmpstr,"(3A)")lrefolder(1:lineLength(lrefolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
c          write(gcmcfolder,"(3A,F5.3)")lrefolder(1:lineLength(lrefolder)),PATHSEP,"b_",adamparam(myadam)
c          write(cmd,"(3A)")"mkdir -p ",gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP
          call system(cmd)        
        else
          write(6,*)"GR assignadamre: Goes here when only doing re gcmc"
c			GCMC folder needs to be formatted correctly
          write(tmpstr,"(4A)")"mkdir -p ",outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(myadam),cmd)
          write(tmpstr,"(3A)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_"
          call writedoublestring(tmpstr,adamparam(myadam),gcmcfolder)
c          write(gcmcfolder,"(3A,F5.3)")outfolder(1:lineLength(outfolder)),PATHSEP,"b_",adamparam(myadam)
c          write(cmd,"(3A)")"mkdir -p ",gcmcfolder(1:lineLength(gcmcfolder)),PATHSEP
          call system(cmd)
        endif 
      endif

      end subroutine


