      subroutine loadgc(filename)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
c#############################################################
c
c     Based upon loadsolvent.F, a routine to load in GC-solutes
c     These are treated as rigid water molecules
c
c     MSB 31/03/09
c
c#############################################################

      character*(*) filename
      character*300 line
      integer FILE
      integer igcs,iatm,ioldres,ngcs
      integer atnum,iresnum
      character*4 atnam
      character*3 iresnam
      double precision ix,iy,iz
      logical ok
      integer lineLength,lgth,linenum
      parameter(FILE=40)
      
      lgth = lineLength(filename)
      
      write(printstring,*) "Opening GC-solute file ",filename(1:lgth)
      call printLine(INFO,printstring)
      
      open(unit=FILE,file=filename,status="old",form='formatted',err=9600)

      rewind(FILE)
      
      igcs = 0
      ngcs = 0
      iatm = 0
      ioldres = 0
      linenum = 0
      
c     read the file in a do-while loop      
100   continue
        linenum = linenum + 1
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       lowercase and split the line
        call lowercase(line)
        call split(line)
        
        if (nwords.lt.1) goto 100
        
        if (words(1)(1:3).eq.'ter') then
c         we move on to the next solute molecule        
          iatm = 0
        else
c         try to read in a GCsolute atom line
          call readPDBLine(line,ok,atnum,atnam,iresnum,iresnam,ix,iy,iz)
          
c         'ok' is true if this line contained a valid pdb atom
          if (.not.ok) then
             if (words(1).eq.'atom' .or. words(1).eq.'hetatm') then
               write(printstring,*) "Could not understand the following PDB line: ",linenum
               call printLine(WARNING,printstring)
               call printLine(WARNING,line)
             endif
            goto 100
          endif

          if (iatm.eq.0 .or. iresnum.ne.ioldres) then
c           this is the first atom of a new GCsolute.
            igcs = NGCsolutes + 1
            if (igcs.gt.MAXSOLVENTS) then
              write(printstring,*) "Exceeded maximum number of GC-solute molecules, ",MAXSOLVENTS
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXSOLVENTS and recompile")
              goto 100
            endif
            
            NGCsolutes = igcs
            ngcs = ngcs + 1
            GCSName(igcs) = iresnam
            GCSRes(igcs) = iresnum
            ioldres = iresnum
            iatm = 0
          endif

          iatm = iatm + 1
          if (iatm.gt.MAXSOLUTEATOMS) then
            write(printstring,*) "Exceeded maximum number of atoms per gc-solute molecule (",
     .                           MAXSOLUTEATOMS,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,"Increase MAXSOLUTEATOMS and recompile")
            goto 100
          endif
          
          NGCSAtoms(igcs) = iatm
          GCSAtomName(igcs,iatm) = atnam
          ZGCSCoords(igcs,iatm,1) = ix
          ZGCSCoords(igcs,iatm,2) = iy
          ZGCSCoords(igcs,iatm,3) = iz

        endif
        
        goto 100
c     end of do-while loop

9600  continue
c     get here if we could not open the file
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - unable to load GCsolute file ",filename(1:lgth)
      call printLine(WARNING,printstring)
      return

9700  continue
c     got here if there was a problem reading a line
      lgth = lineLength(filename)
      write(printstring,*) "WARNING - problem reading line ",linenum," of the GCsolute file ",
     .     filename(1:lgth)
      call printLine(WARNING,printstring)
      call printLine(WARNING,line)
      call printLine(WARNING,"Skipping the rest of the file!")
      return

c     get here when hit end of the file
9800  continue
c     -------------------------------

      write(printstring,*) "Have now loaded ",NGCsolutes," GCsolute molecules (",ngcs,
     .                     " from this file)"
      call printLine(INFO,printstring)
      
      close(FILE)

      return
            
      end
