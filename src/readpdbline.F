      subroutine readPDBLine(line,ok,atnum,atnam,resnum,resnam,ix,iy,iz)
      implicit none
      include 'printer.inc'
c####################################################################
c
c     This routine tries to read in a PDB atom from 'line'. It returns
c     'ok' as .true. if it succeeded. If 'ok' is .true., then the
c     atom number (atnum), name (atnam), residue number (resnum), 
c     name (resnam) and cartesian coordinates (ix,iy,iz) are
c     also correctly returned
c
c     (C) Christopher Woods, October 2004
c
c####################################################################

      character*(*) line
      logical ok
      integer atnum,resnum
      character*4 atnam,resnam
      double precision ix,iy,iz
      integer lgth,lineLength
      
      ok = .false.
      
c     make sure that the line is long enough to read  
      lgth = lineLength(line)
      if (lgth.lt.54) return
      
c     see if the line starts with 'atom' or 'hetatm'
      if (line(1:4).ne.'atom' .and. line(1:6).ne.'hetatm') return
      
c     pdb files are column formatted - this can be quite a problem as
c     this limits the number of atoms to 9999...
      read(line,10,err=9700,end=9800) atnum,atnam,resnam,resnum,ix,iy,iz
10    format(6x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3)

c     now fix the atom and residue names so that they are left justified
      call leftJustify(atnam)
      call leftJustify(resnam)
      
c     make sure that the atom and residue names are valid
      if (atnam.eq.'   ' .or. resnam.eq.'   ') then
        ok = .false.
      else
        ok = .true.
      endif

      return

9700  continue
      call printLine(WARNING,"Error reading PDB line")
      call printLine(WARNING,line)
      return
      
9800  continue
      call printLine(WARNING,"Reached end of line while trying to read PDB data?")
      call printLine(WARNING,line)
      return

      end
