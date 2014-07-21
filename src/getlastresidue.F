      integer function getLastResidue(FILE,line)
      implicit none
c######################################################
c
c     This function scans the pdb file FILE to
c     find which is the last residue in the protein.
c
c     This function assumes that all of the atoms 
c     in a residue are together
c
c     (C) Christopher Woods
c
c######################################################

      integer FILE
      character*(*) line
      character*6 word
      integer resnum
      
      getLastResidue = 0
      
100   continue      
        read(FILE,10,err=9700,end=9700) line
10      format(a300)

c       get the first word in the file
        word = line(1:6)
c       lowercase the word
        call lowerCase(word)

c       if the line starts with 'ter' then return as we have reached the
c       end of the protein
        if (word(1:3).eq.'ter') goto 9700

c       see if the line starts with 'atom' or 'hetatm'
        if (word(1:4).ne.'atom' .and. word(1:6).ne.'hetatm') goto 100
      
c       pdb files are column formatted - this can be quite a problem as
c       this limits the number of atoms to 9999...
        read(line,20,err=9600,end=9600) resnum
20      format(23x,i3)        

        if (resnum.ge.0) getLastResidue = resnum

9600    continue
                
      goto 100
      
9700  continue

      rewind(FILE)
      
      return
      end
