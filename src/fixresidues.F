      subroutine fixresidues(line,check)

      include 'dimensions.inc'
      include 'splitter.inc'
      include 'proteins.inc'
      include 'printer.inc'

      integer isol,iprnt,iprot,ires,iend,istrt,j,i,inflx
      integer getInt
      logical check
      character*300 line

      call split(line)

          if (words(1).eq.'fixbackbone') then
c           fix the protein backbone, format is
c           'fixbackbone iprot all none 1 3 5-10 35-46 15 #comment'
            if (nwords.lt.3) then 
              call printLine(WARNING,"Could not interpret fixbackbone line.")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif

c           get the protein involved                        
            iprot = getInt(words(2),0)
            if (iprot.le.0 .or. iprot.gt.NProteins) then
              write(printstring,"(a,i8)") "Attempting to fix backbone of invalid protein, ",iprot
              call printLine(WARNING,printstring)
            endif
            
c           now loop over all of the remaining words and get the residues involved
            do i=3,nwords
              call getRange(words(i),1,NResidues(iprot),istrt,iend)
              if (istrt.eq.0 .and. iend.eq.0) then
                write(printstring,"(a,i8)") "Moving backbone of all residues of protein ",iprot
                call printLine(INFO,printstring)
                if (.not.check) then
                  do j=1,NResidues(iprot)
                    fixedBac(iprot,j) = .false.
                  enddo
                endif
              else if (istrt.eq.iend .and. istrt.gt.0) then
                write(printstring,"(2(a,i8))") "Fixing backbone of protein ",iprot," residue ",istrt
                call printLine(INFO,printstring)
                fixedBac(iprot,istrt) = .true.
              else if (istrt.gt.0 .and. iend.gt.0) then
                write(printstring,"(3(a,i8))") "Fixing backbone of protein ",iprot," residues: ",istrt," to ",iend
                call printLine(INFO,printstring)
                if (.not.check) then
                  do j=istrt,iend
                    fixedBac(iprot,j) = .true.
                  enddo
                endif
              endif
            enddo

          else if (words(1).eq.'fixresidues') then
c           fix the protein residues, format is
c           'fixresidues iprot all none 1 3 5-10 35-46 15 #comment'
            if (nwords.lt.3) then 
              call printLine(WARNING,"Could not interpret fixresidues line.")
              call join(1,nwords,printstring)
              call printLine(WARNING,printstring)
            endif

c           get the protein involved                        
            iprot = getInt(words(2),0)
            if (iprot.le.0 .or. iprot.gt.NProteins) then
              write(printstring,"(a,i8)") "Attempting to fix residues of invalid protein, ",iprot
              call printLine(WARNING,printstring)
            endif
            
c           now loop over all of the remaining words and get the residues involved
            do i=3,nwords
              call getRange(words(i),1,NResidues(iprot),istrt,iend)
            
              if (istrt.eq.0 .and. iend.eq.0) then
                write(printstring,"(a,i8)") "Unfixing all residues of protein ",iprot
                call printLine(INFO,printstring)
                if (.not.check) then
                  do j=1,NResidues(iprot)
                    fixedRes(iprot,j) = .false.
                  enddo
                endif
              else if (istrt.eq.iend .and. istrt.gt.0) then
                write(printstring,"(2(a,i8))") "Fixing protein ",iprot," residue ",istrt
                call printLine(INFO,printstring)
                fixedRes(iprot,istrt) = .true.
              else if (istrt.gt.0 .and. iend.gt.0) then
                write(printstring,"(3(a,i8))") "Fixing protein ",iprot," residues: ",istrt," to ",iend
                call printLine(INFO,printstring)
                if (.not.check) then
                  do j=istrt,iend
                    fixedRes(iprot,j) = .true.
                  enddo
                endif
              endif
            enddo

c           now work out how many flexible residues there are in this protein
            inflx = 0
            do i=1,NResidues(iprot)
              if (.not.fixedRes(iprot,i)) then 
                inflx = inflx + 1
                ResMove(iprot,inflx) = i
              endif
            enddo
            NResMove(iprot) = inflx
          endif


      return

      end
