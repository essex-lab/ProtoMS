      subroutine tempSwap(iter,stat)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine attempts temperature RE swap moves
c
c     Samuel Genheden, 2013
c
c###################################################################
      integer iter,lami,i,ranInt,choice,aidx,bidx,stat,lineLength
      double precision temp,nrg(MAXTEMPS)
      character*255 cmd
      logical reTest

      stat = 0
      do lami=1,nlambda
    
c       Read all energies	
        if (ntemp(lami).gt.1) then
          do i=1,ntemp(lami)
            call readSPEnergy(energyfiles(lami,i),nrg(i),stat)
            if (stat.ne.0) return  
          enddo        
          write(enediststream(lami),"(I6,50F15.4)")iter,(nrg(i),i=1,ntemp(lami))
        endif

c       Randomly selects odd/even pairs
        choice = ranInt(0,1)
        do i=1,ntemp(lami)-(choice+1),2
          aidx = i + choice
          bidx = i + choice + 1
        
c          if (lami.eq.1) then
c            write(6,*)aidx," ",bidx," ",nrgA," ",nrgB
c          endif
c          write(printstring,"(A,2I6,X,2F5.2)")'Attempting swap of',aidx,bidx,lambdaval(aidx),lambdaval(bidx)
c          call prettyline(INFO,printstring)

c         Perform Metropolis test
          npttry(lami) = npttry(lami) + 1
          if (reTest(pt_beta_t(lami,aidx),pt_beta_t(lami,bidx),nrg(aidx),nrg(bidx))) then
c           Swap restart files
            write(cmd,"(4A)")
     .    "cp ",rstfiles(lami,aidx)(1:lineLength(rstfiles(lami,aidx)))," ",outdir(1:lineLength(outdir))//"/RETITEMP"
            call system(cmd)
            write(cmd,"(4A)")
     .   "mv ",rstfiles(lami,bidx)(1:lineLength(rstfiles(lami,bidx)))," ",rstfiles(lami,aidx)(1:lineLength(rstfiles(lami,aidx)))
            call system(cmd)
            write(cmd,"(4A)")
     .  "mv ",outdir(1:lineLength(outdir))//"/RETITEMP"," ",rstfiles(lami,bidx)(1:lineLength(rstfiles(lami,bidx)))
            call system(cmd)
          
            temp = pttempcurr(lami,aidx)
            pttempcurr(lami,aidx) = pttempcurr(lami,bidx)
            pttempcurr(lami,bidx) = temp
            nptswap(lami) = nptswap(lami) + 1
c            write(printstring,*)'Test passed'
c            call prettyline(INFO,printstring)
c            write(printstring,*)'Swapping: ',rstfiles(aidx)(1:lineLength(rstfiles(aidx))),' ',rstfiles(bidx)(1:lineLength(rstfiles(bidx)))
c            call prettyline(INFO,printstring)
c          else
c            write(printstring,*)'Test failed'
c            call prettyline(INFO,printstring)
          endif
        enddo
      enddo

      end subroutine
