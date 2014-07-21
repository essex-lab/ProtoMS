      subroutine retiSwap(stat)
      include 'dimensions.inc'
      include 'reti.inc'
      include 'printer.inc'
c###################################################################
c
c     This routine attempts lambda RE swap moves
c
c     Samuel Genheden, 2012
c
c###################################################################
      integer i,ranInt,choice,aidx,bidx,stat,lineLength
      double precision lamA,lamfA,lambA,nrgA,nrgfA,nrgbA
      double precision lamB,lamfB,lambB,nrgB,nrgfB,nrgbB,temp
      character*255 cmd
      logical retiTest
     
      stat = 0
c     Randomly selects odd/even pairs
      choice = ranInt(0,1)
      do i=1,nlambda-(choice+1),2
        aidx = i + choice
        bidx = i + choice + 1
        
c       Read RETI energy files  
        call readRETIEnergy(retifiles(aidx),lamA,lamfA,lambA,nrgA,nrgfA,nrgbA,stat)
        if (stat.ne.0) return 
        
        call readRETIEnergy(retifiles(bidx),lamB,lamfB,lambB,nrgB,nrgfB,nrgbB,stat)
        if (stat.ne.0) return
        
c        write(printstring,"(A,2I6,X,2F5.2)")'Attempting swap of',aidx,bidx,lambdaval(aidx),lambdaval(bidx)
c        call prettyline(INFO,printstring)

c       Perform Metropolis test
        nretitry = nretitry + 1
        if (retiTest(lamA,lamfA,lambA,nrgA,nrgfA,nrgbA,lamB,lamfB,lambB,nrgB,nrgfB,nrgbB)) then
c         Swap restart files
          write(cmd,"(4A)")"cp ",rstfiles(aidx,1)(1:lineLength(rstfiles(aidx,1)))," ",outdir(1:lineLength(outdir))//"/RETITEMP"
          call system(cmd)
          write(cmd,"(4A)")
     .  "mv ",rstfiles(bidx,1)(1:lineLength(rstfiles(bidx,1)))," ",rstfiles(aidx,1)(1:lineLength(rstfiles(aidx,1)))
          call system(cmd)
          write(cmd,"(4A)")"mv ",outdir(1:lineLength(outdir))//"/RETITEMP"," ",rstfiles(bidx,1)(1:lineLength(rstfiles(bidx,1)))
          call system(cmd)
          
          temp = lambdacurr(aidx)
          lambdacurr(aidx) = lambdacurr(bidx)
          lambdacurr(bidx) = temp
          nretiswap = nretiswap + 1
c          write(printstring,*)'Test passed'
c          call prettyline(INFO,printstring)
c          write(printstring,*)'Swapping: ',rstfiles(aidx)(1:lineLength(rstfiles(aidx))),' ',rstfiles(bidx)(1:lineLength(rstfiles(bidx)))
c          call prettyline(INFO,printstring)
c        else
c          write(printstring,*)'Test failed'
c          call prettyline(INFO,printstring)
        endif
      enddo
      
      end subroutine

