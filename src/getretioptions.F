      subroutine getRETIOptions
      include 'dimensions.inc'
      include 'reti.inc'
      include 'keyset.inc'
      include 'splitter.inc'
      include 'constants.inc'
      include 'printer.inc'
c#######################################################
c
c     This subroutine gets the values of most of the
c     user control options. The options are sets of
c     key-value pairs.
c
c     The values are set either in a command file
c     (see parsefile.F) or via environmental 
c     variables (see subroutine getenv2).
c
c     Samuel Genheden, November 2012
c     From getoptions.F
c
c#######################################################

      integer getInt,i,j,year,month,day,hour,mins,secs,ranseed,linelength
      double precision getDouble
      logical getOption
      character*25 tempstr
 
      if (getOption('temperature')) then
        retitemp = getDouble(words(1),25.0d0)
        write(printstring,"(A,F5.1,A)")" Metropolis-test temperature ",retitemp," C"
        call prettyLine(INFO,printstring)
      endif

      if (getOption('iterations')) then
        niter = getInt(words(1),1)
        write(printstring,*)" Iterations = ",niter
        call prettyLine(INFO,printstring) 
      endif
   
      if (getOption('retiblocks')) then
        nretiblock = getInt(words(1),2)
        write(printstring,*)" Reti blocks = ",nretiblock
        call prettyLine(INFO,printstring)
      else
        nretiblock = 2
      endif

      if (getOption('ptblocks')) then
        nptblock = getInt(words(1),2)
        write(printstring,*)" Parallel tempering blocks = ",nptblock
        call prettyLine(INFO,printstring)
      else
        nptblock = 0
      endif

      nlambda = 0
      if (getOption('lambdas')) then
        nlambda = nwords
        do i=1,nlambda
          lambdaval(i) = getDouble(words(i),0.0d0)
        enddo
        write(printstring,*)" Will start simulations at the following lambda-values:"
        call prettyLine(INFO,printstring)
        do i=1,nlambda
          write(printstring,'(F8.3)')lambdaval(i)
          call prettyLine(INFO,printstring)        
        enddo
      endif

      do i=1,nlambda
        call writeNumberString("temperatures",i,tempstr)
        write(tempstr,"(2A)")tempstr(1:linelength(tempstr)),"_"
        if (getOption(tempstr)) then
          ntemp(i) = nwords
          do j=1,ntemp(i)
            pttemperature(i,j) = getDouble(words(j),retitemp)
          enddo
          if (pttemperature(i,1).ne.retitemp) pttemperature(i,1) = retitemp
        else 
          ntemp(i) = 1
          pttemperature(i,1) = retitemp
        endif  
        write(printstring,*)" Temperature(s) for lambda ",i," : ",(int(pttemperature(i,j)),j=1,ntemp(i))
        call prettyLine(INFO,printstring)
        do j=1,ntemp(i)
          pttemperature(i,j) = pttemperature(i,j)+DEG2KEL
          pt_beta_t(i,j) = dble(1.0) / (GASR*pttemperature(i,j))
        enddo
      enddo
      retitemp = retitemp + DEG2KEL
      mt_beta_t = dble(1.0) / (GASR*retitemp)

      if (getOption('slavescript')) then
        slave = words(1)
        write(printstring,*)" Slave script = ",slave
        call prettyLine(INFO,printstring)
      endif

      if (getOption('outdir')) then
        outdir = words(1)
        write(printstring,*)" Output directory = ",outdir
        call prettyLine(INFO,printstring)
      endif  


      
c     set the random number seed based on the current time
      call getDateAndTime(year,month,day,hour,mins,secs)
      ranseed = 1+hour+secs*100+mins*10000+day*1000000+month*100000000
      
c     make into an i6 or i7 positive, odd integer...
      if (ranseed.lt.0) ranseed = -ranseed
      if (mod(ranseed,2).eq.0) ranseed = ranseed - 1
      if (ranseed.lt.100000) ranseed = ranseed+110000

c     initialise the random number generator with this seed
      call sgrnd(ranseed)

      end subroutine
