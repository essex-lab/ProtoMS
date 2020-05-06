      subroutine scanForVariableGroups
      include 'dimensions.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'flags.inc'
      include 'printer.inc'
c##########################################################
c
c     This routine finds out whether or not we
c     will be need to do single-copy free energy - the
c     routine doesn't yet scan for variable intramolecular
c     parameters!
c
c     (C) Christopher Woods   
c
c##########################################################

c     If lambdas are the same then definitely false as no
c     need for free energy simulation!
      if (dabs(Lambda-LambdaF).lt.TINY .and. 
     .    dabs(Lambda-LambdaB).lt.TINY) then
c       lambdas are the same, so no need for SingleCopy
        write(printstring,10) Lambda
10      format('Performing simulation at lambda = ',f8.3)
        call prettyLine(INFO,printstring)
        DifferentLambdas = .false.
      else
        DifferentLambdas = .true.
        write(printstring,20) LambdaB,Lambda,LambdaF
20      format('Running simulation in lambda window ',f8.3,
     .         ' <-- ',f8.3,' --> ',f8.3)
        call prettyLine(INFO,printstring)
      endif
      return
      end
