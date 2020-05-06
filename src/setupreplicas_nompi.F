      subroutine setupReplicas_nompi
      include 'dimensions.inc'
      include 'flags.inc'
      include 'printer.inc'

      logical getOption

      doLambdaRE = .false.
      doTemperatureRE = .false.
      doMultiGCMC = .false.
      doMultiJaws2 = .false.

      if (getOption('temperaturere').or.getOption('lambdare').or.
     .     getOption('multigcmc').or.getOption('multijaws2')) then
          call printLine(FATAL,"ProtoMS was compiled without MPI. "//
     .                         "Replica exchange is not available.")
          call closeProgram(-1, "Replica exchange not available.")
      endif

      end
