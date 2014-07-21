      subroutine setDefaultOptions
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'solvents.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'pressure.inc'
c############################################
c
c     This subroutine sets the initial, 
c     default values of all of the options
c
c     (C) Christopher Woods
c
c############################################

      protodebug = .false.
      testenergy = .false.
      prettyprint = .true.
      lambda = ZERO
      lambdaf = ZERO
      lambdab = ZERO
      temperatureC = 25.0
      cut = 15.0
      feather = 0.5
      cuttype = RESIDUEMOLECULE
      prefSample = .false.
      svnPrefSol = 0
      thres = 0.95
      pressure = ZERO
      maxVolChange = ZERO
      coordx = 0
      coordy = 0
      coordz = 0
      gridx = 0
      gridy = 0
      gridz = 0      


      boundary = AGREESOLVENT
      solventBox = .false.
      solventCap = .false.
                             
      return
      end
