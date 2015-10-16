       subroutine savelow(igcs,movnum)
c      A little routine which is activated if an accepted move with theta
c      less than 0.05 occurs. it prints out 
       include 'dimensions.inc'
       include 'metropolis.inc'
       include 'enums.inc'
       include 'simparams.inc'
       include 'constants.inc'
       include 'flags.inc'
       include 'proteins.inc'
       include 'solvents.inc'
       include 'solutes.inc'
       include 'printer.inc'
       include 'parameters.inc'
       include 'quickpot.inc'
       include 'moves.inc'
       include 'gb.inc'
       include 'energies.inc'
       include 'stacks.inc'
       
       double precision X,Y,Z,dist,dist2
       integer igcs,movnum

       X = GCSCentre(igcs,1)
       Y = GCSCentre(igcs,2)
       Z = GCSCentre(igcs,3)

       lowstack = lowstack + 1
       stackGridXL(lowstack,1) = X 
       stackGridYL(lowstack,2) = Y
       stackGridZL(lowstack,3) = Z

       return
       end
      



