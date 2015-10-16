       subroutine savehigh(igcs,movnum)
c      A little routine which is activated if an accepted move with theta
c      greater than 0.95 occurs. it prints out 
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

       highstack = highstack + 1
       stackGridX(highstack,1) = X 
       stackGridY(highstack,2) = Y
       stackGridZ(highstack,3) = Z

       return
       end
      



