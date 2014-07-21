       subroutine savepoint(igcs,movnum)
c      A little routine which is activated if an accepted move with theta
c      greater than 0.95 occurs. It looks through the grid points, finds the 
c      one which is closest, then places a marker on that point
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

       do 70  X = (coordx),(gridx+coordx),space
       do 80  Y = (coordy),(gridy+coordy),space
       do 90  Z = (coordz),(gridz+coordz),space

       dist2 = ( (X-GCSCentre(igcs,1))**2 + (Y-GCSCentre(igcs,2))**2 + (Z-GCSCentre(igcs,3))**2 )
       dist = sqrt(dist2)

       if (dist.lt.0.6) then

       cavCoordDepth = cavCoordDepth + 1
       stackGridX(cavCoordDepth,1) = X 
       stackGridY(cavCoordDepth,2) = Y
       stackGridZ(cavCoordDepth,3) = Z
       endif
      
90     continue
80     continue
70     continue
       return
       end
      



