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
       
       double precision X,Y,Z,dist,dist2x,dist2y,dist2z
       integer igcs,movnum,nx,ny,nz,ix,iy,iz

       nx = ceiling(gridx / space)
       ny = ceiling(gridy / space)
       nz = ceiling(gridz / space)

       do 70  ix = 0, nx-1
         X = coordx + ix * space
         dist2x = (X - GCSCentre(igcs, 1))**2

         do 80  iy = 0, ny-1
           Y = coordy + iy * space
           dist2y = (Y - GCSCentre(igcs, 2))**2

           do 90  iz = 0, nz-1
             Z = coordz + iz * space
             dist2z = (Z - GCSCentre(igcs, 3))**2

             dist = sqrt(dist2x + dist2y + dist2z)
             if (dist.lt.0.6) then
               cavCoordDepth = cavCoordDepth + 1
               stackGridX(cavCoordDepth,1) = X 
               stackGridY(cavCoordDepth,2) = Y
               stackGridZ(cavCoordDepth,3) = Z
             endif

90         continue
80       continue
70     continue
       return
       end
      



