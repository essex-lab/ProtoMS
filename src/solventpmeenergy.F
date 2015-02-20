      subroutine solventPmeEnergy(exclnrg,selfnrg,charges,x,y,z,natoms,maxpmeatoms,initialize)
      include 'dimensions.inc'
      include 'solvents.inc'
      include 'constants.inc'
      include 'simparams.inc'
      include 'enums.inc'
      include 'printer.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'quickpot.inc'
c#######################################################################
c
c     This subroutine calculates necessary energies for 
c     particle-mesh Ewald for all solvent molecules
c
c     (C) Samuel Genheden, October 2013
c
c#######################################################################

      integer natoms,maxpmeatoms
      double precision exclnrg,selfnrg,q2sum
      double precision charges(maxpmeatoms),x(maxpmeatoms),y(maxpmeatoms),z(maxpmeatoms)
      logical initialize
      integer i,j,k
      include 'nbenergyvar2.inc'
      double precision dx,dy,dz,r,rinv,nrg,Charge,qq,qq1,qq2,getLowScale,getHighScale

      q2sum = ZERO
      do i=1,NSolvents
        do j=1,NSvnAtoms(i)
c         Load arrays for mesh calculation
          natoms = natoms + 1
          x(natoms) = ZSvnCoords(i,j,1)/10.0d0
          y(natoms) = ZSvnCoords(i,j,2)/10.0d0
          z(natoms) = ZSvnCoords(i,j,3)/10.0d0

c         set up the charge
          ipar1 = ZSvnPar(i,j)
          ipar2 = ZSvnPar(i,j)
          include 'qinitial.inc'
          charges(natoms) = funcichg1

c         Only do this once, the waters are fixed so the exclusion and self-energy is constant
          if (initialize) then 
c            if (ipar1.lt.0) then
c              qq1 = VariableChg(-ipar1,1)*VariableChg(-ipar1,1)
c              qq2 = VariableChg(-ipar1,3)*VariableChg(-ipar1,3)
c              q2sum = q2sum + getLowScale(Lambda)*qq1+getHighScale(Lambda)*qq2
c            else
              q2sum = q2sum + funcichg1*funcichg1
c            endif
            
            do k=j+1,NSvnAtoms(i)
              ipar1 = ZSvnPar(i,j)
              ipar2 = ZSvnPar(i,k)
              include 'qinitial.inc'
              dx = ZSvnCoords(i,j,1)-ZSvnCoords(i,k,1)
              dy = ZSvnCoords(i,j,2)-ZSvnCoords(i,k,2)
              dz = ZSvnCoords(i,j,3)-ZSvnCoords(i,k,3)
              r = sqrt(dx**2 + dy**2 + dz**2)
              rinv = 1.0d0/r
              qq = funcichg1*funcichg2
              nrg = qq*rinv*erf(pme_alpha*r)
              exclnrg = exclnrg + nrg
            enddo              
          endif       
        enddo
      enddo

      if (initialize) then
        exclnrg = -exclnrg*FAC4PIEPS0
        selfnrg = -q2sum*FAC4PIEPS0*pme_alpha/dsqrt(PI)
      endif

c      write(6,*)q2sum

      return

      end subroutine
