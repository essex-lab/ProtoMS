      subroutine recalculateParameters
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'parameters.inc'
c##############################################################
c
c     This routine recalculates all of the forcefield parameters
c     for the current value of lambda, lambdaf and lambdab
c
c     (C) Christopher Woods
c
c##############################################################

      double precision perturbForcefield
      integer i,j

c#################################################
c#################################################

      if (NVariableCLJ.le.0) return
      do i=1,NVariableCLJ
        VariableChg(i,2) = perturbForceField(VariableChg(i,1),VariableChg(i,3),Lambda)
        VariableSig(i,2) = perturbForceField(VariableSig(i,1),VariableSig(i,3),Lambda)
        VariableEps(i,2) = perturbForceField(VariableEps(i,1),VariableEps(i,3),Lambda)
c       now the forwards and backwards parameters
        VariableChg(i,4) = perturbForceField(VariableChg(i,1),VariableChg(i,3),LambdaB)
        VariableSig(i,4) = perturbForceField(VariableSig(i,1),VariableSig(i,3),LambdaB)
        VariableEps(i,4) = perturbForceField(VariableEps(i,1),VariableEps(i,3),LambdaB)
        VariableChg(i,5) = perturbForceField(VariableChg(i,1),VariableChg(i,3),LambdaF)
        VariableSig(i,5) = perturbForceField(VariableSig(i,1),VariableSig(i,3),LambdaF)
        VariableEps(i,5) = perturbForceField(VariableEps(i,1),VariableEps(i,3),LambdaF)
      enddo

      if (NVarBondParam.gt.0) then
        do i=1,NVarBondParam
          do j=1,2      
            VariableBndParam(i,2,j) = perturbForceField(
     .                               VariableBndParam(i,1,j),
     .                               VariableBndParam(i,3,j),
     .                                                    Lambda)
c           backwards and forwards...
            VariableBndParam(i,4,j) = perturbForceField(
     .                               VariableBndParam(i,1,j),
     .                               VariableBndParam(i,3,j),
     .                                                    LambdaB)
            VariableBndParam(i,5,j) = perturbForceField(
     .                               VariableBndParam(i,1,j),
     .                               VariableBndParam(i,3,j),
     .                                                    LambdaF)
          enddo
        enddo
      endif
      
      if (NVarAngleParam.gt.0) then
        do i=1,NVarAngleParam
          do j=1,2
            VariableAngParam(i,2,j) = perturbForceField(
     .                                   VariableAngParam(i,1,j),
     .                                   VariableAngParam(i,3,j),
     .                                           Lambda)
c           backwards and forwards...
            VariableAngParam(i,4,j) = perturbForceField(
     .                                   VariableAngParam(i,1,j),
     .                                   VariableAngParam(i,3,j),
     .                                           LambdaB)
            VariableAngParam(i,5,j) = perturbForceField(
     .                                   VariableAngParam(i,1,j),
     .                                   VariableAngParam(i,3,j),
     .                                           LambdaF)
          enddo
        enddo
      endif

c     because dihedral parameters are more complex they are dynamically calculated
c     in dihedralEnergy
      
      if (NVarUreyBradleyParam.gt.0) then
        do i=1,NVarUreyBradleyParam
          do j=1,2      
            VariableUByParam(i,2,j) = perturbForceField(
     .                                VariableUByParam(i,1,j),
     .                                VariableUByParam(i,3,j),
     .                                                    Lambda)
c           backwards and forwards...
            VariableUByParam(i,4,j) = perturbForceField(
     .                                VariableUByParam(i,1,j),
     .                                VariableUByParam(i,3,j),
     .                                                    LambdaB)
            VariableUByParam(i,5,j) = perturbForceField(
     .                                VariableUByParam(i,1,j),
     .                                VariableUByParam(i,3,j),
     .                                                    LambdaF)
          enddo
        enddo
      endif
      
      return
      end
