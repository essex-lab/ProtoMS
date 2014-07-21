      subroutine setBackwardsStructure(isol,newsol,tmol)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'solutes.inc'
c#########################################################
c
c     This sets this solute as a backwards free energy
c     structure, that is not used in a normal simulation
c
c     (C) Christopher Woods     14 August 2002
c
c##########################################################

      integer isol,newsol,tmol,i
      double precision perturbForceField
                                   
      NSolBackwards(isol) = newsol
      IgnoreSolute(newsol) = .true.
      
c     assign all of the geometry values to the solute, and rebuild from
c     the new zmatrix
      do i=1,NSolVariable(isol)
        ZSolVal(newsol,SolVariableAtm(isol,i),SolVariableTyp(isol,i)) = 
     .     perturbForceField(SolVariableVal(isol,i,1),SolVariableVal(isol,i,2),lambdaB)
      enddo

      call buildSolute(newsol)
      
      return
      end
