      subroutine biaspotential(igcs,potential)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'stacks.inc'
      include 'solvents.inc'
      include 'constants.inc'

c###################################################################
c     A routine which returns the value of the bias potential for
c     theta sampling moves
c     MSB 02/11/09
c##################################################################

      integer igcs
      double precision Vo,Vc,dGcon,potential
c     Form of the bias potential....
c     E(bias) = [-dG(hyd) + dG(constraint)] theta(igcs)
c     dG(constraint) = -kT.ln(V(const)/V(o))
      Vo = ((1 / AVAGADRO) / 55.56) * 1d27
c     print*, "Vo" ,Vo
c     Value of the constraint - generally a 27A cube...
      Vc = gridx*gridy*gridz

      dGcon = GASR * temperature * log(Vc/Vo)
      potential = (bias - dGcon)*thetaparam(igcs)
      return
      end



