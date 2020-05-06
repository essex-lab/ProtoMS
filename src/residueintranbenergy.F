      
      double precision function residueIntraNBEnergy(iprot,ires,
     .             cnrg,ljnrg,cnrgb,ljnrgb,cnrgf,ljnrgf)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'enums.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'quickpot.inc'
c#############################################################
c
c     This function calculates the intraresidue non-bonded
c     energy for residue ires in protein iprot
c
c     (C) Christopher Woods   13 March 2002
c
c#############################################################

      integer iprot,ires
      double precision cnrg,ljnrg,cnrgb,cnrgf,ljnrgb,ljnrgf
      integer i
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
c     No softcore between residues
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0

      residueIntraNBEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgb = ZERO
      ljnrgb = ZERO
      cnrgf = ZERO
      ljnrgf = ZERO
      softcore = .false.
      ddd = .false.
c     no periodic boundaries in intramolecular NB interactions
      boundx = .false.
      boundy = .false.
      boundz = .false.
            
      funciprot1 = iprot
      funciprot2 = iprot

      if (QuickMove .and. DDDon) ddd = .true.
c     run over all of the nb pairs in this residue
      if (NResNB(iprot,ires).gt.0) then
        do i=1,NResNB(iprot,ires)
          funcat1 = ResNB(iprot,ires,i,1)
          funcat2 = ResNB(iprot,ires,i,2)
          
          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'
            
          cnrg = cnrg + funccnrg
          ljnrg = ljnrg + funcljnrg
          cnrgf = cnrgf + funccnrgf
          cnrgb = cnrgb + funccnrgb
          ljnrgf = ljnrgf + funcljnrgf
          ljnrgb = ljnrgb + funcljnrgb
        enddo
      endif
      
c     run over all of the 1-4 pairs in this residue
      if (NResNB14(iprot,ires).gt.0) then
        do i=1,NResNB14(iprot,ires)
          funcat1 = ResNB14(iprot,ires,i,1)
          funcat2 = ResNB14(iprot,ires,i,2)
          
          include 'nbinitial.inc'
          include 'nbpropro.inc'
          include 'nbenergy2.inc'
            
          cnrg = cnrg + OneOverScl14C*funccnrg
          ljnrg = ljnrg + OneOverScl14LJ*funcljnrg
          cnrgf = cnrgf + OneOverScl14C*funccnrgf
          cnrgb = cnrgb + OneOverScl14C*funccnrgb
          ljnrgf = ljnrgf + OneOverScl14LJ*funcljnrgf
          ljnrgb = ljnrgb + OneOverScl14LJ*funcljnrgb
        enddo
      endif
              
      residueIntraNBEnergy = cnrg + ljnrg
      return
      end
