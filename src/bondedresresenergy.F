      
      double precision function bondedResResEnergy(iprot,sres1,sres2,cnrg,ljnrg, 
     .                                    cnrgb,ljnrgb,cnrgf,ljnrgf,mindist)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'parameters.inc'
      include 'proteins.inc'
      include 'quickpot.inc'
c#######################################################
c
c     This is an internal function only intended to be
c     used by residueResidueEnergy, used to calculate
c     the residue-residue energy between residues within
c     a protein that are only bonded via the backbone...
c
c     (C) Christopher Woods   25 January 2002
c
c########################################################

      integer iprot,ires1,ires2
      integer sres1,sres2,ibndindex,nats1,nats2
      double precision cnrg,ljnrg,cnrgf,cnrgb,ljnrgf,ljnrgb,mindist
      integer i,j,istrt1,iend1,istrt2,iend2
      integer getResResBondIndex
      include 'nbenergyvar2.inc'
      double precision charge,sigma,epslon      
      logical bondmap(MAXSCATOMS+5,MAXSCATOMS+5)
c     JM : The scl values are 1 because there cannot be any soft core
c     between two protein side chains with the current implementation
      double precision scl,sclf,sclb
      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0
      bondedResResEnergy = ZERO
      cnrg = ZERO
      ljnrg = ZERO
      cnrgf = ZERO
      cnrgb = ZERO
      ljnrgf = ZERO
      ljnrgb = ZERO
      softcore = .false.
      ddd = .false.
c     numerically sort the residues
      if (sres1.gt.sres2) then
        ires1 = sres2
        ires2 = sres1
      else
        ires1 = sres1
        ires2 = sres2
      endif

c     get the index of the residue pair
      ibndindex = getResResBondIndex(iprot,ires1,ires2)
      if (ibndindex.le.0) then
        write(printstring,*) "Cannot find the index for bonded residues ",ires1,ires2
        call printLine(WARNING,printstring)
        return
      endif
        
      funciprot1 = iprot
      funciprot2 = iprot

      boundx = .false.
      boundy = .false.
      boundz = .false.
                  
c     loop over all non-bonded pairs between these residues
      if (NResResNB(iprot,ibndindex).gt.0) then
        
c       the array 'ResResNB' contains the *bonded* pairs (to save memory). Invert this
c       to the get the list of nonbonded pairs. The bonded 1-4 pairs are also listed
c       in this array
        istrt1 = BBZid(iprot,ires1,1)
        istrt2 = BBZid(iprot,ires2,1)
        iend1 = BBZid(iprot,ires1,4)
        iend2 = BBZid(iprot,ires2,4)
        
        nats1 = iend1 - istrt1 + 1
        nats2 = iend2 - istrt2 + 1

c       initialise all pairs of atoms as not-bonded              
        do i=1,nats1
          do j=1,nats2
            bondmap(i,j) = .true.
          enddo
        enddo

c       remove all bonded pairs              
        do i=1,NResResNB(iprot,ibndindex)
          bondmap(ResResNB(iprot,ibndindex,i,1)-istrt1+1,
     .            ResResNB(iprot,ibndindex,i,2)-istrt2+1) = .false.
        enddo

c       now loop over all pairs of atoms, and do calculation for nonbonded pairs
        if (QuickMove .and. DDDon) ddd = .true.

        if (ddd) then
c       use unoptimised nbenergy routine to handle implicit solvent
           do i=istrt1,iend1
              do j=istrt2,iend2
                 if (bondmap(i-istrt1+1,j-istrt2+1)) then
                    funcat1 = i
                    funcat2 = j
                    
                    include 'nbinitial.inc'
                    include 'nbpropro.inc'
                    include 'nbenergy2.inc'
                    
                    if (funcmindist.lt.mindist) mindist = funcmindist
                    
                    cnrg = cnrg + funccnrg
                    ljnrg = ljnrg + funcljnrg
                    cnrgf = cnrgf + funccnrgf
                    cnrgb = cnrgb + funccnrgb 
                    ljnrgf = ljnrgf + funcljnrgf
                    ljnrgb = ljnrgb + funcljnrgb
                 endif
              enddo
           enddo
        else
           do i=istrt1,iend1
              do j=istrt2,iend2
                 if (bondmap(i-istrt1+1,j-istrt2+1)) then
                    funcat1 = i
                    funcat2 = j
                    
                    include 'nbinitial.inc'
                    include 'nbpropro.inc'
                    include 'nbenergy_opt.inc'
                    
                    if (funcmindist.lt.mindist) mindist = funcmindist
                    
                    cnrg = cnrg + funccnrg
                    ljnrg = ljnrg + funcljnrg
                    cnrgf = cnrgf + funccnrgf
                    cnrgb = cnrgb + funccnrgb 
                    ljnrgf = ljnrgf + funcljnrgf
                    ljnrgb = ljnrgb + funcljnrgb
                 endif
              enddo
           enddo
           
        endif
        
      endif
      
c     loop over all of the 14-pairs between these residues
      if (ddd) then
        if (NResResNB14(iprot,ibndindex).gt.0) then
          do i=1,NResResNB14(iprot,ibndindex)
            funcat1 = ResResNB14(iprot,ibndindex,i,1)
            funcat2 = ResResNB14(iprot,ibndindex,i,2)
          
            include 'nbinitial.inc'
            include 'nbpropro.inc'
            include 'nbenergy_opt.inc'
              
            if (funcmindist.lt.mindist) mindist = funcmindist
          
            cnrg = cnrg + OneOverScl14C*funccnrg
            ljnrg = ljnrg + OneOverScl14LJ*funcljnrg
            cnrgf = cnrgf + OneOverScl14C*funccnrgf
            cnrgb = cnrgb + OneOverScl14C*funccnrgb 
            ljnrgf = ljnrgf + OneOverScl14LJ*funcljnrgf
            ljnrgb = ljnrgb + OneOverScl14LJ*funcljnrgb
          enddo
        endif
      else
        if (NResResNB14(iprot,ibndindex).gt.0) then
          do i=1,NResResNB14(iprot,ibndindex)
            funcat1 = ResResNB14(iprot,ibndindex,i,1)
            funcat2 = ResResNB14(iprot,ibndindex,i,2)
            
            include 'nbinitial.inc'
            include 'nbpropro.inc'
            include 'nbenergy_opt.inc'
              
            if (funcmindist.lt.mindist) mindist = funcmindist
          
            cnrg = cnrg + OneOverScl14C*funccnrg
            ljnrg = ljnrg + OneOverScl14LJ*funcljnrg
            cnrgf = cnrgf + OneOverScl14C*funccnrgf
            cnrgb = cnrgb + OneOverScl14C*funccnrgb 
            ljnrgf = ljnrgf + OneOverScl14LJ*funcljnrgf
            ljnrgb = ljnrgb + OneOverScl14LJ*funcljnrgb
          enddo
        endif
      endif
      
      bondedResResEnergy = cnrg + ljnrg
      return
      end
