      
      subroutine proteinPmeEnergy(exclnrg,selfnrg,charges,x,y,z,natoms,maxpmeatoms)
      include 'dimensions.inc'
      include 'simparams.inc'
      include 'constants.inc'
      include 'proteins.inc'
      include 'printer.inc'
      include 'enums.inc'
      include 'parameters.inc'
      include 'flags.inc'
      include 'quickpot.inc'
c#######################################################################
c
c     This subroutine calculates necessary energies for 
c     particle-mesh Ewald for all protein residues
c
c     (C) Samuel Genheden, October 2013
c
c#######################################################################

      integer iprot,ires
      integer natoms,maxpmeatoms
      double precision exclnrg,q2sum,selfnrg,exclpronrg,selfpronrg
      double precision charges(maxpmeatoms),x(maxpmeatoms),y(maxpmeatoms),z(maxpmeatoms)
      integer i,j,k,istart,iend,inats
      include 'nbenergyvar2.inc'
      double precision charge      
c     No softcore between residues
      double precision scl,sclf,sclb
      logical take

      scl = 1.0d0
      sclf = 1.0d0
      sclb = 1.0d0

      exclpronrg = ZERO
      selfpronrg = ZERO
      q2sum = ZERO
      softcore = .false.
      ddd = .false.
c     no periodic boundaries in intramolecular NB interactions
      boundx = .false.
      boundy = .false.
      boundz = .false.
            
      do iprot=1,NProteins
        funciprot1 = iprot
        funciprot2 = iprot
        do ires=1,NResidues(iprot)
          istart = BBZid(iprot,ires,1)
          iend = BBZid(iprot,ires,4)
        
          inats = iend - istart + 1
          do i=1,inats
            funcat1 = istart+i-1
            do j=i+1,inats
              funcat2 = istart+j+1
          
              take = .True.
              do k=1,NResNB(iprot,ires)
                if (((ResNB(iprot,ires,i,1).eq.funcat1).and.(ResNB(iprot,ires,i,2).eq.funcat2)).or.((ResNB(iprot,ires,i,2).eq.funcat1).and.(ResNB(iprot,ires,i,1).eq.funcat2))) take = .false.
              enddo
              do k=1,NResNB14(iprot,ires)
                if (((ResNB14(iprot,ires,i,1).eq.funcat1).and.(ResNB14(iprot,ires,i,2).eq.funcat2)).or.((ResNB(iprot,ires,i,2).eq.funcat1).and.(ResNB(iprot,ires,i,1).eq.funcat2))) take = .false.
              enddo


              if (take) then
c               initialise the calculation
                include 'nbinitial.inc'
c               set up for a solute-solute interaction
                include 'nbpropro.inc'
c               set up the charge
                include 'qinitial.inc'
c               actually perform the calculation
                include 'qenergy.inc'

                exclpronrg = exclpronrg + funccnrg
              endif
c           End of second atom loop
            enddo 

c           Now calculate the sum of the charges squared (for self-energy)
            ipar1 = ZPar(funciprot1,funcat1)
            ipar2 = ZPar(funciprot2,funcat2)
c           set up the charge
            include 'qinitial.inc'
            q2sum = q2sum + funcichg1*funcichg1
       
c           Load arrays for mesh calculation
            natoms = natoms + 1
            x(natoms) = ZCoords(funciprot1,funcat1,1)/10.0d0
            y(natoms) = ZCoords(funciprot1,funcat1,2)/10.0d0
            z(natoms) = ZCoords(funciprot1,funcat1,3)/10.0d0
            charges(natoms) = funcichg1
c         End of first atom loop
          enddo   
                      
        enddo
      enddo

      selfpronrg = -q2sum*FAC4PIEPS0*pme_alpha/sqrt(PI)
      exclpronrg = -exclpronrg

      selfnrg = selfnrg + selfpronrg
      exclnrg = exclnrg + exclpronrg
      return
      end
