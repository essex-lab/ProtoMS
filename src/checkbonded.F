      
      integer function checkBonded(iprot,sat1,sat2)
      include 'dimensions.inc'
      include 'constants.inc'
      include 'proteins.inc'
c############################################################
c
c     This function checks to see if the atoms 
c     sat1 and sat2 in protein iprot are bonded
c
c     This function also ensures that all of the implicit
c     bonds/angles/dihedrals between atoms in the protein
c     are added if they haven't been already
c
c     Returns 0 if they are not bonded,
c             1 if they are bonded
c             2 if they are angled
c             3 if they are dihedraled
c
c     (C) Christopher Woods, October 2004
c
c############################################################

      integer iprot,ires1,iat1,ires2,iat2,sres1,sres2,sat1,sat2
      integer istrt1,iend1,istrt2,iend2,i,j,k
      integer iu,iv,im,in,ip

      if (sat1.eq.sat2) then
        checkBonded = 0
        return
      endif
      
      sres1 = Zid2ResID(iprot,sat1,1)
      sres2 = Zid2ResID(iprot,sat2,1)

c     make sure the residue are numerically ordered
      if (sres1.gt.sres2) then
        ires1 = sres2
        ires2 = sres1
        iat1 = sat2
        iat2 = sat1
      else
        ires1 = sres1
        ires2 = sres2
        iat1 = sat1
        iat2 = sat2
      endif

      checkBonded = 0        
c     if the residues are not bonded, then the atoms can't be bonded
      if (ires1.ne.ires2) then
        if (.not.ResResBond(iprot,ires1,ires2)) return
      endif
        
      iu = Zid2ResID(iprot,iat1,2) + 5
      iv = Zid2ResID(iprot,iat2,2) + 5

      if (ires1.eq.ires2) then
c       the atoms are in the same residue - see if the connectivity map says they are bonded
        if (TempRes1Map(iu,iv)) then
c         they are!
          checkBonded = 1
          return
        endif
        
c       how about angled?
        istrt1 = BBZid(iprot,ires1,1)
        iend1 = BBZid(iprot,ires1,4)
        do i=istrt1,iend1
          if (i.ne.iat1 .and. i.ne.iat2) then
            in = Zid2ResID(iprot,i,2) + 5
            if (TempRes1Map(iu,in)) then
              if (TempRes1Map(in,iv)) then
c               iu is bonded to in and in to iv so iu and iv are angled!
                checkBonded = 2
c               make sure this angle has been added to the protein zmatrix
                call addProteinAngle(iprot,iat1,i,iat2,.false.,.false.,ZERO,-1,-1)           
                return
              endif
            endif
          endif
        enddo
        
        do i=istrt1,iend1
          if (i.ne.iat1 .and. i.ne.iat2) then
            in = Zid2ResID(iprot,i,2) + 5
            if (TempRes1Map(iu,in)) then
c             how about dihedraled?                
              do j=istrt1,iend1
                if (j.ne.i .and. j.ne.iat1 .and. j.ne.iat2) then
                  ip = Zid2ResID(iprot,j,2) + 5
                  if (TempRes1Map(in,ip) .and. TempRes1Map(ip,iv)) then
c                   iu is bonded to in, in to ip and ip to iv, so 
c                   iu and iv are dihedraled          
                    CheckBonded = 3
c                   make sure that this dihedral has been added to the protein zmatrix
                    call addProteinDihedral(iprot,iat1,i,j,iat2,.false.,.false.,ZERO,-1,-1)
                    return
                  endif
                endif
              enddo
            endif
          endif
        enddo
      else
c       the atoms are in different residues
c       loop over all of the interresidue atoms and see if they are bonded

c       this code will not work if a group of atoms forms a ring between residues, 
c       as there is the chance that a pair of atoms that are angled may be found
c       to be dihedraled first...

        istrt1 = BBZid(iprot,ires1,1)
        iend1 = BBZid(iprot,ires1,4)
        istrt2 = BBZid(iprot,ires2,1)
        iend2 = BBZid(iprot,ires2,4)
        do i=istrt1,iend1
          im = Zid2ResID(iprot,i,2) + 5
          do j=istrt2,iend2
            in = Zid2ResID(iprot,j,2) + 5
c           this requires ires1 and ires2 to be numerically ordered
            if (TempResResMap(im,in)) then
c             the residues are bonded through i and j
              if (i.eq.iat1 .and. j.eq.iat2) then
                checkBonded = 1
                return
              else if (i.eq.iat1) then
c               if iat2 is bonded to j, then this is an angle
                if (TempRes2Map(in,iv)) then
                  CheckBonded = 2
                  call addProteinAngle(iprot,iat1,j,iat2,.false.,.false.,ZERO,-1,-1)
                  return
                endif
c               now get all atoms bonded to j
                do k=istrt2,iend2
                  if (k.ne.j) then
                    ip = Zid2ResID(iprot,k,2) + 5
                    if (TempRes2Map(in,ip)) then
                      if (TempRes2Map(ip,iv)) then
                        CheckBonded = 3
                        call addProteinDihedral(iprot,iat1,j,k,iat2,.false.,.false.,ZERO,-1,-1)
                        return
                      endif
                    endif
                  endif
                enddo
              else if (j.eq.iat2) then
c               if iat1 is bonded to i then this is an angle
                if (TempRes1Map(im,iu)) then
                  CheckBonded = 2
                  call addProteinAngle(iprot,iat1,i,iat2,.false.,.false.,ZERO,-1,-1)
                  return
                endif
C               now get all atoms bonded to i
                do k=istrt1,iend1
                  if (k.ne.i) then
                    ip = Zid2ResID(iprot,k,2) + 5
                    if (TempRes1Map(im,ip)) then
                      if (TempRes1Map(ip,iu)) then
                        CheckBonded = 3
                        call addProteinDihedral(iprot,iat1,k,i,iat2,.false.,.false.,ZERO,-1,-1)
                        return
                      endif
                    endif
                  endif
                enddo
              else
c               to be valid, i must be bonded to iat1, and j must be
c               bonded to iat2
                if (TempRes1Map(im,iu).and.TempRes2Map(in,iv)) then
                  CheckBonded = 3
                  call addProteinDihedral(iprot,iat1,i,j,iat2,.false.,.false.,ZERO,-1,-1)
                  return
                endif
              endif
            endif
          enddo
        enddo
      endif
      
      checkBonded = 0
      return
      end
