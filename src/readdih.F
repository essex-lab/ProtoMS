      subroutine readDih(FILE,line,changemode,endfile,ndih)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'enums.inc'
c##############################################################
c
c     This routine reads in the dihedral parameters. These
c     are quite different to the bond/angle/urey 
c     parameters as *multiple* dihedral parameters can
c     be assigned to each dihedral! (this is the CHARMM model).
c
c     Each dihedral term has the following functional
c     form;
c
c     E_dih(phi) = k1 * ( 1 + k2*( cos(k3*phi - k4) ))
c
c     k1 should have units of kcal mol-1
c     k2 is typically equal to +1 or -1 and is unitless
c     k3 is typically a positive integer and is unitless
c     k4 should be given in units of *degrees*
c
c     Each dihedral can have multiple dihedral terms assigned to
c     it. This is in contrast to OPLS or AMBER that have a
c     single term for each dihedral, but this term is composed
c     of multiple cosine terms, e.g. for OPLS
c
c     E_dih(phi) = v0 + 0.5*( v1*(1.0 + cos(phi+f1) + 
c                             v2*(1.0 - cos(2*phi+f2)) + 
c                             v3*(1.0 + cos(3*phi+f3)))
c
c     This dihedral can be represented via four dihedral terms;
c     (1) k1 = v0, k2 = 0, k3 = 0, k4 = 0
c     (2) k1 = 0.5*v1, k2 = +1, k3 = 1.0, k4 = f1
c     (3) k1 = 0.5*v2, k2 = -1, k3 = 2.0, k4 = f2
c     (4) k1 = 0.5*v3, k2 = +1, k3 = 3.0, k4 = f3
c
c     Note that in OPLS k4 is given in units of PI/2 radians (e.g. 1.0 = 90 degrees)
c
c     While this scheme is more complicated than OPLS or AMBER, this
c     CHARMM scheme is much more flexible, allowing for many dihedral
c     terms to be applied, thus allowing complete freedom for the 
c     parameters to fit the desired potential surface.
c
c     Each dihedral term is indexed by a parameter number. A dihedral
c     is specified by a set of four AMBER atom types, e.g.
c
c     par 1 a1 a2 a3 a4
c     par 2 b1 b2 b3 b4
c     par 3 c1 c2 c3 c4
c     par 4 d1 d2 d3 d4
c
c     atm HC CT CT HC 1 2 3 4
c
c     says that the dihedral between atoms with AMBER atom types
c     HC-CT-CT-HC has a potential made of four dihedral terms, formed
c     with parameters (a1,a2,a3,a4), (b1,b2,b3,b4), (c1,c2,c3,c4) and
c     (d1,d2,d3,d4)
c
c##################################################################

      integer FILE
      character*(*) line
      logical changemode,endfile
      integer ndih,id
      integer findDihedralIndex,readInt
      double precision readDouble
      character*2 iatm1,iatm2,iatm3,iatm4
      integer ipar,npar,i
      double precision ik1,ik2,ik3,ik4
      logical ok
      integer idih(MAXDIHTERMSPERDIHEDRAL)
      
      changemode = .false.
      endfile = .false.
      ndih = 0
      
100   continue
        read(FILE,10,err=9700,end=9800) line
10      format(a300)

c       do not lowercase as AMBER types are case-sensitive
        
        if (line(1:1).eq.'#') goto 100
        
        call split(line)
        if (nwords.lt.2) goto 100
        
        call lowercase(words(1))
        
        if (words(1).eq.'mode') then
          changemode = .true.
          call lowercase(words(2))
          return
        else if (words(1).eq.'term') then
c         read in the dihedral cosine term, the format is;
c         "term ipar ik1 ik2 ik3 ik4"
          if (nwords.lt.6) then
            call printLine(WARNING,"Not enough parameters for dihedral parameter line")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the parameters
          ipar = readInt(words(2))
          ik1 = readDouble(words(3))
          ik2 = readDouble(words(4))
          ik3 = readDouble(words(5))
c         read in ik4 and convert from degrees to radians          
          ik4 = readDouble(words(6)) * PI/ONEEIGHTY
          
          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative dihedral terms!")
            call printLine(WARNING,line)
            goto 100
          else if (ipar.gt.MAXDIHTERMS) then
            write(printstring,*) "Parameter ",ipar," is too large (",MAXDIHTERMS,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
            call printLine(WARNING,"Increase MAXDIHTERMS and recompile")
            goto 100
          endif
          
          if (GotDihTerm(ipar)) then
            write(printstring,20) ipar,DihTerm(ipar,1),DihTerm(ipar,2),ik1,ik2
20          format("Overwriting dihedral term ",i4,"(",f8.3,",",f8.3,"...) with (",
     .                 f8.3,",",f8.3,"...)")
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif
          
          DihTerm(ipar,1) = ik1
          DihTerm(ipar,2) = ik2
          DihTerm(ipar,3) = ik3
          DihTerm(ipar,4) = ik4
          
          GotDihTerm(ipar) = .true.
          
          ndih = ndih + 1
        else if (words(1).eq.'par') then
c         this is a dihedral parameter - this is formed as a 
c         combination of dihedral terms, format is;
c         par ipar iterm1 iterm2 iterm3 ... up to MAXDIHTERMS  #comment
c         Each dihedral parameter must have at least one term
          
          if (nwords.lt.3) then
            call printLine(WARNING,"Not enough terms in this dihedral parameter!")
            call printLine(WARNING,line)        
            goto 100
          endif

c         read in the parameters - work out the number of parameters - this is the rest of the
c         contents of the line, up to any comment (indicated by a '#')
          npar = 0
          do i=3,nwords
            if (words(i)(1:1).eq.'#') then
c             from here on in we have a comment! - break out of the loop
              goto 200
            endif
            npar = npar + 1
          enddo
200       continue
            
          if (npar.eq.0) then
            call printLine(WARNING,"You must specify some terms for this dihedral parameter!")
            call printLine(WARNING,line)
            goto 100
          else if (npar.gt.MAXDIHTERMSPERDIHEDRAL) then
            write(printstring,*) "Exceeded maximum number of terms per dihedral (",
     .                           MAXDIHTERMSPERDIHEDRAL,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
            call printLine(WARNING,"Increase MAXDIHTERMSPERDIHEDRAL and recompile")
            goto 100
          endif
               
c         read in the dihedral parameter
          ipar = readInt(words(2))
          if (ipar.le.0) then
            call printLine(WARNING,"Cannot have zero or negative dihedral parameter!")
            call printLine(WARNING,line)
            return
          else if (ipar.gt.MAXDIHPARAM) then
            write(printstring,*) "Exceeded maximum number of dihedral parameters (",
     .                           MAXDIHPARAM,")"
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
            call printLine(WARNING,"Increase MAXDIHPARAM and recompile")
            goto 100
          else if (GotDihParam(ipar)) then
            write(printstring,*) "Overwriting dihedral parameter ",ipar
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          endif        

          GotDihParam(ipar) = .true.
                      
          ok = .true.
          do i=1,npar
            idih(i) = readInt(words(i+2))
c           ensure that all of these terms are within bounds!
            if (idih(i).lt.0) then
              write(printstring,*) "Cannot refer to a negative dihedral term! ",idih(i)
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
              ok = .false.
            else if (idih(i).gt.MAXDIHTERMS) then
              write(printstring,*) "You cannot refer to a dihedral term above MAXDIHTERMS, ",MAXDIHTERMS
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
              ok = .false.
            else if (.not.GotDihTerm(idih(i))) then
              write(printstring,*) "Referring to a dihedral term that has not yet been loaded! ",idih(i)
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
              ok = .false.
            endif
          enddo
          
          if (.not. ok) then
            call printLine(WARNING,
     .       "There were errors with this dihedral parameter so it will be set as a dummy!")
            NDihTerms(ipar) = 0
          else
            if (NDihTerms(ipar).gt.0) then
              write(printstring,*) "Overwriting dihedral parameter ",ipar
              call printLine(WARNING,printstring)
              call printLine(WARNING,line)
            endif
                    
c           save the terms for this dihedral
            NDihTerms(ipar) = npar
            do i=1,npar
              DihTermList(ipar,i) = idih(i)
            enddo
          endif
        else if (words(1).eq.'atm') then
c         read in amber atoms and associated parameter. Format is;
c         atm iatm1 iatm2 iatm3 iatm4 ipar
          
          if (nwords.lt.6) then
            call printLine(WARNING,"Incorrect format of dihedral atom line!")
            call printLine(WARNING,line)
            goto 100
          endif
          
c         read in the atoms
          iatm1 = words(2)(1:2)
          iatm2 = words(3)(1:2)
          iatm3 = words(4)(1:2)
          iatm4 = words(5)(1:2)
          
          ipar = readInt(words(6))
          
          call lowercase(words(2))
          call lowercase(words(3))
          call lowercase(words(4))
          call lowercase(words(5))
          
c         dihedral with dummy atoms must be zero
          if (words(2)(1:2).eq.'dm' .or. words(3)(1:2).eq.'dm' .or. 
     .        words(4)(1:2).eq.'dm' .or. words(5)(1:2).eq.'dm') then
            if (ipar.ne.0) then
              call printLine(WARNING,"Cannot have non-zero dihedral parameter for dihedral containing dummies!")
              call printLine(WARNING,line)
            endif
            goto 100
          endif
          
          call sortDihedral(iatm1,iatm2,iatm3,iatm4)
          
          id = findDihedralIndex(iatm1,iatm2,iatm3,iatm4)
          
          if (id.ne.0) then
c           we are overwriting an existing parameter
            write(printstring,30) iatm1,iatm2,iatm3,iatm4
30          format("Overwriting dihedral parameter for ",a2,"-",a2,"-",a2,"-",a2)
            call printLine(WARNING,printstring)
            call printLine(WARNING,line)
          else
c           load the new dihedral term into the array
            id = NDihParam + 1
            if (id.gt.MAXDIHPARAM) then
              write(printstring,*) "Exceeded maximum number of dihedral parameters (",
     .                             MAXDIHPARAM,")"
              call printLine(WARNING,printstring)
              call printLine(WARNING,"Increase MAXDIHPARAM and recompile")
              goto 100
            endif
              
            DihName(id,1) = iatm1
            DihName(id,2) = iatm2
            DihName(id,3) = iatm3
            DihName(id,4) = iatm4
            NDihParam = id            
          endif

          if (streamHandle(DETAIL).ge.0) then
            write(printstring,40) iatm1,iatm2,iatm3,iatm4,ipar
40          format("Setting dihedral parameter between AMBER atoms ",a,"-",a,"-",a,"-",a," to ",i5)
            call printLine(DETAIL,printstring)
          endif          

          DihParamIndex(id) = ipar
          
          ndih = ndih + 1
        else
          call printLine(WARNING,"Unrecognised dihedral line!")
          call printLine(WARNING,line)
        endif
c     end of DO-WHILE loop
      goto 100
9700  continue
      call printLine(WARNING,"Error reading parameter file (DIHEDRAL)!")
9800  continue
      endfile = .true.
      return
      
      end
