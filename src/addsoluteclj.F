      subroutine addSoluteCLJ(isol,state,factor)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'constants.inc'
      include 'splitter.inc'
      include 'parameters.inc'
      include 'solutes.inc'
      include 'enums.inc'
c##########################################################
c
c     This subroutine add CLJ parameters add for solute
c     isol based on the already loaded parameters
c     at the moment it will only modify the charges
c
c     (C) Samuel Genheden, May 2013
c
c##########################################################
      integer isol,state
      double precision factor,ichg,isig,ieps,jchg
      double precision charge
      integer i,ipar,jpar,ipar1,ipar2,addCLJ,getCLJ

      do i=1,NSolAtoms(isol)
        ipar = ZSolPar(isol,i)
        if (ipar.lt.1) then
          write(printstring,*)"Cannot deal with solutes that already have variable parameters!"  
          call printline(WARNING,printstring)
          write(printstring,*)" Will not add any additional parameters"  
          call printline(WARNING,printstring)
          return
        endif
      enddo

      do  i=1,NSolAtoms(isol)
        ipar = CLJType(ZSolPar(isol,i))
        if (state.lt.3) then
          jpar = ipar + NSolAtoms(isol) + 1
          do while ((GotCLJ(jpar,1).eqv..true.).and.(jpar.lt.MAXCLJ)) 
            jpar = jpar + NSolAtoms(isol) + 1
          enddo
          if (jpar.ge.MAXCLJ) then
            write(printstring,*)" Could not find a free parameter index",jpar,ipar  
            call printline(FATAL,printstring)
            return          
          endif
          ChgDB(jpar) = ChgDB(ipar)*factor
          AtmDB(jpar) = AtmDB(ipar)
          SigDB(jpar) = SigDB(ipar)
          EpsDB(jpar) = EpsDB(ipar)
          ProtonNumDB(jpar) = ProtonNumDB(ipar)
          GotCLJ(jpar,1) = .true. 
        else
          ChgDb(ipar) = ChgDB(ipar)*factor
        endif

        if (state.eq.1) then
          ipar1 = addCLJ(jpar)
          ipar2 = addCLJ(ipar)
        else if (state.eq.2) then
          ipar1 = addCLJ(ipar)
          ipar2 = addCLJ(jpar)
        endif
              
        if (state.lt.3) then
          ZSolPar(isol,i) = getCLJ(ipar1,ipar2)
          write(printstring,"(A,I2,A,F7.3,F7.3,A)")
     .   "New parameters for atom ",i," charges = (",VariableChg(-ZSolPar(isol,i),1),VariableChg(-ZSolPar(isol,i),3),")"
        else
         write(printstring,"(A,I2,A,F7.3,F7.3,A)")
     .   "New parameters for atom ",i," charges = (",Charge(ZSolPar(isol,i)),Charge(ZSolPar(isol,i)),")"
        endif
        call printline(INFO,printstring)

      enddo
      
      if (state.eq.3) call calculateCLJPairs

      end subroutine
