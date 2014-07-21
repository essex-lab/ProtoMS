      double precision function ureyBradleyEnergy(ityp,ival,iuby,
     .                                     ureyBradleyEnergyB,ureyBradleyEnergyF)
      include 'dimensions.inc'
      include 'printer.inc'
      include 'proteins.inc'
      include 'constants.inc'
      include 'solutes.inc'
      include 'flags.inc'
      include 'enums.inc'
      include 'parameters.inc'
c#######################################################
c
c     This is the function used to evaluate the energy
c     of urey-bradley term iuby in ityp. It returns the 
c     energy in kcal mol-1, or ZERO if the Urey-Bradley term
c     doesn't exist...
c
c     (C) Christopher Woods
c
c#######################################################

      integer iprot,iuby,isol,ityp,iparam,ival,ires1,ires2
      double precision ir,ik,ureyBradleyEnergyB,ureyBradleyEnergyF,iubyval
      double precision iubyvalf,iubyvalb,irf,ikf,irb,ikb
      logical ivariable
      
c     UbyEnergy uses harmonic potential...
c     UbyEnergy = 0.5 * k * (R - R0)**2      
                  
      UreyBradleyEnergy = ZERO
      UreyBradleyEnergyF = ZERO
      UreyBradleyEnergyB = ZERO                  
                       
      ivariable = .false.
      if (ityp.eq.PROTEIN) then
        iprot = ival
        ires1 = Zid2ResID(iprot,UreyBradley(iprot,iuby,1),1)
        ires2 = Zid2ResID(iprot,UreyBradley(iprot,iuby,2),1)
c       don't bother evaluating bonds which are within fixed residues   
        if ( FixedRes(iprot,ires1) .and. FixedRes(iprot,ires2) ) then
           return
        endif
        iparam = UreyBradley(iprot,iuby,3)
        iubyval = UreyBradleySize(iprot,iuby)
        iubyvalf = iubyval
        iubyvalb = iubyval
      else if (ityp.eq.SOLUTE) then
        isol = ival
        iparam = SolUreyBradley(isol,iuby,3)
        iubyval = SolUreyBradleySize(isol,iuby)
        if (NSolVariable(isol).gt.0.and.DifferentLambdas) then
          ivariable = .true.
          iubyvalf = SolUreyBradleySize(NSolForwards(isol),iuby)
          iubyvalb = SolUreyBradleySize(NSolBackwards(isol),iuby)
        else
          iubyvalf = iubyval
          iubyvalb = iubyval
        endif
      else 
        isol = 0
        iprot = 0
        iparam = 0
        iubyval = ZERO
        iubyvalf = iubyval
        iubyvalb = iubyval
        write(printstring,*) "Unknown system type requested ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug")
        return
      endif
      
      if (iparam.lt.0) then
        ik = VariableUByParam(-iparam,2,1)
        ir = VariableUByParam(-iparam,2,2)
      else if (iparam.eq.0) then
        ureyBradleyEnergy = ZERO
        ureyBradleyEnergyF = ZERO
        ureyBradleyEnergyB = ZERO
        return
      else 
        ik = UByParam(iparam,1)
        ir = UByParam(iparam,2)
      endif

      ureyBradleyEnergy = ik * (iubyval-ir)**2

      if (DifferentLambdas) then
        if (iparam.gt.0 .and. .not.ivariable) then
          ureyBradleyEnergyf = ureyBradleyEnergy
          ureyBradleyEnergyb = ureyBradleyEnergy
        else
          if (iparam.lt.0) then
c           backwards
            ikb = VariableUByParam(iparam,4,1)
            irb = VariableUByParam(iparam,4,2)
c           forwards
            ikf = VariableUByParam(iparam,5,1)
            irf = VariableUByParam(iparam,5,2)
          else
            ikf = ik
            ikb = ik
            irf = ir
            irb = ir
          endif
          ureyBradleyEnergyB = ikb * (iubyvalb-irb)**2
          ureyBradleyEnergyF = ikf * (iubyvalf-irf)**2
        endif
      else
        ureyBradleyEnergyF = ureyBradleyEnergy
        ureyBradleyEnergyB = ureyBradleyEnergy
      endif 
      
      return
      
      end
