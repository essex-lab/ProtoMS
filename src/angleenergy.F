      double precision function angleEnergy(ityp,ival,iang,angleEnergyB,angleEnergyF)
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
c     of angle iang in ityp ival. It returns the 
c     energy in kcal mol-1, or ZERO if the ang
c     doesn't exist...
c
c     (C) Christopher Woods   23 January 2002
c
c#######################################################

      integer iprot,iang,ival,isol,iparam,ityp,ires1,ires2,ires3
      double precision ir,ik,angleEnergyB,angleEnergyF,iangval   
      double precision iangvalf,iangvalb,irf,ikf,irb,ikb
      logical ivariable

c     AngEnergy uses harmonic potential...
c     AngEnergy = 0.5 * k * (theta - theta0)**2      

      AngleEnergy = ZERO
      AngleEnergyF = ZERO
      AngleEnergyB = ZERO
                       
      ivariable = .false.
      if (ityp.eq.PROTEIN) then
        iprot = ival
c       don't bother evaluating angle if angle within fixed residues
        ires1 = Zid2ResID(iprot,Angle(iprot,iang,1),1)
        ires2 = Zid2ResID(iprot,Angle(iprot,iang,2),1)
        ires3 = Zid2ResID(iprot,Angle(iprot,iang,3),1)
        if (FixedRes(iprot,ires1).and.FixedRes(iprot,ires2).and.
     .         FixedRes(iprot,ires3)) then
           return
        endif
        iparam = Angle(iprot,iang,4)
        iangval = AngleSize(iprot,iang)
        iangvalf = iangval
        iangvalb = iangval
      else if (ityp.eq.SOLUTE) then
        isol = ival
        iparam = SolAngle(isol,iang,4)
        iangval = SolAngleSize(isol,iang)
        if (NSolVariable(isol).gt.0 .and. DifferentLambdas) then
          ivariable = .true.
          iangvalf = SolAngleSize(NSolForwards(isol),iang)
          iangvalb = SolAngleSize(NSolBackwards(isol),iang)
        else
          iangvalf = iangval
          iangvalb = iangval
        endif
      else
        iprot = 0
        isol = 0
        iparam = 0
        iangval = ZERO                                     
        iangvalf = iangval
        iangvalb = iangval
        write(printstring,*) "Unknown system type requested - ",ityp
        call printLine(WARNING,printstring)
        call printLine(WARNING,"This is probably a programming bug!")
        return
      endif

      if (iparam.lt.0) then
C       variable parameter angle
        ik = VariableAngParam(-iparam,2,1)
        ir = VariableAngParam(-iparam,2,2)
      else if (iparam.eq.0) then
        angleEnergy = ZERO
        angleEnergyB = ZERO
        angleEnergyF = ZERO
        return
      else
        ik = AngParam(iparam,1)
        ir = AngParam(iparam,2)
      endif

      angleEnergy = ik  * (iangval-ir)**2
     
      if (DifferentLambdas) then
        if (iparam.gt.0 .and. .not.ivariable) then
          angleEnergyB = angleEnergy
          angleEnergyF = angleEnergy
        else
          if (iparam.lt.0) then
            ikb = VariableAngParam(-iparam,4,1)
            irb = VariableAngParam(-iparam,4,2)
            ikf = VariableAngParam(-iparam,5,1)
            irf = VariableAngParam(-iparam,5,2)
          else
            ikf = ik
            ikb = ik
            irf = ir
            irb = ir
          endif
          angleEnergyB = ikb  * (iangvalb-irb)**2
          angleEnergyF = ikf  * (iangvalf-irf)**2
        endif
      else
        angleEnergyB = angleEnergy
        angleEnergyF = angleEnergy
      endif
      
c$$$      write(printstring,*) 'iang ',iang,' ik ',ik,' ir ',ir,
c$$$     &  ' ikb ',ikb,' irb ',irb,' ikf ',ikf,' irf ',irf
c$$$      call printLine(DEBUG,printstring)
c$$$
c$$$      write(printstring,*) ' iangval ',iangval,' iangvalF ',iangvalf,
c$$$     &  ' iangvalB ',iangvalb
c$$$      call printLine(DEBUG,printstring)
c$$$      write(printstring,*) 'nrg ',angleEnergy,' nrgf ',angleEnergyF,
c$$$     &  ' nrgB ',angleEnergyB
c$$$      call printLine(DEBUG,printstring)
      return
      
      end
