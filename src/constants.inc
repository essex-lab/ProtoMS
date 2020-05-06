
c#############################################
c
c  This file contains the physical constants
c  that are used in ProtoMC
c
c#############################################

c physical constants - for units of kcal mol-1, deg C, etc.
c     pi, and multiples of pi
      double precision PI
      parameter(PI=3.14159265358979323846264D+00)
      double precision TWOPI
      parameter(TWOPI=PI*2.0D+00)
      double precision PIOTWO
      parameter(PIOTWO=PI*0.5D+00)

c     quick and easy access to useful numbers
      double precision ZERO,ONE,TWO,THREE,FOUR,NINETY,ONEEIGHTY,THREESIXTY
      double precision HALF,THIRD,QUARTER
      double precision SMALL,LARGE,TINY,HUGE
      parameter(ZERO=0.D+00)
      parameter(ONE=1.D+00)
      parameter(TWO=2.D+00)
      parameter(THREE=3.D+00)
      parameter(FOUR=4.D+00)
      parameter(NINETY=90.D+00)
      parameter(ONEEIGHTY=180.D+00)
      parameter(THREESIXTY=360.D+00)
      parameter(HALF=ONE/TWO)
      parameter(THIRD=ONE/THREE)
      parameter(QUARTER=ONE/FOUR)
      parameter(SMALL=0.001D+00)
      parameter(TINY=1.D-8)
      parameter(LARGE=99999.0D+00)
      parameter(HUGE=1.D+8)

c     convert between Angstroms and metres
      double precision ANGST2METRE
      parameter(ANGST2METRE=1.0D-10)
      double precision METRE2ANGST
      parameter(METRE2ANGST=1.0D+10)

c     convert kJ, J and kcal
      double precision CAL2J
      parameter(CAL2J=4.184D+00)
      double precision J2CAL
      parameter(J2CAL=1.0D+00/CAL2J)
      double precision KCAL2J
      parameter(KCAL2J=1000.0D+00*CAL2J)
      double precision J2KCAL
      parameter(J2KCAL=1.0D+00/KCAL2J)

c     degrees to kelvin
      double precision DEG2KEL
      parameter(DEG2KEL=273.15D+00)

c     avagadros number      
      double precision AVAGADRO
      parameter(AVAGADRO=6.02214179d23)
      double precision ONEOVERAVAGADRO
      parameter(ONEOVERAVAGADRO=1.0D+00/AVAGADRO)
c     Boltzmann constant in joules K-1
      double precision KBJ
      parameter(KBJ=1.380650324d-23)
c     Boltzmann constant in kcal K-1
      double precision KB
      parameter(KB=KBJ*J2KCAL)
c     Gas constant in kcal mol-1 K-1
      double precision GASR
      parameter(GASR=KB*AVAGADRO)
c     Gas constant in J mol-1 K-1
      double precision GASRJ
      parameter(GASRJ=GASR*KCAL2J)
c     Planck constant in kg m2 s-1
      double precision PLANCK
      parameter (PLANCK=6.626068D-34)

c     atmospheres to pascals
      double precision ATM2PASCAL
      parameter(ATM2PASCAL=101325D+00)
c     volume in angstroms^3 to volume in metre^3 per mol
      double precision ANGSTVOL2METPERMOL
      parameter(ANGSTVOL2METPERMOL=ANGST2METRE**3*AVAGADRO)
c     atmospheres to internal pressure units - when multiplied by volume
c     in A^3 will return energy in kcal mol-1
      double precision P2INTERNALP
      parameter(P2INTERNALP=J2KCAL*ATM2PASCAL*ANGSTVOL2METPERMOL)
      
c     Convert unit charges to Coulombs
      double precision E2C
      parameter(E2C=1.60217646263d-19)
c     convert Coulombs to unit charges
      double precision C2E
      parameter(C2E=1.0D+00 / E2C)

c     epsilon0 in Coulomb^2 J-1 m-1
      double precision EPS0CJM
      parameter(EPS0CJM=8.854187817d-12)
c     epsilon0 in UnitCharge^2 J-1 A-1 mol
      double precision EPS0J
      parameter(EPS0J=EPS0CJM*C2E*ONEOVERAVAGADRO*C2E*ANGST2METRE)
c     epsilon0 in UnitCharge^2 kcal-1 A-1 mol
      double precision EPSILON0
      parameter(EPSILON0=EPS0J/J2KCAL)
c     4PiEps0 in UnitCharge^2 kcal-1 A-1 mol
      double precision FOURPIEPS0
      parameter(FOURPIEPS0=4.0D+00*PI*EPSILON0)
c     1 / 4PiEps0 in UnitCharge-1 kcal A mol-1
      double precision FAC4PIEPS0
      parameter(FAC4PIEPS0=1.0D+00/FOURPIEPS0)
