c-----------------------------------------------------
c  Variables needed by inlined function nbenergy.inc
c
c     DO NOT TOUCH THESE VARIABLES IN YOUR CODE!!!

c     these are the input variables..
      integer funciprot1,funciprot2
      integer funcisol1,funcisol2
      integer funcisvn1,funcisvn2
      integer funcigcs1,funcigcs2
      integer funcat1,funcat2
      integer funcitmp,funcitmp1,funcitmp2
      integer ipar1,ipar2

c     these are the variables that are useful in your subroutine(output)
      double precision funcnbnrg,funccnrg,funcljnrg,funcmindist
      double precision funcnbnrgf,funcnbnrgb,funcdistf2,funcdistb2,funcdistf6,funcdistb6
      double precision invfuncdistf,invfuncdist6f,invfuncdist12f
      double precision invfuncdistb,invfuncdist6b,invfuncdist12b
      double precision funccnrgf,funccnrgb,funcljnrgf,funcljnrgb
      
c     do not use these variables...
      double precision funcichg1,funcichg2,funcisig1,funcisig2,funcieps1,funcieps2
      double precision funcichg1f,funcichg2f,funcisig1f,funcisig2f,funcieps1f,funcieps2f
      double precision funcichg1b,funcichg2b,funcisig1b,funcisig2b,funcieps1b,funcieps2b
      double precision invfuncdist,funcdist2,funcdist6,invfuncdist6,invfuncdist12
      double precision funcdx,funcdy,funcdz
      double precision funcsqrteps,funcsig3
      double precision funcsqrtepsf,funcsig3f
      double precision funcsqrtepsb,funcsig3b
      double precision funcdxf,funcdyf,funcdzf,funcdxb,funcdyb,funcdzb
c     JM mod
      double precision qiqj,qiqjf,qiqjb
      double precision funcaij,funcaijf,funcaijb
      double precision funcbij,funcbijf,funcbijb
      double precision softdist,softdistf,softdistb
      double precision sigat1,sigat2,thet
      logical boundx,boundy,boundz,funcdogeom,funcsigsclf,funcsigsclb
      logical softcore,ddd,pertchg                                
      double precision fac,facf,facb
      double precision shift,shiftf,shiftb
      common /NBCOMMON/ funcnbnrg,funccnrg,funcljnrg,funcmindist,
     .                  funcnbnrgf,funcnbnrgb,funccnrgf,funccnrgb,funcljnrgf,funcljnrgb,
     .                  funcichg1,funcichg2,funcisig1,funcisig2,funcieps1,funcieps2,
     .                  invfuncdist6f,invfuncdist12f,invfuncdist6b,invfuncdist12b,
     .                  funcichg1f,funcichg2f,funcisig1f,funcisig2f,funcieps1f,funcieps2f,
     .                  funcichg1b,funcichg2b,funcisig1b,funcisig2b,funcieps1b,funcieps2b,
     .                  invfuncdist,funcdist2,funcdist6,invfuncdist6,invfuncdist12,
     .                  funcdxf,funcdyf,funcdzf,funcdxb,funcdyb,funcdzb,funcdistf2,
     .                  funcdistb2,funcdistf6,funcdistb6,funcdx,funcdy,funcdz,funcsqrteps,funcsig3,
     .                  funcsqrtepsf,funcsig3f,qiqj,qiqjf,qiqjb,softdist,
     .                  softdistf,softdistb,sigat1,sigat2,
     .                      fac,facf,facb,shift,shiftf,shiftb,thet,
     .                  funcaij,funcaijf,funcaijb,funcbij,funcbijf,funcbijb,
     .                  funcsqrtepsb,funcsig3b,funcitmp,funcitmp1,funcitmp2,
     .                  funciprot1,funciprot2,funcisol1,funcisol2,funcisvn1,funcisvn2,funcigcs1,funcigcs2,
     .                  funcat1,funcat2,ipar1,ipar2,
     .                  boundx,boundy,boundz,funcdogeom,funcsigsclf,funcsigsclb,softcore,ddd,pertchg
