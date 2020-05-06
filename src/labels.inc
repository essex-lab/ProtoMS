c#############################################
c
c  This file contains the variables
c  describing labels (or 'id')
c
c#############################################
      integer MAXID
      parameter(MAXID=1000) 
      integer LabelID(MAXID) 
      integer LabelType(MAXID) 
      integer LabelMolecule(MAXID) 
      character*4 LabelResidueName(MAXID) 
      integer LabelResidueNumber(MAXID)
      character*4 LabelName(MAXID)
      common /ILABELDATA/ LabelID,LabelType,LabelMolecule,LabelName,
     .                    LabelResidueName,LabelResidueNumber
    



