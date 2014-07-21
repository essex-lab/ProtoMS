#!/usr/bin/env python2.3

import os,sys
try:
    inFile = sys.argv[1]
except IndexError:
    print "USAGE is %s ProtoMC1.0 template file" % (sys.argv[0])
    sys.exit(-1)

stream = open(inFile,'r')
buffer = stream.readlines()
stream.close()
inAtoms = False
inResidue = False
backbone = ['N','CA','C','O','HN','HA','HA1']
FirstresLines = []
MiddleresLines = []
LastresLines = []
print "mode override\n"
for line in buffer:
    if line.startswith('#'):
        continue
    if line.startswith('RESIDUE'):
        inResidue = True
        elems = line.split()
        residue = elems[1]
        continue
    if line.startswith('ATOMS'):
        inAtoms= True
        continue
    if line.startswith('BONDS'):
        inAtoms= False
        # Now print data
        print "residue %s %s" % (residue,'first')
        print '\n'.join(FirstresLines)
        print ' '
        print "residue %s %s" % (residue,'middle')
        print '\n'.join(MiddleresLines)
        print ' '
        print "residue %s %s" % (residue,'last')
        print '\n'.join(LastresLines)
        print ' '
##         print "chain aacenter%s " % residue
##         print '\n'.join(MLines)
##         print "chain aanterm%s " % residue
##         print '\n'.join(NTLines)
##         print "chain aacterm%s " % residue
##         print '\n'.join(CTLines)
        #sys.exit(-1)
        # clear
        FirstresLines = []
        MiddleresLines = []
        LastresLines = []
        MLines = []
        NTLines = []
        CTLines = []
        continue
    if inAtoms and inResidue:
        elems = line.split()
        atom = elems[0]
        if (atom == 'HN'):
            line = "atom HN1 %d" % (int(elems[2]))
            FirstresLines.append(line)
            line = "atom HN2 %d" % (int(elems[2]))
            FirstresLines.append(line)
            line = "atom HN3 %d" % (int(elems[2]))
            FirstresLines.append(line)
            line = "atom HN %d " % (int(elems[1]))
            MiddleresLines.append(line)
            line = "atom HN %d " % (int(elems[3]))
            LastresLines.append(line)              
        elif (atom == 'O'):
            line = "atom O %d" % (int(elems[3]))
            LastresLines.append(line)            
            line = "atom OT %d "%  (int(elems[3]))
            LastresLines.append(line)
            line = "atom O %d" % (int(elems[2]))
            FirstresLines.append(line)
            line = "atom O %d" % (int(elems[1]))
            MiddleresLines.append(line)            
        else:
            line = "atom %s %d " % (atom,int(elems[2]))
            FirstresLines.append(line)
            line = "atom %s %d " % (atom,int(elems[1]))
            MiddleresLines.append(line)
            line = "atom %s %d " % (atom,int(elems[3]))
            LastresLines.append(line)        

#        print line.strip()
        
        

