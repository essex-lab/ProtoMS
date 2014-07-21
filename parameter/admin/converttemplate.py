#!/bin/env python2

# This script is used to provide a first approximation of a translation
# of a ProtoMS 1 template into a ProtoMS 2 template. The result will
# still need to be significantly hand-modified, though this script does
# remove most of the donkey work!

import sys

restemps = {}
soltemps = {}
svntemps = {}

skipatoms = ('n','ca','c','hn','o','ha')

class ResTemplate:
    
    def __init__(self,nam):
        self.name = nam
        
        self.atoms = []
        self.zmats = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        
        #add to the global dictionary
        restemps[nam] = self

    def output(self):
        print "mode template"
        print "residue %s" % self.name
        print "info rotate 3.0 translate 0.5"
        print "backbone first aanterm middle aacenter last aacterm single aasingle"
        for atm in self.atoms:
            print "atom  %3s  %3s  %3s  %3s" % (atm)

        for zmat in self.zmats:
            print "zmat  %3s  %5.3f  %5.1f  %6.1f" % (zmat)

        for bond in self.bonds:
            print "bond  %3s  %3s  flex %5.2f" % (bond)
            
        for angle in self.angles:
            print "angle  %3s  %3s  %3s  flex %5.2f" % (angle)
            
        for dihedral in self.dihedrals:
            print "dihedral  %3s  %3s  %3s  %3s  flex %5.2f" % (dihedral)

def readFile(f):
    return open(f,"r").readlines()
    
def checkMode(words):
    if (words[0].lower() == "residue"):
        #new residue template
        tmpl = ResTemplate(words[1])
        return (True,"residue",tmpl)
    elif (words[0].lower() == "solute"):
        return (True,None,None)
    elif (words[0].lower() == "solvent"):
        return (True,None,None)
    else:
        return (False,None,None)
    
def checkSubMode(words):
    if (words[0].lower() == "atoms"):
        return (True,"atoms")
    elif (words[0].lower() == "bonds"):
        return (True,"bonds")
    elif (words[0].lower() == "angles"):
        return (True,"angles")
    elif (words[0].lower() == "dihedrals"):
        return (True,"dihedrals")
    else:
        return (False,None)
    
lines = readFile(sys.argv[1])

mode = None
tmpl = None
submode = None

for line in lines:
    if (line.find("#") == 0): continue
    
    words = line.split()
    if (len(words) <= 0): continue
    
    (change,m,t) = checkMode(words)
    if (change):
        mode = m
        tmpl = t
        submode = None
        continue

    if (mode is None): continue
    elif (mode == "residue"):
        #reading in a resdiue to 'tmpl' object - first
        #get the mode we are in within the residue template
        (change,s) = checkSubMode(words)
        if (change):
            submode = s
            continue
        
        if (len(line) < 26): continue
        
        if (submode == "atoms"):
            if (len(words) < 13):
                print "Cannot recognise line ",line.rstrip()
                continue
            
            #skip backbone atoms
            if (words[0].lower() in skipatoms): continue
            
            tmpl.atoms.append((words[0],words[4],words[7],words[10]))
            tmpl.zmats.append((words[0],float(words[6]), \
                              float(words[9]),float(words[12])))
        elif (submode == "bonds"):
            #skip implicit bonds
            if (line[18:19] == '0'): continue
            else:
                tmpl.bonds.append((line[2:5],line[6:9],float(line[26:])))
        elif (submode == "angles"):
            #skip implicit angles
            if (line[18:19] == '0'): continue
            else:
                tmpl.angles.append((line[2:5],line[6:9],line[10:13], \
                                   float(line[26:])))
        elif (submode == "dihedrals"):
            #skip implicit dihedrals
            if (line[18:19] == '0'): continue
            else:
                tmpl.dihedrals.append((line[2:5],line[6:9], \
                                       line[10:13],line[14:17], \
                                       float(line[26:])))
            
# now write out each residue template in the new format
keys = restemps.keys()
keys.sort()

for key in keys:
    restemps[key].output()
    print " "
    
