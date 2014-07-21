#!/bin/env python

import sys
import string

lines = open(sys.argv[1],"r").readlines()

chains = {}
residues = {}

cur = None
for line in lines:
    words = line.split()
    if (len(words) < 2): continue
    
    if (words[0] == "chain"):
        chains[words[1]] = {}
        cur = chains[words[1]]
        continue
    elif (words[0] == "residue"):
        residues[words[1]] = {}
        cur = residues[words[1]]
        continue
        
    if (cur is None): continue
    
    if (words[0] == "atom"):
        cur[words[1]] = [words[2],words[3]]
        
lines = open(sys.argv[2],"r").readlines()

res = None
for line in lines:
    line = line.rstrip()

    words = line.split()
    if (len(words) < 2): 
        print line
        continue
    elif (line[0:1] == '#'): 
        print line
        continue
        
    if (words[0] == "chain"):
        res = words[1]
        if (res not in chains):
            print "#WARNING: %s not in overrides" % res
            res = None
        else:
            res = chains[res]
            
        print line
        continue
    elif (words[0] == "residue"):
        res = words[1]
        if (res not in residues):
            print "#WARNING: %s not in overrides" % res
            res = None
        else:
            res = residues[res]
            
        print line
        continue
        
    if (res is None):
        print line
        continue
        
    if ( (words[0] != "atom") and (words[0] != "bbatom")): 
        print line
        continue
        
    if (words[0] == "bbatom"):
        if (words[2] not in res):
            print "#WARNING: %s not in residue!" % words[2]
            print line
            continue
        
        atm = res[words[2]]
        print line," %4s %4s" % (atm[0],atm[1])
    elif (words[0] == "atom"):
        if (words[1] not in res):
            print "#WARNING: %s not in residue!" % words[1]
            print line
            continue

        atm = res[words[1]]
        print "%s %4s  %4s %4s  %s" % (words[0],words[1],atm[0],atm[1],string.join(words[2:]," "))
