#!/bin/env python

## This script takes a set of old-style ProtoMC parameter files
## (bond,angle,dihedral,clj,template) and converts them into a
## single new style parameter file.

import os
import sys
import string

###
### Functions used in the script
###

def readlines(filename):
    try:
        if (filename is None): return []
        f = open(filename,"r")
        lines = f.readlines()
        return lines
    except:
        print "Problem reading file '%s'" % filename
        return []

def readBond(filename):
    """Read all of the angle file and get the angle parameters"""
    lines = readlines(filename)

    params = {}
    atms = []
    
    stage = 1
    for line in lines:
        line = line.rstrip()
        if line[0:1] == '#': continue
	elif line[0:3].lower() == "end": continue
    
        words = line.split()
	if (len(words) < 1):
	    stage = 2
	    continue
    
        if stage == 1:    
            #format is par  k  r  ! comment
            (ipar,k,r) = words[0:3] 
            comment = string.join(words[4:]," ")
            
            params[int(ipar)] = {'k' : k, 'r' : r, 'comment' : comment}
        else:
            #format is XX-XX par
            atom = line[0:5]
            atoms = atom.split("-")
            par = int(line.split(" ")[-1])
            atms.append({"atms" : atoms, "par" : par})
    
    return (params,atms)
    
def readAngle(filename):
    """Read all of the angle file and get the angle parameters"""
    lines = readlines(filename)
    
    #get the parameters...
    params = {}
    atms = []

    stage = 1      
    for line in lines:
        line = line.rstrip()
        if (line[0:1] == "#"): continue
        words = line.split()
        if (len(words) < 1):
            stage = 2
            continue
        
	if stage == 1:
            #format is par  k  theta  ! comment
            (ipar,k,theta) = words[0:3]
            comment = string.join(words[4:]," ")
        
            params[int(ipar)] = { 'k' : k, 'theta' : theta, 'comment' : comment}
        else:
            #format is XX-XX-XX  par
            if (line[0:3].lower()  == "end"): continue
            atom = line[0:8]
            atoms = atom.split("-")
            par = int(line.split(" ")[-1])
            atms.append({ "atms" : atoms, "par" : par})
                
    return (params,atms)


def readDihedral(filename, mode='OPLS'):
    """Read all of the dihedral file and get the dihedral parameters"""
    lines = readlines(filename)
    
    #get the parameters...
    params = {}
    atms = []

    stage = 1      
    for line in lines:
        line = line.rstrip()
        if (line[0:1] == "#"): continue
        words = line.split()
        if (len(words) < 1):
            stage = 2
            continue

	if (stage ==1):
            if mode == 'AMBER':
                # format is par v0 f0 v1 f1 v2 f2 v3 f3 v4 f4 ! comment
                (par,v0,f0,v1,f1,v2,f2,v3,f3,v4,f4) = words[0:11]
                comment = string.join(words[12:]," ")
                f0 = float(f0) * 90.0
                f1 = float(f1) * 90.0
                f2 = float(f2) * 90.0
                f3 = float(f3) * 90.0
                f4 = float(f4) * 90.0
                v0 = float(v0)
                v1 = float(v1) 
                v2 = float(v2) 
                v3 = float(v3)
                v4 = float(v4)
        
                # save params, converting phases to degrees (from units of 90 degs)
                # and v1/v2/v3/v4
                params[int(par)] = { 'v0' : v0, 'f0' : f0, \
                                     'v1' : v1, 'f1' : f1, \
                                     'v2' : v2, 'f2' : f2, \
                                     'v3' : v3, 'f3' : f3, \
                                     'v4' : v4, 'f4' : f4, \
                                     'comment' : comment}
            else:
                # format is par v0 f0 v1 f1 v2 f2 v3 f3  ! comment
                (par,v0,f0,v1,f1,v2,f2,v3,f3) = words[0:9]
                comment = string.join(words[11:]," ")
                f0 = float(f0) * 90.0
                f1 = float(f1)*90.0
                f2 = float(f2)*90.0
                f3 = float(f3)*90.0
                v1 = float(v1)*0.5
                v2 = float(v2)*0.5
                v3 = float(v3)*0.5
        
                # save params, converting phases to degrees (from units of 90 degs)
                # and v1/v2/v3 divided by two
                params[int(par)] = { 'v0' : v0, 'f0' : f0, \
                                     'v1' : v1, 'f1' : f1, \
                                     'v2' : v2, 'f2' : f2, \
                                     'v3' : v3, 'f3' : f3, \
                                     'comment' : comment}
        else:
            #format is XX-XX-XX-XX par scl comment
            if (line[0:3].lower()  == "end" ): continue
            elems = line.split(" ")
            if (len(elems) == 1): continue
            atom = line[0:11]
            atoms = atom.split("-")
            if (len(atoms) != 4):
                atoms = ["  ","  ","  ","  "]
            par = int(line[12:16])
            if mode == 'AMBER':
                comment = line[24:]
            else:
                comment = line[20:]
            atms.append({ "atms" : atoms, "par" : par, 'comment' : comment})
                
    return (params,atms)
    
def readCLJ(filename, mode='OPLS'):
    """Read in all of the charge/Lennard Jones parameters"""
    lines = readlines(filename)
    
    params = {}
    rstosig = 0.8908987181
    
    i = 0
    for line in lines:
        i = i + 1
        line = line.rstrip()
        #skip the first two lines
        if (i <= 2): continue
        #format is par protnum atm charge sigma epsilon comments
        if line.startswith('END'): continue
        words = line.split()
	if (words[0].lower() == "end"): break
        par = int(words[0])
        #skip comment lines
        if (par == 999): continue
        # get the parameters
        (prot,atm,chg,sig,eps) = words[1:6]
        if mode == 'AMBER':
            sig = float(sig) * rstosig * 2.00
        #get the comment
        # comment = line[43:]
        comment = line[40:]
        params[par] = { 'atm' : atm, 'protnum' : prot, \
                        'charge' : chg, 'sigma' : sig, \
                        'epsilon' : eps, 'comment' : comment }
    
    return params
    
def readTemplate(filename):
    """Read in all of the templates"""
    lines = readlines(filename)
    print "#Need to implement conversion of template files!"

def outputBonds(params,atms):
    """Output bond parameters in new format"""
    
    #output the bond parameters first
    print "\nmode bond"
    print "# U(r) = k(r-r0)**2"
    print "#parameter k(kcal mol-1 A-2) r0(A) comment"
    for par in params:
        param = params[par]
        print "par %5s   %9.3f     %8.4f  #%s" % (par,float(param["k"]),float(param["r"]),param["comment"])

    #now output the atms
    print "#atm atm1 atm2 parameter"
    doneatoms = {}
    for atm in atms:
        atoms = atm["atms"]
        key = string.join(atoms,"-")
        if (key in doneatoms):
            print >> sys.stderr, "WARNING! Duplicated bond atoms! ",key
            continue
        else: doneatoms[key] = 1
        print "atm %3s  %3s    %5d" % (atoms[0],atoms[1],atm["par"])

def outputAngles(params,atms):
    """Output angle parameters in new format"""
    
    #output the angle parameters first
    print "\nmode angle"
    print "# U(theta) = k(theta-theta0)**2"
    print "#parameter k(kcal mol-1 deg-2) theta0(deg) comment"
    for par in params:
        param = params[par]
        print "par %5s   %9.3f     %8.4f  #%s" % (par,float(param["k"]),float(param["theta"]),param["comment"])

    #now output the atms
    print "#atm atm1 atm2 atm3 parameter"
    doneatoms = {}
    for atm in atms:
        atoms = atm["atms"]
        key = string.join(atoms,"-")
        if (key in doneatoms): 
            print >> sys.stderr, "WARNING! Duplicated angle atoms! ",key
            continue
        else: doneatoms[key] = 1
        
        print "atm %3s  %3s  %3s    %5d" % (atoms[0],atoms[1],atoms[2],atm["par"])

def outputDihedrals(params,atms, mode='OPLS'):
    """Output dihedral parameters in new format"""
    
    #output the dihedral terms first - this is more complex as these
    #parameters are split into individual cosine terms
    print "\nmode dihedral"
    print "# U(phi) = k1*( 1.0 + k2*(cos[k3*phi + k4]) )"
    print "#term k1(kcal mol-1) k2 k3 k4(deg) #comment"
    trans = {}
    i = 1
    if mode == 'AMBER':
        # Need to make sense of all this stuff...
        for par in params:
            trans[par] = "%5d %5d %5d %5d %5d" % (i,i+1,i+2,i+3,i+4)
            dih = params[par]
## VERY WRONG !
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i,float(dih["v0"]),0.0,0.0,float(dih["f0"]))
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+1,float(dih["v1"]),1.0,1.0,float(dih["f1"]))
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+2,float(dih["v2"]),-1.0,2.0,float(dih["f2"]))
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+3,float(dih["v3"]),1.0,3.0,float(dih["f3"]))
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+4,float(dih["v4"]),-1.0,4.0,float(dih["f4"]))
## A BIT WRONG
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i,float(dih["v0"]),0.0,0.0,0.0)
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+1,float(dih["v1"]),1.0,1.0,0.0)
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+2,float(dih["v2"]),-1.0,2.0,0.0)
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+3,float(dih["v3"]),1.0,3.0,0.0)
##             print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+4,float(dih["v4"]),-1.0,4.0,0.0)
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i,float(dih["v0"]),0.0,0.0,float(dih["f0"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+1,float(dih["v1"]),1.0,1.0,float(dih["f1"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+2,float(dih["v2"]),1.0,2.0,float(dih["f2"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+3,float(dih["v3"]),1.0,3.0,float(dih["f3"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+4,float(dih["v4"]),1.0,4.0,float(dih["f4"]))            
            i = i+5
    else:
        for par in params:
            trans[par] = "%5d %5d %5d %5d" % (i,i+1,i+2,i+3)
            dih = params[par]
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i,float(dih["v0"]),0.0,0.0,float(dih["f0"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+1,float(dih["v1"]),1.0,1.0,float(dih["f1"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+2,float(dih["v2"]),-1.0,2.0,float(dih["f2"]))
            print "term %5d  %8.3f %8.3f %8.3f %8.3f" % (i+3,float(dih["v3"]),1.0,3.0,float(dih["f3"]))
            i = i+4

    #now the parameters
    print "#par  term1  term2 etc..  #comment"
    for par in params:
        dih = params[par]
        print "par %5d %s  # %s" % (par,trans[par],dih["comment"])

    #now output the atms
    print "#atm atm1 atm2 atm3 atm4 parameter #comment"
    doneatms = {}
    for atm in atms:
        atoms = atm["atms"]
        key = string.join(atoms,"-")
        if (key in doneatms):
            print >> sys.stderr, "WARNING! Duplicated dihedral atoms ",key
            continue
        else: doneatms[key] = 1
        par = atm["par"]
        print "atm  %3s  %3s  %3s  %3s %5d      #%s" % \
            (atoms[0],atoms[1],atoms[2],atoms[3],\
             par,atm["comment"])

def outputCLJs(params):
    """Print out all of the charge/Lennard Jones parameters"""
    print "\nmode clj"
    print "# Energy is sum of coulomb and Lennard Jones terms"
    print "#parameter atm proton-num charge(|e|) sigma(A) epsilon(kcal mol-1) #comment"
    for par in params:
        clj = params[par]
        print "par %5s   %2s   %3s   %9.4f %9.4f %9.4f  #%s" % \
              (par,clj["atm"],clj["protnum"],float(clj["charge"]),float(clj["sigma"]), \
               float(clj["epsilon"]),clj["comment"])
    

###
### Parse the options and read the files
###

import getopt

longopts = ["bond=","angle=","dihedral=","clj=","template=","ff="]

try:
    opts = getopt.getopt(sys.argv[1:],[],longopts)
except getopt.GetoptError:
    print "Error parsing options!"
except:
    pass
   
files = {}
files["bond"] = None
files["angle"] = None
files["dihedral"] = None
files["clj"] = None
files["template"] = None
files['ff'] = 'AMBER'

for option in opts[0]:
    typ = option[0].lstrip("--")
    files[typ] = option[1]
   
#read in all of the parameters
(bond,bndatms) = readBond(files["bond"])
(angle,angatms) = readAngle(files["angle"])
(dihedral,dihatms) = readDihedral(files["dihedral"], mode=files['ff'])
cljs = readCLJ(files["clj"], mode=files['ff'])

readTemplate(files["template"])

#now output the parameters to stdout in the correct format
print "#Automatically generated new style ProtoMC template file"
outputBonds(bond,bndatms)
outputAngles(angle,angatms)
outputDihedrals(dihedral,dihatms, mode=files['ff'])
outputCLJs(cljs)
