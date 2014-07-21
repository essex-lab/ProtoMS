#!/bin/env python

import sys
import os
import string
import re

#dictionary of all of the functions
funcs = {}
#dictionary of all of the subroutines
subs = {}
#dictionary of all of the routines
routines = {}

#dictionary of all of the global variables
variables = {}
#dictionary of parameters
parameters = {}

#bc string
bcstr = []

class Variable:
    """This class holds all of the information about a variable"""

    def __init__(self,name,typ,filename=None):
        self.name = name.strip().lstrip().rstrip()
        self.type = typ
        self.dimensions = None
        self.value = None
        self.origval = None
        
        variables[name.lower()] = self
        
        if (filename is None):
            self.file = None
        else:
            self.file = filename

    def setValue(self,val):
        val = val.lower().replace("**","^")
        self.value = val
        self.origval = self.value
        fullval = self.expandValue()
        self.value = fullval
        
        fullval = "%s = %s;" % (self.name.lower().replace("\\",""),fullval)

        bcstr.append(fullval)

    def expandValue(self):
        """Expand this parameter's value"""
        if (self.value.isdigit()): return self.value
        
        echo = "%s %s" % (string.join(bcstr," "),self.value.lower().replace("**","^"))

        #if (echo.lstrip().rstrip() == ""):
        #    return 1.0
    
        cmd  = "echo '%s' | bc 2> /dev/null" % echo

        pipe = os.popen(cmd,"r")
        line = pipe.readline().rstrip("\n").strip().rstrip().lstrip()
        
        if (line.isdigit()): return line
        
        return 1.0

    def setDimensions(self,dims):
        self.dimensions = dims

    def size(self):
        base = 32
        if (self.type == "logical"): base = 1
        elif (self.type.find("character") == 0):
            tmp = self.type.split("*")
            if (len(tmp) > 1):
                base = int(tmp[-1])*8
            else:
                base = 8
        elif (self.type == "double precision"):
            base = 64

        #now get the dimensions
        if (self.dimensions is None): return base
        
        for dim in self.dimensions:
            #try to convert the dimension to numbers
            if dim.isdigit():
                base = base * int(dim)
            elif dim.lower() in parameters:
                base = base * int(parameters[dim.lower()].value)
            else:
                #see if this is a combo
                dims = dim.split("+")
                total = 0
                for dim in dims:
                    if dim.isdigit():
                        total = total + int(dim)
                    elif dim.lower() in parameters:
                        total = total + int(parameters[dim.lower()].value)
                        
                if (total > 0):
                    base = base * total
            
        return base

    def stringDimensions(self):
        dims = []
        if (self.dimensions is None):
            return "1"
            
        for dim in self.dimensions:
            if (dim.isdigit()):
                dims.append(dim)
            else:
                dims.append("%s" % (dim))

        return string.join(dims,", ")

    def numberDimensions(self):
        dims = []
        if (self.dimensions is None):
            return "1"
            
        for dim in self.dimensions:
            if dim.isdigit():
                dims.append(str(dim))
            elif dim.lower() in parameters:
                dims.append(str(parameters[dim.lower()].value))
            else:
                #see if this is a combo
                ds = dim.split("+")
                tmp = None
                for dim in ds:
                    if (dim.isdigit()):
                        if tmp is None: tmp = str(dim)
                        else: tmp = "%s + %s" % (tmp,str(dim))
                    elif dim.lower() in parameters:
                        if tmp is None: tmp = str(parameters[dim.lower()].value)
                        else: tmp = "%s + %s" % (tmp,str(parameters[dim.lower()].value))
                    else:
                        if tmp is None: tmp = str(dim)
                        else: tmp = "%s + %s" % (tmp,dim)

                dims.append(tmp)
                
        return string.join(dims,", ")

class Routine:
    """This class holds all of the information about a function or subroutine
    
    Attributes:
        __isfunc - whether or not this is a function
        __name   - the name of this routine
        __calls  - list of all of the routines that this routine calls
        __calledby - list of all of the routines that call this routine
        __reads  - list of all of the global variables that this routine reads
        __writes - list of all of the global variables that this routine modifies
        __comment - the comment associated with this routine
        __args   - arguments to the subroutine
        __returns - return type of the function
        lines    - range of lines that this routine is on
        file     - the file that this routine is from
        intrinsic - whether or not this is an intrinsic routine
    
    """
    def __init__(self,name,func=False,returns=None):
        self.__isfunc = func
        self.__name = name.replace("_","\_")
        self.__calls = []
        self.__calledby = []
        self.__reads = []
        self.__writes = []
        self.__args = []
        self.__comment = None
        self.__returns = returns
        self.lines = [0,0]
        self.file = None
        self.intrinsic = False
        
        routines[self.__name.lower()] = self
        
        if (func): funcs[self.__name.lower()] = self
        else: subs[self.__name.lower()] = self
        
    def isFunction(self):
        return self.__isfunc
        
    def isSubroutine(self):
        return not self.__isfunc

    def name(self):
        return self.__name

    def declaration(self):
        if (self.__isfunc):
            if (len(self.__args) == 0):
                return "%s function %s()" % (self.__returns,self.__name)
            else:
                return "%s function %s(%s)" % (self.__returns,self.__name, \
                                                string.join(self.__args,",").rstrip(","))
        else:
            if (len(self.__args) == 0):
                return "subroutine %s()" % (self.__name)
            else:
                return "subroutine %s(%s)" % (self.__name, \
                                                string.join(self.__args,",").rstrip(","))
    
    def calledBy(self):
        return self.__calledby
        
    def calls(self):
        return self.__calls

    def addCalledBy(self,name):
        name = name.lower()
    
        if (name in self.__calledby): return
        
#        if (not name in routines):
#            print "Cannot find called-by routine ",name
            
        self.__calledby.append(name)

    def addCalls(self,name):
        name = name.lower()
        
        if (name in self.__calls): return
        
        if (not name in routines):
            #print "Cannot find called routine ",name
            return
            
        self.__calls.append(name)
        routines[name].addCalledBy(self.__name)

    def addReads(self,name):
        name = name.lower()

        if (name not in variables): return
        
        val = variables[name]
        if (val in self.__reads): return
        
        self.__reads.append(val)
        
    def addWrites(self,name):
        name = name.lower()

        if (name not in variables): return

        val = variables[name]
        if (val in self.__writes): return
        
        self.__writes.append(val)
        addReads(self,name)                

    def cleanReadWrite(self):
        """Cleans lists, as if writes variable then must also read variable - 
           remove the double entry"""
        for val in self.__writes:
            if val in self.__reads: self.__reads.remove(val)

    def reads(self):
        return self.__reads
        
    def writes(self):
        return self.__writes
    
    def setArgs(self,line):
        """Try to find contents of '(...)'"""
        words = line.split("(")
        if len(words) < 2: return
        
        #args should be in words[1] - split into variables
        words = words[1].strip().rstrip(")").split(",")
        self.__args = words
        
    def __get_comment(self):
        return self.__comment
    
    def __set_comment(self,text):
        self.__comment = text.replace("#"," ")
        
    comment = property(__get_comment, __set_comment, None, "Comment")

def runCmd(cmd):
    """Runs command 'cmd' and returns the stdout as a list of lines"""
    
    pipe = os.popen(cmd,"r")
    lines = pipe.readlines()
    pipe.close()
    
    return lines
    
def readFile(f):
    """ Reads the file 'f' and returns a list of all of the lines"""
    
    try:
        lines = open(f,"r").readlines()
        return lines;
    except IOError:
        print "Cannot open %s" % f
        return []

def isComment(line):
    if (line[0:1].lower() == "c"): return True
    else: return False
    
def isContinue(line):
    if isComment(line): return False
    elif (len(line) < 6): return False
    elif (line[6:7] == "."): return True
    else: return False
    
def routineName(word):
    name = word.split("(")[0]
#    name = name.replace("_",'\_')
    return name
    
def parseInclude(f):
    """ Parses a fortran include file"""

    lines = readFile(f)

    i = 0
    nlines = len(lines)
    
    while (i < nlines):
        line = lines[i]
        
        #see if we need to join together continued lines
        if (i < len(lines)-1 and isContinue(lines[i+1])):
            line = line = lines[i+1][6:]
            i = i + 1
        
        #increment line    
        i = i + 1    
        
        #skip all comment lines
        if (isComment(line)): continue
        
        try:
            words = line.split()
            #skip short lines
            if (len(words) < 1): continue
    
            #if the first word is 'common' then skip the line
            if (words[0].lower() == "common"): 
                continue
            elif (words[0].lower().find("parameter") != -1):
                #get the parameter
                param = line.split("(")[1].split("=")[0]
                value = line.split("=")[1].rstrip(")\n")
                if (param == "MAXATOMS"): value = value + ")"
                
                #try to find this variable
                if (param.lower() in variables):
                    v = variables[param.lower()]
                    parameters[param.lower()] = v
                    
                    #expand the value
                    v.setValue(value)
                
                #now add this parameter
            elif (words[0].lower() == "integer" or words[0].lower() == "logical" or \
                  words[0].lower().find("character") == 0 or \
                  words[0].lower() == "double" and words[1].lower() == "precision"):
                #this is a global variable (or global variables!)
                valtyp = words[0].lower()
                strt = 1
                
                if (words[0].lower() == "double"):
                    valtyp = "double precision"
                    strt = 2
                    
                #get all of the variables on the line
                line = string.join(words[strt:]," ")
                #split by ')' - this splits arrays up by arrays
                mvals = line.split(")")
                for mval in mvals:
                    dims = mval.split("(")
                    vals = dims[0].split(",")
                    
                    for val in vals:
                        val = val.strip().rstrip().lstrip()
                        var = Variable(val,valtyp,f)
                        
                    #if dims[1] then the last variable was an array
                    if (len(dims) == 1): continue
                    
                    dims = dims[1].split(",")
                    
                    variables[val.lower()].setDimensions(dims)
                        
        except IndexError:
            print "Index error!"
            return
        except KeyError, (e):
            print "Key error for line ",line.rstrip()
            print e,e.args
            continue
        

if (__name__ == "__main__"):
    
    #get all the files in the current directory
    files = os.listdir("./")
    files.sort()

    #read the dimensions.inc file
    parseInclude("dimensions.inc")

    #parse the include files first
    for f in files:
        if (f.endswith(".inc")):
            #skip 'nb*.inc' files
            if (f.find("nb") == 0 or f == "dimensions.inc"): continue
            #this is an include file!
            parseInclude(f)

    #now print out the memory usage
    total = 0.0    
    totals = {}
    print "There are %d variables" % len(variables)
    for v in variables:
        variable = variables[v]
        sz = variable.size()/(1024*8.0)
        szstr = "   %8.2f kB   " % sz
        total = total + sz
        if (sz > 5500):
            szstr = "** %8.2f MB **" % (float(sz) / (1024.0))
        elif (sz > 1024):
            szstr = "   %8.2f MB   " % (float(sz) / (1024.0))
            
        if (sz >= 1):
            tmp = "%s\t%s(%s) %s (%s)" % (szstr,variable.name, \
                      variable.stringDimensions(),variable.type,variable.numberDimensions())
                
            if (sz in totals):
                totals[sz] = "%s\n%s" % (totals[sz],tmp)
            else:
                totals[sz] = tmp

    keys = totals.keys()
    keys.sort()
    for key in keys:
        print totals[key]
                                         
    total = total / (1024.0)
    print "\nNote that the above values are estimates. Based on these\n" \
          "values, ProtoMS would require %4.0f MB if fully loaded." % total

