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
        self.name = name.replace("_","\_").strip().lstrip().rstrip()
        self.type = typ
        self.dimensions = None
        self.value = None
        self.origval = None
        
        if (filename is None):
            self.file = None
        else:
            self.file = filename

    def setValue(self,val):
        self.value = val.replace("(","").replace(")","")
        self.origval = self.value
        fullval = self.expandValue()
        self.value = fullval
        
        fullval = "%s = %s;" % (self.name.lower().replace("\\",""),fullval)
        
        bcstr.append(fullval)

    def expandValue(self):
        """Expand this parameter's value"""
        if (self.value.isdigit()): return self.value
        
        cmd  = "echo '%s %s' | bc" % (string.join(bcstr," "),self.value.lower())
        pipe = os.popen(cmd,"r")
        line = pipe.readline()
        
        if (line.isdigit()): return line
        
        return 1

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
        
        if (not name in routines):
            print "Cannot find called-by routine ",name
            
        self.__calledby.append(name)

    def addCalls(self,name):
        name = name.lower()
        
        if (name in self.__calls): return
        
        if (not name in routines):
            print "Cannot find called routine ",name
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
    
    lines = open(f,"r").readlines()
    return lines;

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
                value = line.split("=")[1].split(")")[0]
                
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
                        variables[val.lower()] = Variable(val,valtyp,f)
                        
                    #if dims[1] then the last variable was an array
                    if (len(dims) == 1): continue
                    
                    dims = dims[1].split(",")
                    
                    variables[val.lower()].setDimensions(dims)
                        
                  
                #see if the previous lines were comments - if so then
                #assume these comments are associated with this variable
                #j = i - 2
                #comment = None
                #while (j >= 0):
                #    if (isComment(lines[j])):
                #        if (comment is None):
                #            comment = lines[j][1:].strip().lstrip().rstrip()
                #        else:
                #            comment = lines[j][1:].strip().lstrip().rstrip() + " " + comment
                #        j = j - 1
                #    else:
                #        break

        except IndexError:
            print "Index error!"
            return
        except KeyError, (e):
            print "Key error for line ",line.rstrip()
            print e,e.args
            continue
    

def parseFortran(f,mode=1):
    """Parses a fortran source file, using the specified mode (mode 1 is used
       to get all of the routine names and comments, while mode 2 is
       used to perform an in-depth scan of the source code one all of the
       function names have been resolved."""

    lines = readFile(f)
    
    routine = None
    
    i = 0
    nlines = len(lines)
    
    while(i < nlines):
        line = lines[i]
        
        #see if this line continues onto the next line - if so then
        #join the lines together
        if (i < len(lines)-1 and isContinue(lines[i+1])):
            line = line + lines[i+1][6:]
            i = i + 1
        
        i = i + 1
        try:
            #split the line into words
            words = line.split()
            if (len(words) <= 0): continue
        
            #if we are not in a routine then we need to find one!
            if (routine is None):
                if (isComment(line) or isContinue(line)): continue
                
        
                if (words[0].lower() == "subroutine" or words[0].lower() == "program"):
                    #we have found a new subroutine! - get the name
                    name = routineName(words[1])
                    if (mode == 1):
                        routine = Routine(name)
                        routine.setArgs(line)
                    else:
                        routine = routines[name.lower()]
                elif (words[1].lower() == "function"):
                    #we have found a new function! - get the name
                    name = routineName(words[2])
                    if (mode == 1):
                        routine = Routine(name,True,words[0].lower())
                        routine.setArgs(line)
                    else:
                        routine = routines[name.lower()]
                elif (words[2].lower() == "function"):
                    #we have found a double precision function!
                    name = routineName(words[3])
                    if (mode == 1):
                        routine = Routine(name,True,"double precision")
                        routine.setArgs(line)
                    else:
                        routine = routines[name.lower()]

                if (mode == 1):
                    if (not routine is None):
                        #we have just created the routine - save the file and first line number
                        routine.file = f
                        routine.lines[0] = i+1

            else:
                if (isComment(line)):
                    #only do this in mode 1
                    if (mode == 1): continue
                    
                    #do not need to look at comments if this routine already has a comment
                    if (not routine.comment is None): continue
                     
                    #see if this is a routine comment
                    if (line.find("#############") != -1):
                        #read the lines to get the comment
                        j = i+1
                        comment = None
                        while (j < nlines):
                            if (not isComment(lines[j]) or lines[j].find("############") != -1): break
                            
                            l = lines[j][1:].strip().lstrip().rstrip()
                            
                            if (comment is None):
                                comment = l + "\\\\ "
                            else:
                                comment = comment + l + "\\\\ "
                                
                            j = j + 1
                        
                        #set the line number
                        i = j
                        
                        #set the routine's comment
                        if (not comment is None):
                            routine.comment = comment.replace("#"," ").replace("_","\_")
                                                
                        continue

                #if the first word is 'end' then the routine is ended
                if (words[0].lower() == "end"): 
                    #save the last line of the routine
                    if (mode == 1): routine.lines[1] = i+1
                    routine = None

                #if we are in mode 1 then we have finished our work
                if (mode == 1): continue

                #look at each line - what is it doing? - is a subroutine called?
                if (line.lower().find("call") != -1):
                    for j in range(0,len(words)):
                        if (words[j].lower() == "call"):
                            name = routineName(words[j+1])
                            routine.addCalls(name.lower())
                            break
                        
                if (words[0].lower() == "logical" or \
                    words[0].lower() == "integer" or \
                    words[0].lower().find("character") != -1):
                 
                    #this is a variable declaration - if this variable
                    #has the same name as a function, then this routine is assumed
                    #to call this function.
                    vals = string.join(words[1:]," ").split(",")
                    for val in vals:
                        if (val.lower() in funcs):
                            routine.addCalls(val)
                            
                elif (words[0].lower() == "double" and words[1].lower() == "precision"):
                
                    #this is a double precision variable
                    vals = string.join(words[2:]," ").split(",")
                    for val in vals:
                        if (val.lower() in funcs):
                            routine.addCalls(val)
                            
                else:
                    #one variable is being assigned to another - open the line up,
                    #(e.g. replace "(" with " ( ", "+" with " + " etc.
                    
                    rexp = re.compile(",|\+|\-|\)|\(|\*|\/|\.")
                    line = rexp.sub(" ",line)                    

                    #split the line into words
                    words = line.split()
                    #check each word to see if it is a global variable
                    for word in words:
                        if (word.lower() in variables):
                            routine.addReads(word.lower())

        except IndexError:
            print "Index error!"
            return
        except KeyError, (e):
            print "Key error for line ",line.rstrip()
            print e,e.args
            continue
            
def writeVariables(vals,f):
    """Write the variables in the list to a table in the filehandle f"""
    
    if (len(vals) == 0):
        f.write("Nothing.\\\\ \n")
        return

    sortvals = []
    
    for val in vals:
        sortvals.append(val.name)
        
    sortvals.sort()
        
    f.write('\\begin{longtable}{l l l}\n')
    
    n = len(sortvals)
    nrows = (n/3)
    
    if (3*nrows < n): nrows = nrows + 1
    
    for i in range(0,nrows):
        (x,y,z) = (i*3,i*3+1,i*3+2)
        
        for j in (x,y):
            if (j < n):
                f.write("%s & " % sortvals[j])
            else:
                f.write(" & ")
        
        if (z < n):
            f.write("%s \\\\ \n" % sortvals[z])
        else:
            f.write(" \\\\ \n")        
        
    f.write('\end{longtable}')
            
def writeRoutineTable(names,f):
    """Write the routines in the list to a table in the filehandle f"""
    
    if (len(names) == 0):
        f.write("Nothing.\\\\ \n")
        return

    names.sort()
    
    n = len(names)
    nrows = (n / 3)
    
    if (3*nrows < n): nrows = nrows + 1
    
    f.write('\\begin{longtable}{l l l}\n')
    for i in range(0,nrows):
        (x,y,z) = (i*3,i*3+1,i*3+2)
        
        for j in (x,y):
            if (j < n):
                if (names[j].lower() in routines): name = routines[names[j].lower()].name()
                else: name = names[j]
            
                f.write("\\hyperlink{%s}{%s} & " % (name,name))
            else:
                f.write(" & ")
                
        if (z < n):
            if (names[z].lower() in routines): name = routines[names[z].lower()].name()
            else: name = names[z]
            f.write("\\hyperlink{%s}{%s} \\\\ \n" % (name,name))
        else:
            f.write("\\\\ \n")

    f.write('\end{longtable}\n')

            
def writeTexRoutine(routine,f):
    """Write out the routine 'routine' to filehandle f, starting on a newpage if newpage is true"""

    #start each routine on a new page    
    f.write("\\newpage \n")
    
    #write out the declaration, and hyperlink target    
    f.write("\hypertarget{%s}{\Large{\\textbf{%s}}}" % (routine.name(), \
                               routine.declaration().replace(",",", ")))
    
    f.write("\n\n\\normalsize\n")
    
    #write out the file that holds this routine
    f.write("Declared in lines %d to %d of \\file{%s}.\n" % (routine.lines[0],routine.lines[1],routine.file))

    f.write("\\newline\n")

    #now write out the comment
    f.write("\\rule{6in}{0.5mm}\\newline\n")
    if (routine.comment is None):
        f.write("\\texttt{No source level documentation is available.\\\\}\n")
    else:
        f.write("\\texttt{%s}\n" % routine.comment)
    f.write("\\rule{6in}{0.5mm}\\newline\\newline\n")

    #now write out all the routines that call this one
    f.write("This is called by; \\\\ \n")
    writeRoutineTable(routine.calledBy(),f)

    #now write out all of the routines that this one calls
    f.write("This calls; \\\\ \n")
    writeRoutineTable(routine.calls(),f)

    #clean the reads/write list
    routine.cleanReadWrite()

    #now write out all of the variables that this routine writes
    f.write("This uses the following global variables; \\\\ \n")
    writeVariables(routine.reads(),f)
    

def writeTexOutput(filedir):
    """Write out all of the latex output to 'filename'"""
    f = open("%s/autodocs.tex" % filedir,"w")
    
    f.write('\section{Subroutines and Functions}\label{sec:devman-routines}\n')

    #write out the names of all of the subroutines and functions in a 3-column table
    f.write('\pms is composed of the following subroutines and functions.\n\n')
    
    writeRoutineTable(routines.keys(),f)

    names = subs.keys()
    names.sort()
    for name in names:
        writeTexRoutine(subs[name],f)
        
    names = funcs.keys()
    names.sort()
    for name in names:
        writeTexRoutine(funcs[name],f)

    f.close()

def writeDotOutput(filedir,topfunc,texfile,depth=None):

    if (not topfunc.lower() in routines): 
        print "Cannot find ",topfunc
        return

    texfile.write('\\begin{figure}[H]\centering\n')
    texfile.write("\\includegraphics[height=5in]{%s-dot}\n" % topfunc)

    if (depth is None):
        texfile.write("\caption{Call graph for %s.}" % topfunc)
    else:
        texfile.write("\caption{Call graph for %s (truncated at level %d).}" % (topfunc,depth))
        
    texfile.write("\\label{fig:devman-callgraph-%s}\n" % topfunc)
    texfile.write('\end{figure}\n')
    
    routine = routines[topfunc.lower()]

    f = open("%s/%s.dot" % (filedir,routine.name()),"w")
    
    f.write("digraph objtree {\nratio=fill\nsize=\"6.0,5.0\"\n")

    done = {}

    if (len(routine.calledBy()) > 0):
       for r in routine.calledBy():
           f.write("%s -> %s\n" % (routines[r].name(),routine.name()))

    writeDotFunc(routine,f,done,depth)
    
    f.write("}\n")
    f.close()

def writeDotFunc(routine,f,done=None,depth=None):
    
    if (done is None):
        done = {}
        
    if (not depth is None):
        if (depth <= 0): return done
        
    if (routine in done): return done
    else: done[routine] = 1
    
    if (len(routine.calls()) > 0):
       if (not depth is None): 
           depth = depth - 1

       for r in routine.calls():
           f.write("%s -> %s\n" % (routine.name(),routines[r].name()))
           done = writeDotFunc(routines[r],f,done,depth)

    return done
    

if (__name__ == "__main__"):
    
    #get all the files in the current directory
    files = os.listdir("./")
    files.sort()

    #parse the include files first
    for f in files:
        if (f.endswith(".inc")):
            #skip 'nb*.inc' files
            if (f.find("nb") == 0): continue
            #this is an include file!
            parseInclude(f)

    #now parse the fortran files
    for f in files:
        if (f.endswith(".F")):
            #this is a fortran file!
            parseFortran(f)

    #now parse each fortran file performing a deep read of the source code
    # (this can only be performed once all of the files have been read in mode 1)
    for f in files:
        if (f.endswith(".F")):
            #this is a fortran file!
            parseFortran(f,2)
            
    print "There are %d subroutines and %d functions" % (len(subs),len(funcs))

    #see if we should write this all out
    if (len(sys.argv) < 2): sys.exit(0)

    writeTexOutput(sys.argv[1])
    
    f = open("%s/callgraphs.tex" % sys.argv[1],"w")
    
    callgraphs = {"totalEnergy": None, \
                  "residueEnergy": None, \
                  "soluteEnergy": None, \
                  "solventEnergy": None, \
                  "residueMove": None, \
                  "soluteMove": None, \
                  "solventMove": None, \
                  "volumeMove": None, \
                  "readParFile": 2, \
                  "readBnd": None, \
                  "readAng": None, \
                  "readDih": None, \
                  "readUB": None, \
                  "readCLJ": None, \
                  "readTmpl": None, \
                  "assignSystem": 2, \
                  "assignSoluteTemplates": None, \
                  "assignSolventTemplates": None, \
                  "buildProtein": None, \
                  "loadProtein": None, \
                  "loadSolute": None, \
                  "loadSolvent": None, \
                  "simulation": 1, \
                  "getOption": None, \
                  "simulate": 2, \
                  "equilibrate": 2, \
                  "readRestart": None, \
                  "printPDB": None, \
                 }
    
    f.write("""\section{Call Graphs}\label{sec:devman-callgraphs}
               This section contains call graphs for key subroutines and functions
               (henceforth called routines). The call graph for a routine
               shows which routines call it, and which routines it calls.
               The call graph will then show which routines all called routines
               call. This graph can get very complicated for a top-level
               routine, so higher level graphs are truncated after 
               a certain depth (e.g. a depth of one would only show the 
               routines called by this routine, while a depth of two would
               show all routines call by this routine, or called by routines
               called by this routine. Fortran 77 does not allow recursion,
               so the control flow through a call graph will always be from
               top to bottom.\n\nThe call graphs for the following routines
               are shown; \n\n \\begin{tabularx}{6in}{l l l l l l}""")

    keys = callgraphs.keys()
    keys.sort()
    
    nkeys = len(keys) / 3
    if (3*nkeys < len(keys)): nkeys = nkeys + 1
    
    for i in range(0,nkeys):
        j = 3*i

        if (j >= len(keys)): break
        graph = keys[j]
        f.write("%s & Figure \\ref{fig:devman-callgraph-%s}" % (graph,graph))
        
        j = j + 1
        if (j >= len(keys)):
            f.write(" & & ")
        else:
            graph = keys[j]
            f.write(" & %s & \\ref{fig:devman-callgraph-%s} " % (graph,graph))

        j = j + 1
        if (j >= len(keys)):
            f.write(" & & \\\\ \n")
        else:
            graph = keys[j]
            f.write(" & %s & \\ref{fig:devman-callgraph-%s} \\\\ \n" % (graph,graph))

    f.write("\end{tabularx}\n")

    f.write('\\begin{landscape}\n')
    f.write('\pagestyle{empty}\n')

    for graph in callgraphs:
        writeDotOutput(sys.argv[1],graph,f,callgraphs[graph])
    
    f.write('\end{landscape}\n')
    f.write('\pagestyle{fancy}\n')

    f.write("""\section{Dimensions and Memory Usage}\label{sec:devman-dimensions}
               This section lists the default values of the dimensions and parameters used in
               the program, and an approximation of the memory used by the 
               global arrays.""")

    params = parameters.keys()
    params.sort()
    
    f.write("""\n\n\\begin{longtable}{|l|l|}\hline\n
                Parameter & Default Value \\\\ \hline \n""")
    
    for param in params:
        p = parameters[param]
        f.write("\hypertarget{%s}{%s} & %s \\\\ \n" % (p.name.lower(),p.name,p.origval))
        
    f.write("""\hline\end{longtable}""")
            
            
    f.write("""\n\n\\begin{longtable}{|l|l|r|}\hline\n
                 Global Variable & Type & Memory Usage \\\\ \hline \n""")
    
    total = 0.0    
    for v in variables:
        variable = variables[v]
        sz = variable.size() / 1024.0
        szstr = "%8.2f kB" % sz
        total = total + sz
        if (sz > 5500):
            szstr = "\\textbf{%8.2f MB}" % (float(sz) / 1024.0)
        elif (sz > 1024):
            szstr = "%8.2f MB" % (float(sz) / 1024.0)
            
        if (sz > 300):
            dims = []
            for dim in variable.dimensions:
                if (dim.isdigit()):
                    dims.append(dim)
                else:
                   dims.append("\hyperlink{%s}{%s}" % (dim.lower(),dim))
                
            f.write("\parbox{3in}{%s(%s)} & %s & %s \\\\ \hline \n" % (variable.name, \
                      string.join(dims,", "),variable.type,szstr))
                                         
    total = total / 1024.0
    f.write("\hline\end{longtable}\n\n Note that the above values are estimates. Based on these \
               values, \pms would require %4.0f MB if fully loaded (in practise the actual \
               required memory is much less than this)." % total)

    f.close()            

