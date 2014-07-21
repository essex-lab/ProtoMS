#!/usr/bin/python

import sys,os,pwd
import datetime,socket

def writeLine(line):
    """Write an ID line to the HEADER stream"""
    
    lgth = len(line)
    
    istrt = 0
    idfile.write("      write(printstring,*)\n")
    
    while (istrt < lgth):
        iend = istrt + 64
        if (iend > lgth): iend = lgth 
        
        idfile.write("     . \"%s\"" % line[istrt:iend])
        istrt = istrt + 64
        if (istrt < lgth): idfile.write(",\n")
        else: idfile.write("\n")

    idfile.write("      call printBoxText(HEADER,printstring)\n")

try:
    filename = sys.argv[1]
except:
    print "USAGE: getid.py outputfile"
    sys.exit(0)

#open the include file
try:
    idfile = open(filename,"w")
except:
    print "Cannot open id file %s for writing" % idfile
    sys.exit(0)

#now write down the date that this code was compiled, and who has compiled it
#user = os.getlogin()
user = pwd.getpwuid(os.getuid())[0]
today = datetime.datetime.today()
host = socket.gethostname()

writeLine("Compiled on %s by %s at %s." % (host,user,today))

try:
    lines = os.popen("svn info ../","r").readlines()
except:
   writeLine("Subversion information not available on the computer "\
             "that compiled this version of ProtoMC.")
else:
    svninfo = {}
    for line in lines:
        words = line.split()
        if (len(words) < 2): continue
        
        if (words[0].find("URL") != -1):
            svninfo["url"] = words[1]
        elif (words[0].find("Revision") != -1):
            svninfo["rev"] = words[1]
    
    try:
        writeLine("Subversion version %s. If you have access to the " \
                  "subversion respository you can extract the source code "\
                  "for this binary by typing" % svninfo["rev"])
        writeLine("'svn co -r %s %s'." % (svninfo["url"],svninfo["rev"]))
    except:
        writeLine("Subversion information not available on the computer "\
                  "that compiled this version of ProtoMS.")
    
