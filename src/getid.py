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
user = pwd.getpwuid(os.getuid())[0]
today = datetime.datetime.today()
host = socket.gethostname()

writeLine("Compiled on %s by %s at %s." % (host,user,today))

