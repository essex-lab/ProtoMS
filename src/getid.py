#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import pwd
import datetime
import socket


def writeLine(line):
    """Write an ID line to the HEADER stream"""

    lgth = len(line)

    istrt = 0
    idfile.write("      write(printstring,*)\n")

    while istrt < lgth:
        iend = istrt + 64
        if iend > lgth:
            iend = lgth

        idfile.write('     . "%s"' % line[istrt:iend])
        istrt = istrt + 64
        if istrt < lgth:
            idfile.write(",\n")
        else:
            idfile.write("\n")

    idfile.write("      call printBoxText(HEADER,printstring)\n")


try:
    filename = sys.argv[1]
except IndexError:
    print("USAGE: getid.py outputfile")
    sys.exit(0)

# open the include file
try:
    idfile = open(filename, "w")
except:
    print("Cannot open id file %s for writing" % idfile)
    sys.exit(0)

# now write down the date that this code was compiled, and who has compiled it
user = pwd.getpwuid(os.getuid())[0]
today = datetime.datetime.today()
today = today.replace(microsecond=0)
host = socket.gethostname()

writeLine("Compiled on %s by %s at %s." % (host, user, today))

try:
    identify = os.popen("hg identify", "r").readline()
    identify = identify.strip("\n ")
    path = os.popen("hg paths default", "r").readline()
    path = path.strip("\n ")
except:
    writeLine("Mercurial information not available.")
else:
    if identify.startswith("abort") or path == "not found!":
        writeLine("Mercurial information not available.")
    else:
        writeLine("Revision: %s" % identify)
        if "@" in path:
            writeLine("From: %s" % path[path.find("@") + 1 :])
