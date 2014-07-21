
import sys

for file in sys.argv[1:]:
    lines = open(file, "r").readlines()
    FILE = open(file, "w")

    for line in lines:
        line = line.replace("\t", "      ")
        print >>FILE,line,

    FILE.close()
