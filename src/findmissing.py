
import sys

for file in sys.argv[1:]:
    ofile = file.replace(".F", ".o")

    try:
        open(ofile, "r")
    except:
        print "Missing %s" % file
