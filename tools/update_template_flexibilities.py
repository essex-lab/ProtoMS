import numpy as np
import sys
sys.path.append ( '/home/chris/protoms_tmp/tools' )
import simulationobjects as sim
import argparse

parser = argparse.ArgumentParser ( description = 'Script to parse move histograms from ProtoMS dihedral tuning'
                                   'simulations to generate new template files.' )
parser.add_argument ( '-a', '--accept_file', required = True, help = 'simulation accept file containing histograms' )
parser.add_argument ( '-t', '--tem_file', required = True, help = 'template file to update with new flexibilities' )
parser.add_argument ( '-o', '--output', required = True, help = 'output template file' )
args = parser.parse_args()

with open ( args.accept_file ) as f:
    while not f.next().startswith ( 'Move histograms' ):
        pass
    hists = []
    try:
        while 1:
            while not f.next().startswith ( 'Dihedral' ):
                pass
            f.next()
            hists.append ( [ float ( f.next() ) for i in range ( 20 ) ] )
            
    except StopIteration:
        pass

tem = sim.TemplateFile ( )
tem.read ( args.tem_file )
dihs = [ i for i in tem.templates[0].connectivity if i.type == 'dihedral' ]

try:
    s = 0.5 + 0.5 / len ( hists )
    hists = np.array ( hists )
    for j, hist in enumerate ( hists ):
        c = np.cumsum ( hist )
        for i in range ( len ( c ) ):
            c[i] = c[i] / ( i + 1.0 )

        for i in range ( len ( c ) ):
            if c[i] < (0.6): #**(1/4.)):
                break
        print ( i + 1 ) * s
        dihs[j].flex = ( i + 1 ) * s
except ZeroDivisionError:
    print("No Rotatable dihedrals template not modified")
    
tem.write(args.output)
