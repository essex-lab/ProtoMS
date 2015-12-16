"""This test is to check the Sampling MC moves."""

import nose
import unittest
#import argparse
import os
import sys
import subprocess

import tools
from tools import simulationobjects as sim

from subprocess import call

# ---------------------------------------------
# ProtoMS Equilibration Sampling MC moves test
# ---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
ref_dir = proto_env + "/tests/energies/"

#Energy components stored by SnapshotResults objects
comps = [ 'interaction_energies', 'internal_energies', 'capenergy', 'extraenergy' ]

class TestSampling(unittest.TestCase):

    """Test for component energies output by ProtoMS."""

    def test_energies(self):


        if call("$PROTOMSHOME/build/protoms3 run_t4p.cmd", shell=True):
            raise Exception ( "Energies test simulation failed!" )

        self.compare_results ( 'out_t4p/results', '%s/out_t4p/results' % ref_dir,
                               't4p' )


        if call("$PROTOMSHOME/build/protoms3 run_t3p.cmd", shell=True):
            raise Exception ( "Energies test simulation failed!" )

        self.compare_results ( 'out_t3p/results', '%s/out_t3p/results' % ref_dir,
                               't3p' )

    def compare_results ( self, res, ref, label ):
        results = sim.SnapshotResults ()
        with open ( res ) as f:
            results.parse ( f )

        ref_results = sim.SnapshotResults ()
        with open ( ref ) as f:
            ref_results.parse ( f )

        for i in comps:
            comp = getattr ( results, i )
            ref_comp = getattr ( ref_results, i )
            # where comp is a dict
            try:
                for key in comp:
                    for e, ref_e in zip ( comp[key], ref_comp[key] ):
                        if e.curr != ref_e.curr:
                            raise Exception ( "%s energy of %s in %s of %s simulation is incorrect"
                                              % ( e.type, key, i, label ) )
                continue
            except TypeError:
                pass

            # where comp is an EnergyResults object
            try:
                if comp.curr != ref_comp.curr:
                    raise Exception ( "%s energy of in %s simulation is incorrect"
                                      % ( comp.type, label ) )
                continue
            except AttributeError:
                pass

            # where comp is an integer
            if comp != ref_comp:
                raise Exception ( "%s in %s simulation is incorrect" % ( i, label ) )

# Entry point for nosetests or unittests
if __name__ == '__main__':
    unittest.main()
    nose.runmodule()
