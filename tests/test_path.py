#!/usr/bin/python2.7

"""
Tests for ProtoMS python module (protoms.py): Setup tests
"""

#---------------------------------------------
# Imports
#---------------------------------------------

import os
import unittest
import sys
import nose
import nose.tools as nt
import subprocess
import re
#---------------------------------------------
# Setup tests
#---------------------------------------------

# Storing PROTOMSHOME environment variable to a python variable.
proto_env = os.environ["PROTOMSHOME"]
proto_path = os.popen("pwd").read()
proto_path = re.sub('\\n$', '', proto_path)


class TestSetUp(unittest.TestCase):

# Test if PROTOMSHOME is set properly. PROTOMSHOME must be set after protoms install at HOME directory.
    def setUp(self):
        super(TestSetUp, self).setUp()
        print("Setting up PROTOMSHOME check test case.") 

    def tearDown(self):
        super(TestSetUp, self).tearDown()
        print("Cleaning up PROTOMSHOME check test case.")

    def test_protoms_path(self):
        
        try:
            os.path.exists(proto_path)
            os.getenv("PROTOMSHOME")
        
        except KeyError, e:
            print e
            print("PROTOMSHOME is not set.")

     
        try:
            self.assertEqual(proto_env,proto_path,msg=None)

        except:
            print("PROTOMSHOME is not set properly.")


if  __name__ == '__main__':
    unittest.main()
