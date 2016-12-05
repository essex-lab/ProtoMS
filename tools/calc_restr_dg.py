# Authors: Richard Bradshaw
#          Ana Cabedo Martinez
#          Chris Cave-Ayland
#          Samuel Genheden
#          Gregory Ross

"""
Program to calculate free energies using TI, BAR and MBAR
"""

import os
import logging

import numpy as np

import simulationobjects
import calc_ti
import calc_bar
import pms2pymbar

import matplotlib.pyplot as plt

logger = logging.getLogger('protoms')

import argparse

# Setup a parser of the command-line arguments
parser = argparse.ArgumentParser(description="Program to calculate free energy of releasing a Harmonic restraint using the Zwanzig equation.")
parser.add_argument('-f','--file',help='the instantaneous results file at the desired lambda value.')
parser.add_argument('-s','--skip',type=int,help="the number of blocks to skip to calculate the free energy differences in one window. default is 0. Skip must be greater or equal to 0",default=0)
parser.add_argument('-m','--max',type=int,help="the upper block to use. default is 99999 which should make sure you will use all the available blocks. max must be greater or equal to 0",default=99999)
parser.add_argument('-t','--temperature',type=float,help="the simulation temperature in degrees. Default is 25.0",default=25.0)
args = parser.parse_args()

# Setup the logger
logger = simulationobjects.setup_logger()

# Fix negative values of skip and max
if args.max < 0 :
  args.max = 99999
if args.skip <= 0 :
  args.skip = -1

RT = 1.9872041*(args.temperature+273.15)/1000

rf = simulationobjects.ResultsFile ()
rf.read(filename=args.file)
rf.make_series()

print -RT * np.log ( np.exp ( rf.series.harmonic.curr / RT ).mean() )  



  
