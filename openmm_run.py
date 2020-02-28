from __future__ import print_function
import argparse
import sys
import os

from omm_readinputs import *
from omm_readparams import *
from omm_vfswitch import *
from omm_barostat import *
from omm_restraints import *
from omm_hmr import *
from omm_rewrap import *
from omm_run import *

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=True)
parser.add_argument('-p', dest='psffile', help='Input CHARMM PSF file', required=True)
parser.add_argument('-c', dest='crdfile', help='Input CHARMM CRD file', required=True)
parser.add_argument('-t', dest='toppar', help='Input CHARMM-GUI toppar stream file', required=True)
parser.add_argument('-b', dest='sysinfo', help='Input CHARMM-GUI sysinfo stream file (optional)', default=None)
parser.add_argument('-icrst', metavar='RSTFILE', dest='icrst', help='Input CHARMM RST file (optional)', default=None)
parser.add_argument('-irst', metavar='RSTFILE', dest='irst', help='Input restart file (optional)', default=None)
parser.add_argument('-ichk', metavar='CHKFILE', dest='ichk', help='Input checkpoint file (optional)', default=None)
parser.add_argument('-opdb', metavar='PDBFILE', dest='opdb', help='Output PDB file (optional)', default=None)
parser.add_argument('-orst', metavar='RSTFILE', dest='orst', help='Output restart file (optional)', default=None)
parser.add_argument('-ochk', metavar='CHKFILE', dest='ochk', help='Output checkpoint file (optional)', default=None)
parser.add_argument('-odcd', metavar='DCDFILE', dest='odcd', help='Output trajectory file (optional)', default=None)
parser.add_argument('-hmr', dest='hmr', help='Apply Hydrogen Mass Repartitioning (optional)', action='store_true', default=False)
parser.add_argument('-rewrap', dest='rewrap', help='Re-wrap the coordinates in a molecular basis (optional)', action='store_true', default=False)
args = parser.parse_args()

urun = OpenMM_Run()
urun.run_with_args(**args)
