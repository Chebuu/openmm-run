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

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *


class OpenMM_Run:
    def __init__(self):
        pass

    def load_params(self, inpfile, toppar, psffile, crdfile, sysinfo=None):
        print("Loading parameters")
        psf = read_psf(psffile)
        self.inputs = read_inputs(inpfile)
        self.params = read_params(toppar)
        self.crd = read_crd(crdfile)
        if sysinfo:
            self.psf = read_box(psf, args.sysinfo)
        else:
            self.psf = gen_box(psf, crd)
        return self.params

    def build_system(self):
        self.system = self.psf.createSystem(self.params, nonbondedMethod=self.inputs.coulomb, 
                                                        nonbondedCutoff=self.inputs.r_off*nanometers,
                                                        constraints=self.inputs.cons,
                                                        ewaldErrorTolerance=self.inputs.ewald_Tol)
        return self.system

    def build_integator(self, hmr=None):
        if self.inputs.vdw == 'Force-switch': self.system = vfswitch(self.system, self.psf, self.inputs)
        if self.inputs.pcouple == 'yes':      self.system = barostat(self.system, self.inputs)
        if self.inputs.rest == 'yes':         self.system = restraints(self.system, self.crd, self.inputs)
        if hmr:                               self.system = HydrogenMassRepartition(self.system, self.psf)
        self.integrator = LangevinIntegrator(inputs.temp*kelvin, self.inputs.fric_coeff/picosecond, 
                                                                self.inputs.dt*picoseconds)
        return self.integrator

    def set_platform(self, useCUDA=False):
        def set_CPU():
            top_platform = 0
            n_platforms = Platform.getNumPlatforms()
            for i in range(0, n_platforms-1):
                curr_platform = Platform.getPlatform(top_platform)
                curr_speed = curr_platform.getSpeed()
                idx_platform = Platform.getPlatform(i)
                idx_speed = idx_platform.getSpeed()
                if idx_speed > curr_speed:
                    top_platform = i
            platform = Platform.getPlatform(top_platform)
            return platform
        if useCUDA: 
            try:
                self.platform = Platform.getPlatformByName('CUDA')
                self.prop = dict(CudaPrecision='single')
            except:
                self.platform = set_CPU()
                self.prop = None
        else:
            self.platform = set_CPU()
            self.prop = None
        return self.platform

    def build_context(self, icrst=None, irst=None, ichk=None):
        self.simulation = Simulation(self.psf.topology, self.system, self.integrator, self.platform, self.self.prop)
        self.simulation.context.setPositions(self.crd.positions)
        if None not in (icrst, irst, ichk):
            if icrst:
                charmm_rst = read_charmm_rst(icrst)
                self.simulation.context.setPositions(charmm_rst.positions)
                self.simulation.context.setVelocities(charmm_rst.velocities)
                self.simulation.context.setPeriodicBoxVectors(charmm_rst.box[0], charmm_rst.box[1], charmm_rst.box[2])
            if irst:
                with open(irst, 'r') as f:
                    self.simulation.context.setState(XmlSerializer.deserialize(f.read()))
            if ichk:
                with open(ichk, 'rb') as f:
                    self.simulation.context.loadCheckpoint(f.read())
        else:
            raise ValueError("At least one argument to 'icrst', 'irst', or 'ichk' must be passed.")
        return self.simulation

    def rewrap(self, simulation=None):
        if simulation is None: simulation=self.simulation
        self.simulation = rewrap(simulation)
        return self.simulation

    def get_system_energy():
        return self.simulation.context.getState(getEnergy=True).getPotentialEnergy()

    def minimize_energy(self, quiet=False):
        if not quiet: print("\nInitial system energy: %s" % self.get_system_energy())
        if self.inputs.mini_nstep > 0:
            if not quiet: print("\nEnergy minimization: %s steps" % self.inputs.mini_nstep)
            self.simulation.minimizeEnergy(tolerance=self.inputs.mini_Tol*kilojoule/mole, maxIterations=self.inputs.mini_nstep)
        if not quiet: print("\nFinal system energy: %s" % self.get_system_energy())

            
    def gen_initial_velocities(self):
        if self.inputs.gen_seed:
            self.simulation.context.setVelocitiesToTemperature(self.inputs.gen_temp, self.inputs.gen_seed)
        else:
            self.simulation.context.setVelocitiesToTemperature(self.inputs.gen_temp)

    def run(self, dcd_file=None, pdb_file=None, rst_file=None, chk_file=None, quiet=False):
        if self.inputs.nstep > 0:
            if not quiet: print("\nMD run: %s steps" % self.inputs.nstep)
            if self.inputs.nstdcd > 0:
                reporter = DCDReporter(dcd_file, self.inputs.nstdcd)
                self.simulation.reporters.append(reporter)
            self.simulation.reporters.append(
                StateDataReporter(sys.stdout, self.inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
                                remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
            )
            # Simulated annealing?
            if self.inputs.annealing == 'yes':
                interval = self.inputs.interval
                temp = self.inputs.temp_init
                for i in range(self.inputs.nstep):
                    self.integrator.setTemperature(temp*kelvin)
                    self.simulation.step(1)
                    temp += interval
            else:
                self.simulation.step(self.inputs.nstep)
        else:
            pass
        if rst_file:
            state = self.simulation.context.getState( getPositions=True, getVelocities=True )
            with open(rst_file, 'w') as f:
                f.write(XmlSerializer.serialize(state))
        if chk_file:
            with open(chk_file, 'wb') as f:
                f.write(self.imulation.context.createCheckpoint())
        if pdb_file:
            crd = self.simulation.context.getState(getPositions=True).getPositions()
            PDBFile.writeFile(self.psf.topology, crd, open(pdb_file, 'w'))
        return self.simulation

    def run_with_args(inpfile, psffile, crdfile, toppar, hmr=False, rewrap=False, sysinfo=None, icrst=None, ichk, opdb, orst, ochk, odcd, ):
        self.load_params()
        self.build_system()
        self.build_integator()
        self.set_platform()
        self.build_context()
        self.rewrap()
        self.minimize_energy()
        self.gen_initial_velocities()
        self.run()
        

