from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import os

from parmed.gromacs import GromacsGroFile, GromacsTopologyFile

# Thermostat parameters
temperature = 300 * kelvin
pressure = 1.0 * atmosphere
barofreq = 100
friction = 0.1 * picosecond
timestep = 0.001 * picoseconds

gro = GromacsGroFile.parse("../build/pack.gro")
top = GromacsTopologyFile("../build/topol.top")
# PBCs
top.box = gro.box

barostat = MonteCarloBarostat(pressure, temperature)
system = top.createSystem(
    rigidWater=True,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * nanometer,
    constraints=HBonds,
)
system.addForce(barostat)

integrator = LangevinIntegrator(temperature, friction, timestep)

simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
print("Minimizing...")
simulation.minimizeEnergy()

state = simulation.context.getState(
    getEnergy=True, getForces=True, getVelocities=True, getPositions=True
)
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))

position = state.getPositions()
PDBFile.writeFile(simulation.topology, position, open("min.pdb", "w"))

print("Equilibrating...")
simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters = []
dcdfile = "npt.dcd"
logfile = "npt.csv"
chkfile = "npt.chk"

simulation.reporters.append(DCDReporter(dcdfile, 10000))
simulation.reporters.append(CheckpointReporter(chkfile, 10000))
simulation.reporters.append(
    StateDataReporter(
        logfile,
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        speed=True,
    )
)
print("Simulating...")
simulation.step(500000)
