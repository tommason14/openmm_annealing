from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from time import gmtime, strftime
from datetime import datetime
import os
from parmed.gromacs import GromacsGroFile, GromacsTopologyFile


def anneal(
    integrator, simulation, initial=300, final=10, jump=10, step=10000,
):
    """
    Change temperature from the initial to the final temperature, altering
    the temperature by a specified amount, jump, over a number of timesteps, step. 
    Initial, final and temperature increments are assumed to be in kelvin.
    """
    num_intervals = int(abs(final - initial) / jump)
    if final < initial:
        for i in range(num_intervals + 1):
            integrator.setTemperature((initial - (jump * i)) * kelvin)
            simulation.step(step)
    else:
        for i in range(num_intervals + 1):
            integrator.setTemperature((initial + (jump * i)) * kelvin)
            simulation.step(step)


# Thermostat parameters
temperature = 300 * kelvin
friction = 0.1 * picosecond
timestep = 0.001 * picoseconds

gro = GromacsGroFile.parse("../build/pack.gro")
top = GromacsTopologyFile("../build/topol.top")
# PBCs
top.box = gro.box

system = top.createSystem(
    rigidWater=True,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * nanometer,
    constraints=HBonds,
)

integrator = LangevinIntegrator(temperature, friction, timestep)

simulation = Simulation(top.topology, system, integrator)
with open("../npt/npt.chk", "rb") as f:  # checkpoints are binary files
    simulation.context.loadCheckpoint(f.read())

simulation.reporters = []
dcdfile = "anneal.dcd"
logfile = "anneal.csv"
chkfile = "anneal.chk"

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

for _ in range(5):
    anneal(integrator, simulation, 300, 10)
    anneal(integrator, simulation, 10, 300)

# last cooling step over a longer period of time, so spend longer at each temperature over the
# cooling range
anneal(integrator, simulation, 300, 10, step=50000)

# save final stucture
state = simulation.context.getState(getPositions=True)
positions = state.getPositions()
with open("final_structure.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)
