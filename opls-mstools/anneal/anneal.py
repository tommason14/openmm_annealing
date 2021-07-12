#!/usr/bin/env python3

import sys
import argparse
import simtk.openmm as mm
from simtk.openmm import app
import ommhelper as oh
from ommhelper.unit import *
import os

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "-n",
    "--cycles",
    help="Number of annealing cycles, default = 5",
    type=int,
    default=5,
)
parser.add_argument(
    "-tmax",
    "--cycle-max",
    help="Maximum temperature reached during annealing cycles, in kelvin. Default = 300 K",
    type=float,
    default=300,
)
parser.add_argument(
    "-tmin",
    "--cycle-min",
    help="Minimum temperature reached during annealing cycles, in kelvin. Default = 10 K",
    type=float,
    default=10,
)
parser.add_argument(
    "-tfinal",
    "--final-temp",
    help="Final temperature to cool to, in kelvin. Default = 10 K",
    type=float,
    default=10,
)
parser.add_argument("--dt", type=float, default=0.001, help="step size in ps")
parser.add_argument("--gro", type=str, default="conf.gro", help="gro file")
parser.add_argument("--psf", type=str, default="topol.psf", help="psf file")
parser.add_argument("--prm", type=str, default="ff.prm", help="prm file")
parser.add_argument("--cpt", type=str, help="load checkpoint")
args = parser.parse_args()


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


def gen_simulation(
    gro_file="conf.gro",
    psf_file="topol.psf",
    prm_file="ff.prm",
    dt=0.001,
    T=300,
    restart=None,
):
    print("Building system...")
    gro = oh.GroFile(gro_file)
    psf = oh.OplsPsfFile(psf_file, periodicBoxVectors=gro.getPeriodicBoxVectors())
    prm = app.CharmmParameterSet(prm_file)
    system = psf.createSystem(
        prm,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2 * nm,
        constraints=app.HBonds,
        rigidWater=True,
        verbose=True,
    )
    is_drude = any(type(f) == mm.DrudeForce for f in system.getForces())

    ### apply TT damping for CLPol force field
    donors = [atom.idx for atom in psf.atom_list if atom.attype == "HO"]
    if is_drude and len(donors) > 0:
        print("Add TT damping between HO and Drude dipoles")
        ttforce = oh.CLPolCoulTT(system, donors)
        print(ttforce.getEnergyFunction())

    print("Initializing simulation...")
    if is_drude:
        print("Drude Langevin thermostat: 5.0 /ps, 20 /ps")
        integrator = mm.DrudeLangevinIntegrator(
            T * kelvin, 5.0 / ps, 1 * kelvin, 20 / ps, dt * ps
        )
        integrator.setMaxDrudeDistance(0.02 * nm)
    else:
        print("Langevin thermostat: 1.0 /ps")
        integrator = mm.LangevinIntegrator(T * kelvin, 1.0 / ps, dt * ps)

    _platform = mm.Platform.getPlatformByName("CUDA")
    _properties = {"CudaPrecision": "mixed"}
    sim = app.Simulation(psf.topology, system, integrator, _platform, _properties)
    if restart:
        sim.loadCheckpoint(restart)
        sim.currentStep = (
            round(sim.context.getState().getTime().value_in_unit(ps) / dt / 10) * 10
        )
        sim.context.setTime(sim.currentStep * dt)
    else:
        sim.context.setPositions(gro.positions)
        sim.context.setVelocitiesToTemperature(T * kelvin)

    # For each reporter, check if we should append to the file - allows for restart files to
    # be used in another directory
    append_dcd = "dump.dcd" in os.listdir(".")
    sim.reporters.append(
        app.DCDReporter("dump.dcd", 10000, enforcePeriodicBox=False, append=append_dcd)
    )
    sim.reporters.append(oh.CheckpointReporter("cpt.cpt", 10000))

    append_gro = "dump.gro" in os.listdir(".")
    sim.reporters.append(
        oh.GroReporter("dump.gro", 1000, logarithm=True, append=append_gro)
    )
    # if appending to dumps, also append here
    sim.reporters.append(
        oh.StateDataReporter(
            sys.stdout, 1000, box=False, volume=True, append=append_dcd
        )
    )
    if is_drude:
        append_drude = "T_drude.txt" in os.listdir(".")
        sim.reporters.append(
            oh.DrudeTemperatureReporter("T_drude.txt", 10000, append=append_drude)
        )
    return integrator, sim


if __name__ == "__main__":
    oh.print_omm_info()
    integrator, simulation = gen_simulation(
        gro_file=args.gro,
        psf_file=args.psf,
        prm_file=args.prm,
        dt=args.dt,
        T=args.cycle_max,
        restart=args.cpt,
    )

    print(
        f"Oscillating temperature between {args.cycle_max} K and {args.cycle_min} K over {args.cycles} cycles"
    )

    for _ in range(args.cycles):
        anneal(integrator, simulation, args.cycle_max, args.cycle_min)
        anneal(integrator, simulation, args.cycle_min, args.cycle_max)

    # last cooling step over a longer period of time, so spend longer at each temperature over the
    # cooling range
    anneal(integrator, simulation, args.cycle_max, args.final_temp, step=50000)

    # save final stucture
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    with open("final_structure.pdb", "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions, f)
