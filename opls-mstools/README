# Structure generation

Creation of simulation files is handled via the [mstools package](https://github.com/z-gong/ms-tools),
so ensure this is installed before running `create_openmm.py`. The input files also use the ommhelper module in the same package.

On MonARCH, add the following line to your `~/.bashrc` and then source the file again:
```
export PYTHONPATH=$PYTHONPATH:~/p2015120004/apps/openmm-opls:~/p2015120004/apps/openmm-opls/ms-tools
```

Many types of labelled geometry files cna be used including xyz, zmat or pdb formats. Forcefield files must be accepted by mstools, so an [fftool](https://github.com/paduagroup/fftool)-compatible format is recommended. If virtual sites are
required (i.e. for TIP4P or SWM4-NDP water models), the zfp style must be used as the fftool format doesn't support virtual sites.
An fftool forcefield can be converted to a zfp style with the ffconv script distributed with the mstools package, to allow for combination with other files that include virtual site terms. For example: `ffconv.py oplsaa.ff alpha.ff -o oplsaa.zfp`. See the fftool-to-zfp folder for example files.

To create a box of water with toluene, adding virtual sites (`-v`) and drude particles (`-d`), use the following command:
```
create_openmm.py -f toluene.zmat TIP3P.zmat -n 6 600 -ff toluene.zfp SWM4-NDP.zfp -b 30 -d -v
```

Charges are assigned from the forcefield files after the extra particles have been added. 3 files are generated:

- conf.gro, containing residue information and coordinates
- topol.psf, containing topology information
- ff.prm, containing parameters for all bonded and non-bonded terms

These files are used in the subsequent simulations.

# NpT equilibration

After creating a system with an approximate box size, energy minimisation and NpT equilibration can be performed by running [npt/run-bulk.py](npt/run-bulk.py). 
An example SLURM script for the MonARCH cluster has been provided. Note that GPUs are assumed to be available as the script sets the computing platform to 'CUDA'.

# Annealing

After equilibration, annealing is performed using [anneal/anneal.py](anneal/anneal.py) with the checkpoint file from the output of the NpT equilibration. Passing in the `tmax`, `tmin` and `n`/`cycles` parameters tells OpenMM to oscillate temperature between `tmax` and `tmin` for `n` cycles, and then setting the `tfinal` parameter allows the system to cool to the desired final temperature. Then a pdb file named `final_structure.pdb` is written.

This pdb contains all drude particles and virtual sites. To remove these, use `scripts/remove_additional_particles.tcl` with the vmd program like so:
```
vmd -dispdev text -e scripts/remove_additional_particles.tcl -pdb final_structure.pdb
```
This creates a file named `cleaned.pdb` that contains only real atoms.

# Plotting the output

If you are running polarisable simulations, the temperature is saved to T_drude.txt, with separate columns for the atom and drude particle temperatures. To plot the atom temperatures, run `scripts/openmm_plot.sh T_drude.txt` and select 'T_COM', or include the option by running `echo T_COM | openmm_plot.sh T_drude.txt`. The same script will work with the log files containing potential energies and densities.
