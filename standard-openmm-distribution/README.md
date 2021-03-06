# Using built-in tools 

Simulated annealing in the NVT ensemble of a box of 500 TIP4P water molecules.

General workflow would be to equilibrate the simulation in the NPT ensemble prior to annealing.

# Structure generation

In this scenario, the Gromacs `gmx insert-molecules` program was used to generate the initial 
configuration, but any reasonable tool can be used. The input scripts can be modified so that
Gromacs, CHARMM or Amber file formats can be read into OpenMM.
Note that the gromacs parsers built into OpenMM do not recognise virtual sites, so the parmed 
gromacs readers were used instead.

# Annealing

An `anneal` function has been included in the [input file](anneal/anneal.py) so that looping
can be used to generate temperature fluctuations. For example:

```
for _ in range(5):
    anneal(integrator, simulation, 300, 10)
    anneal(integrator, simulation, 10, 300)
```

will lower the temperature from 300 K to 10 K, and then raise it again to 300 K over 5 cycles.
This can be observed when plotting the data in
[anneal/anneal.csv](anneal/anneal.csv):
![](images/temp.png)

Here is the effect on potential energy:
![](images/pot.png)
