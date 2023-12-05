# `hamdb/examples`
This directory contains the code used to simulate monomer/metal combinations.
 - `xyzs`: initial molecular structures guessed using the Avogadro software package.
 - `funcs.py`: code that calls `pyscf` to perform geometry optimization, hessian analysis, etc.
 - `main.py`: uses the functions in `funcs.py` and organizes the results into reaction rates, energies, etc.
 - `submit.slurm`: example submission script used to run `main.py` on the MIT SuperCloud HPC cluster.
 - `example_slurm*`: example output from running the submission script.
 - `exact_env.txt`: shows the exact versions of all packages used (for reproducibility), as listed by `conda`. Don't use this to try and run the code for research purposes, since it is platform specific - you can just install `pyscf` with `pip install pyscf`.


# Usage Notes
0. The only requirement to run this code should be `pyscf`, which will also install `numpy` (the only other dependency used here). You should be able to install `pyscf` with `pip install pyscf`, or else check their documentation.
1. Make a directory called `logfiles` for `pyscf` to write its logs too - if you don't, it will just complain that the directory doesn't exist and then exit.
