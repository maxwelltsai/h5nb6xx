# What is `H5nb6xx`

Admittedly, `H5nb6xx` is a strange name combining `HDF5` (https://en.wikipedia.org/wiki/Hierarchical_Data_Format) and `NBODY6++` (https://github.com/nbodyx/Nbody6ppGPU/). It is an `AMUSE` (https://github.com/amusecode/amuse) gravitational dynamics interface, which instead of integrating an N-body system, it reads the precalculated HDF5 data output of `NBODY6++` and returns the dynamical particle data.

# Why `H5nb6xx`

The development of this interface is originally motivated to couple the star cluster dynamics and planetary system dynamics. Planetary systems and star clusters opearate in very different dynamical timescales, which leads to a hierarchical problem. `H5nb6xx` allows the user to precalculate the star cluster dyanmics, and then use the recorded stellar data to perturb planetary systems. This also makes the star cluster data exactly reproducible in subsequent simulations.

# Installation
Clone this project to the `/src/amuse/community/` directory of the `AMUSE` code, and compile it.

# Publications

- Astrophysical paper: http://adsabs.harvard.edu/abs/2017arXiv170603789C

- Technical paper: http://adsabs.harvard.edu/abs/2015ApJS..219...31C
