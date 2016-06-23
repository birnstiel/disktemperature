# Semenov et al. 2003 opacities

This directory contains a python wrapper for the [Fortran 77 code](http://www2.mpia-hd.mpg.de/home/henning/Dust_opacities/Opacities/Code/opacity.zip) of [Semenov et al. 2003](http://dx.doi.org/10.1051/0004-6361:20031279).

Importing the module should take care of the following things:

- The source files need to reside in the directory, which should be happening when executing `make get`.
- The fortran code and fortran wrapper then need to be compiled with `make`.