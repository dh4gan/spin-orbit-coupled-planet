# Flux Patterns on Planets in Spin-Orbit Resonances

This code (written in C++) calculates the surface flux received by a planet in orbit of a star
with a user defined set of orbital elements and spin period. The code outputs flux as a function of time,latitude and longitude for a single body (with a user-defined Keplerian orbit).

This work originally saw publication in

Brown et al (2014), "Photosynthetic Potential of Planets in 3:2 Spin-Orbit Resonances", International Journal of Astrobiology, Volume 13, Issue 4, pp. 279-289



## Using the Code

The code was written using the Eclipse CDT plugin, and hence the Makefile to compile it must be generated using that setup.  

Once compiled, the code is executed using

`> ./spinorbit_coupled_planet`

The input parameters are read from

`spinorbit_coupled_planet.params`

Outputs can be plotted using the python scripts in `dh4gan/plot-spin-orbit-coupled-planet` 
(Developed in Python 2.7, requires numpy 1.8, matplotlib 1.3)

