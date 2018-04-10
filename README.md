# twostars
Analytic Insolation Calculator for Circumbinary Planets

----------------------------------------------
Written by:

Siegfried Eggl  2015 03 03

----------------------------------------------
Last modified:

2018 04 10

----------------------------------------------
Language:

Fortran 90

----------------------------------------------
References:

Georgakarakos, N., & Eggl, S. (2015). Analytic orbit propagation for transiting circumbinary planets. The Astrophysical Journal, 802(2), 94.

Popp, M., & Eggl, S. (2017). Climate variations on Earth-like circumbinary planets. Nature communications, 8, 14957.

----------------------------------------------
Background:

Circumbinary planets experience potentially large changes in the amount of incident starlight (insolation). This is due to the constantly changing distances between the stars and the planet.
The module at hand incorporates an analytic orbit propagator (Georgakarakos & Eggl 2015) that describes the orbit evolution of circumbinary planets as well as routines to evaluate the amount and spectral distribution of the insolation on an extended rotating planet. Stellar eclipses are included.
In Popp & Eggl (2017) this module has been coupled with a (proprietary) general circulation model (GCM) to investigate the impact of variable insolation conditions on the climate of Earth-like circumbinary planets.


----------------------------------------------
Documentation:

A detailed documentation of the module and its subroutines is provided in the wiki of this repository.


----------------------------------------------
License: 

GNU General Public License v3.0