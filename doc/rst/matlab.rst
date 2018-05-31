Matlab functions
=================
A collection of the MATLAB functions used to simulated surfactant-covered droplets in two dimensional Stokes flow. Code still under development.


Main and initialization
------------------------
These functions initialize the domain, set up the simulation and timestep it forward. 

.. automodule:: main 
   :members:

Stokes
------
Here the velocity given by the linear Stokes equations is computed. We use Gauss-Legendre quadrature, but the methods can use others with some modifications. Matrix-vector multiplications are computed through FMM for speed, and special quadrature is used for highly accurate results. 

.. automodule:: stokes
   :members:

Surfactants
------------
We solve the convection-diffusion equation for insoluble surfactants through a spectral method. Also some data needed for the surfactant computations especially are computed here. 

.. automodule:: surfact
   :members:

Special Quadrature
--------------------
Function going between equidistant discretization points and a panel-based Gauss-Legendre discretization can be found here. The main function for computing the special quadrature can be found among the MEX-files. 

.. automodule:: specquad
   :members:

