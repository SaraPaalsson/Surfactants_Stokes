Introduction
==============

This is an implementation of a boundary integral solver for Stokes flow in two dimensions. It simulates surface tension driven flow for droplets and bubbles, with the option of adding also an imposed far-field flow. The method handles different viscosity ratios between the drops and the surrounding fluid, and includes the option to add (insoluble) surfactants on the droplet interfaces. 

Most of the routines are written in MATLAB. The most time consuming part is written in C and C++ and MEXed for speed. The computations are sped up using SIMD instructions and OpenMP parallellisation.


-----------------

To get started, first build the MEX-files using CMake, by 

.. code-block:: bash
	
	$ cd build
	$ cmake ..
	$ make

To run the code, execute e.g. ``template_run``, which uses parameters set up in the directory ``Indata``. 

Testing
---------
To be added in update.     
