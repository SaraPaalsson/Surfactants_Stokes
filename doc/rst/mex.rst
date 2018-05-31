Mex interface
===============

.. function:: mex_mygmresg_el(m,rhs,tol,bnrm2,@do_printout,f_handle,varargin)

   Solves (I+K)x=rhs using GMRES. 

   Returns:
     **x** -- solution vector

     **i** -- number of iterations needed

   :param m: maximum number of iterations
   :param rhs: right hand side
   :param tol: tolerance
   :param bnrm2: norm of the right hand side
   :param @do_printout: sets if printout desired or not
   :param f_handle: function to compute matrix-vector product.


.. function:: mex_M1((z,q,tol,levmax,numthreads)

    Fast Multipole code for particles in the plane designed for use in Matlab. At each particle z(k) the sum q(m)/(z(m)-z(k)) is computed over all m!=k, where q(m) are the charges of the particles.
  
    The Matlab syntax is
 
      [M1,maxlevel] = mex_M1(z,q,tol,levmax,numthreads)
 
    :param z: positions of the particles, complex
    :param q: the charges of the particles, real or complex
    :param tol: the requested accuracy(a bit pessimistic, tol = 1e-13 is sufficient for machine epsilon accuracy)
    :param levmax: the max number of box subdivisions before we resort to slow direct evaluation. This is to be able to handle (albeit slowly) highly non-uniform distributions of particles
    :param numthreads: the number of computational threads to use 
    :param M1: the output sums
    :param maxlevel: the final depth of the refinement recursion
 
    By Rikard Ojala, 2011-03-11.

.. function:: mex_M1M4(z,zp,dens,tol,levmax,numthreads)

   Fast Multipole code for particles in the plane designed for use in Matlab. At each particle z(k) (corresponding to a discretization point on a parameterized boundary in the complex plane) the sums are computed simultaneously for all m != k, where dens is the density and zp is the derivative of the boundary with respect to parameter.

  The Matlab syntax is
 
    [M1,M4,maxlevel] = mex_M1M4(z,zp,dens,tol,levmax,numthreads)

   Returns:  
     **M1,M4** -- the output sums

     **maxlevel** -- the final depth of the refinement recursion

  :param z: positions of the boundary particles, complex
    
  :param zp: tangents of the boundary, complex
    
  :param dens: density, real/complex
    
  :param tol: the requested accuracy(a bit pessimistic, tol = 1e-13 is sufficient for machine epsilon accuracy)
    
  :param levmax: the max number of box subdivisions before we resort to slow direct evaluation. This is to be able to handle (albeit slowly) highly non-uniform distributions of particles.
    
  :param numthreads: the number of computational threads to use 
  
  By Rikard Ojala, 2011-03-23.

.. function:: mex_M1M5(z,zp,dens,tol,levmax,numthreads)

   Fast Multipole code for particles in the plane designed for use in Matlab. At each particle z(k) (corresponding to a discretization point on a parameterized boundary in the complex plane) the sums are computed simultaneously for all m != k, where dens is the density and zp is the derivative of the boundary with respect to parameter.
 
   The Matlab syntax is
 
     [M1,M5,maxlevel] = mex_M1M5(z,zp,dens,tol,levmax,numthreads)
 
   Returns:
     **M1,M5** -- the output sums

     **maxlevel** -- the final depth of the refinement recursion

  :param z: positions of the boundary particles, complex
  :param zp: tangents of the boundary, complex
  :param dens: density, real/complex
  :param tol: the requested accuracy(a bit pessimistic, tol = 1e-13 is sufficient for machine epsilon accuracy)
  :param levmax: the max number of box subdivisions before we resort to slow direct evaluation. This is to be able to handle (albeit slowly) highly non-uniform distributions of particles
  :param numthreads: the number of computational threads to use 
 
  By Rikard Ojala, 2012-04-04.

.. function:: mex_Trap2GL(fmtau,fmtaup,fmtaupp,T)

   Goes from function defined on equidistant interface discretization points (on a panel) to a function on Gauss-Legendre discretization points. Called by matlab function Trap2GL 

   Returns:
     **z,zp,zpp** -- G.-L. points and their first and second derivatives

  :param fmtau: function defined on equidistant points on panel
  :param fmtaup: function first derivative
  :param fmtaupp: function second derivative
  :param T: G.-L. nodes on panel


.. function:: mex_dospecquad3((z,zp,W,pe,bubble,idx,beta,4)

   Compute special quadrature modifications, both for the integral equation if needed, and for the velocity evaluation.

   Called from MATLAB by

     [modifs,imodifs,rmodifs,cidx] = mex_dospecquad3(z,zp,W,pe,bubble,idx,beta,4);

