Introduction      {#pageIntro}
============

- @subpage pageOverview
- @subpage pageLicense
- @subpage pageReference
- @subpage pageChangeLog

<!-- ************************************************************ -->
@page pageOverview Overview

The Pole EXpansion and Selected Inversion method (%PEXSI) is a fast
method for evaluating certain [selected elements](@ref defSelectedElem)
of a matrix function.  %PEXSI is highly scalable on distributed memory
parallel machines. 

Given a sparse square matrix \f$A\f$ and a certain function
\f$f(\cdot)\f$, the basic idea of %PEXSI is to
expand \f$f(A)\f$ using a small number of rational functions (pole expansion) 
\f[
f(A) \approx \sum_{l=1}^{P} \omega_l(A-z_l I)^{-1},
\f]
and to efficiently evaluate \f$f(A)_{i,j}\f$ by evaluating selected
elements \f$(A-z_l I)^{-1}_{i,j}\f$ (selected inversion).

The currently supported form of \f$f(\cdot)\f$ include:

- \f$f(z)=z^{-1}\f$: Matrix inversion.  Since the matrix inversion is
  already represented as a single term of rational function (pole), no
  pole expansion is needed.  The selected inversion method can be
  significantly faster than directly inverting the matrix and then
  extract the selected elements of the inverse.
  For only using %PEXSI to evaluate selected
  elements of \f$A^{-1}\f$, see @ref pagePselinvComplex for an example.

- \f$f(z)=\frac{2}{1+e^{\beta (z-\mu)}}\f$: Fermi-Dirac function.  This can be
  used as a "smeared" matrix sign function at \f$z=\mu\f$, without
  assuming a spectral gap near \f$z=\mu\f$.  This can be used for
  evaluating the electron density for electronic structure calculation.
  See @ref pagePEXSISolve for an example using %PEXSI for electronic
  structure calculation. 

  @image html FermiDirac.png "Red: Fermi-Dirac function. Black: Matrix sign function" 
  

For sparse matrices, the %PEXSI method can be more efficient than the widely used
[diagonalization method](@ref defDiagonalization) for evaluating matrix
functions, especially when a relatively large number of eigenpairs are
needed to be computed in the diagonalization method.  
%PEXSI can also be used to compute the matrix functions associated with
generalized eigenvalue problems, i.e.

\f[
f(A,B):= V f(\Lambda) V^{-1} \approx \sum_{l=1}^{P} \omega_l(A-z_l B)^{-1},
\f]
where \f$V,\Lambda\f$ are defined through the generalized eigenvalue
problem
\f[
A V = B V \Lambda.
\f]

%PEXSI is most advantageous when a large number of processors are
available, due to the two-level parallelism.  For accurate evaluation,
the pole expansion usually takes around \f$P\approx 80\f$ poles.  All
the \f$80\f$ matrices can be inverted independently among different
groups of processors.  The parallel selected inversion method (PSelInv,
which is included in %PEXSI) can scale well to \f$256\sim 1024\f$ or
more
processors depending on the sparsity of the problem, and the system
size.  Therefore it is most advantageous to use %PEXSI when more than
1000 processors are available.  
For some problems we have also observed that it can be
advantageous to use %PEXSI using hundreds to thousands of processors.

@anchor defDiagonalization
**Diagonalization method** 

The diagonalization method evaluates a matrix function \f$f(A)\f$ by
\f[
f(A) = V f(\Lambda) V^{-1},
\f]
where the orthonormal matrix \f$V=[v_1,\ldots,v_N]\f$, and the diagonal matrix
\f$\Lambda=\mathrm{diag}(\lambda_1,\ldots,\lambda_N)\f$ are defined through the eigenvalue problem
\f[
A V = V \Lambda.
\f]
It is often the case that not all eigenvalues \f$\{\lambda_i\}\f$ are
needed to be computed, depending on the value of \f$f(\lambda_i)\f$.  

@anchor defSelectedElem 
**Selected elements** 

For structurally symmetric matrices (i.e. \f$A_{i,j}\ne 0\f$ implies
\f$A_{j,i}\ne 0\f$), we define the selected
elements of a matrix \f$B\f$ with respect to a matrix \f$A\f$ as the set
\f$\{B_{i,j}\vert A_{i,j}\ne 0\}\f$.

A commonly used case in %PEXSI is the selected elements of
\f$A^{-1}\f$, which corresponds to the set \f$\{A^{-1}_{i,j}\vert A_{i,j}\ne 0\}\f$.


<!-- ************************************************************ -->
@page pageLicense License

%PEXSI is distributed under BSD license (modified by Lawrence Berkeley
National Laboratory).

    PEXSI Copyright (c) 2012 The Regents of the University of California,
    through Lawrence Berkeley National Laboratory (subject to receipt of 
    any required approvals from U.S. Dept. of Energy).  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    (1) Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
    (2) Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
    (3) Neither the name of the University of California, Lawrence Berkeley
    National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    You are under no obligation whatsoever to provide any bug fixes, patches, or
    upgrades to the features, functionality or performance of the source code
    ("Enhancements") to anyone; however, if you choose to make your Enhancements
    available either publicly, or directly to Lawrence Berkeley National
    Laboratory, without imposing a separate written license agreement for such
    Enhancements, then you hereby grant the following license: a non-exclusive,
    royalty-free perpetual license to install, use, modify, prepare derivative
    works, incorporate into other computer software, distribute, and sublicense
    such enhancements or derivative works thereof, in binary and source code form.


<!-- ************************************************************ -->
@page pageReference References

L. Lin, A. Garcia, G. Huhs and C. Yang, SIESTA-PEXSI: Massively parallel
method for efficient and accurate ab initio materials simulation without
matrix diagonalization, [<a href="http://arxiv.org/abs/1405.0194">arXiv</a>]

M. Jacquelin, L. Lin and C. Yang, PSelInv -- A Distributed Memory
Parallel Algorithm for Selected Inversion : the Symmetric Case
[<a href="http://arxiv.org/abs/1404.0447">arXiv</a>]

L. Lin, M. Chen, C. Yang and L. He, Accelerating atomic
orbital-based electronic structure calculation via pole expansion
and elected inversion, J. Phys. Condens. Matter 25, 295501, 2013 
[<a href="http://dx.doi.org/10.1088/0953-8984/25/29/295501">journal</a>]

L. Lin, C. Yang, J. Meza, J. Lu, L. Ying and W. E, SelInv -- An
algorithm for selected inversion of a sparse symmetric matrix, ACM
Trans. Math. Software 37, 40, 2011
[<a href="http://doi.acm.org/10.1145/1916461.1916464">journal</a>]

L. Lin, C. Yang, J. Lu, L. Ying and W. E, A Fast  Parallel
algorithm for selected inversion of structured sparse matrices with
application to 2D electronic structure
calculations, SIAM J. Sci. Comput. 33, 1329, 2011 
[<a href="http://dx.doi.org/10.1137/09077432X">journal</a>]

L. Lin, J. Lu, L. Ying, R. Car and W. E, Fast algorithm for
extracting the diagonal of the inverse matrix with application to
the electronic structure analysis of metallic systems, 
Commun. Math. Sci. 7, 755, 2009
[<a href ="http://projecteuclid.org/euclid.cms/1256562822">journal</a>]

L. Lin, J. Lu, L. Ying and W. E, Pole-based approximation of the
Fermi-Dirac function, Chin. Ann. Math. 30B, 729, 2009 
[<a href="http://dx.doi.org/10.1007/s11401-009-0201-7">journal</a>]


<!-- ************************************************************ -->
@page pageChangeLog Change Log

- v0.6.0 (03/11/2014)
  - First release of %PEXSI.
  - Version integrated with the SIESTA package for Kohn-Sham density
    functional theory (KSDFT) calculation.
  - Parallel selected inversion for complex symmetric matrices.
  - Estimate the density of state profile via inertia counting.
  - Compute the density of states and local density of states.

- v0.7.0 (05/24/2014)
  - Use PPEXSIPlan to coordinate the computation, and allows the code to
    be used for C/C++/FORTRAN.
  - Templated implementation and support for both real and complex arithmetic.
  - New interface routines for FORTRAN based on ISO_C_BINDING (FORTRAN
    2003 and later).
  - Basic interface for KSDFT calculation, with a small number of input
    parameters and built-in heuristic strategies.
  - Expert interface for KSDFT calculation, providing full-control of
    the heuristics. 
  - Symbolic factorization can be reused for multiple calculations.
  - Enhanced error estimate for the pole expansion using energy as a
    guidance.

- v0.7.1 
  - Bug fix: PPEXSIPlanInitialize specifics the input according to
    mpirank instead of outputFileIndex.
