.. _pageTutorial:

Tutorial      
---------

.. _pagePEXSIPlan:

Using plans and generating log files
=====================================

PEXSI is written in C++, and the subroutines cannot directly interface
with other programming languages such as C or FORTRAN.  To solve
this problem, the PEXSI internal data structure is handled using a
datatype :ref:`PPEXSIPlan <PPEXSIPlan>`.  The idea and the usage of 
PPEXSIPlan is similar to `fftw_plan` in the
FFTW (http://www.fftw.org/~fftw/fftw3_doc/Using-Plans.html)
package.

In PEXSI, a matrix (or more accurately, its inverse) is generally
referred to as a "pole". The factorization and selected inversion
procedure for a pole is computed in parallel using `numProcRow *
numProcCol` processors.
 
When only selected inversion (PSelInv) is used, it is recommended to
set the `mpisize` of the communicator `comm` to be just `numProcRow * numProcCol`.
 
When PEXSI is used to evaluate a large number of inverse matrices
such as in the electronic structure calculation, it is best to set
`mpisize` to be `numPole*numProcRow*numProcCol`, where `numPole` is the
number of poles can be processed in parallel. 

Starting from v1.0, when `PPEXSIDFTDriver2` is used, it is best set
`mpisize` to be `numPoint*numPole*numProcRow*numProcCol`, where
`numPoint` is the number of PEXSI evaluations that can be performed in
parallel. When `mpisize < numPoint*numPole*numProcRow*numProcCol`, 
`PPEXSIDFTDriver2` will first parallel over the `numProcRow*numProcCol` 
and `numPoint`.

The output information is controlled by the `outputFileIndex` variable,
which is a local variable for each processor. For instance, if this 
index is 1, then the corresponding processor will output to the file
`logPEXSI1`.  If outputFileIndex is negative, then this processor 
does NOT output logPEXSI files.

.. note::

  - Each processor must output to a **different** file if outputFileIndex
    is non-negative.  
  - When many processors are used, it is **not recommended** for all
    processors to output the log files. This is because the IO takes time
    and can be the bottleneck on many architecture. A good practice is to
    let the master processor output information (generating `logPEXSI0`) or 
    to let the master processor of each pole to output the information. 

.. _PPEXSISelInvRealSymmetricMatrix:

Parallel selected inversion for a real symmetric matrix
========================================================


The parallel selected inversion routine for a real symmetric matrix can
be used as follows. This assumes that the mpisize of `MPI_COMM_WORLD` is
`nprow * npcol`. 

::

    #include  "c_pexsi_interface.h"
    ...
    {
      /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
      ...;
      /* Initialize PEXSI. 
       * PPEXSIPlan is a handle communicating with the C++ internal data structure */
      PPEXSIPlan   plan;
      
      plan = PPEXSIPlanInitialize( 
          MPI_COMM_WORLD, 
          nprow,
          npcol,
          mpirank, 
          &info );
    
      /* Tuning parameters of PEXSI. The default options is reasonable to
       * start, and the parameters in options can be changed.  */
      PPEXSIOptions  options;
      PPEXSISetDefaultOptions( &options );
    
      /* Load the matrix into the internal data structure */
      PPEXSILoadRealHSMatrix( 
          plan, 
          options,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          AnzvalLocal,
          1,     // S is an identity matrix here
          NULL,  // S is an identity matrix here
          &info );
    
      /* Factorize the matrix symbolically */
      PPEXSISymbolicFactorizeRealSymmetricMatrix( 
          plan,
          options,
          &info );
    
      /* Main routine for computing selected elements and save into AinvnzvalLocal */
      PPEXSISelInvRealSymmetricMatrix (
          plan,
          options,
          AnzvalLocal,
          AinvnzvalLocal,
          &info );
    
      ...;
      /* Post processing AinvnzvalLocal */
      ...; 
    
      PPEXSIPlanFinalize(
          plan,
          &info );
    } 

This routine computes the selected elements of the matrix 
:math:`A^{-1}=(H - z S)^{-1}` in parallel.  The input matrix :math:`H`
follows the :ref:`Distribute CSC format <secDistCSC>`, defined through the 
variables `colptrLocal`,`rowindLocal`, `HnzvalLocal`.  The input matrix 
:math:`S` can be omitted if it is an identity matrix and by setting 
`isSIdentity=1`. If :math:`S` is not an identity matrix, the nonzero 
sparsity pattern is assumed to be the same as the nonzero sparsity 
pattern of :math:`H`.  Both `HnzvalLocal` and `SnzvalLocal` are double 
precision arrays.  

An example is given in `examples/driver_pselinv_real.c`, which evaluates the
selected elements of the inverse of the matrix saved in
`examples/lap2dr.matrix`.  
See also :ref:`PEXSI Real Symmetric Matrix <PPEXSISelInvRealSymmetricMatrix>`
for detailed information of its usage.



.. _pagePselinvComplex:

Parallel selected inversion for a complex symmetric matrix
===========================================================


The parallel selected inversion routine for a complex symmetric matrix
is very similar to the real symmetric case. An example is given in
`examples/driver_pselinv_complex.c`. See also :ref:`PEXSI Real Symmetric
Matrix <PPEXSISelInvRealSymmetricMatrix>`
for detailed information of its usage.

.. _pagePselinvRealSymmetricUnsym:

Parallel selected inversion for a real unsymmetric matrix
==========================================================

The parallel selected inversion routine for a real unsymmetric matrix can
be used as follows. This assumes that the size of `MPI_COMM_WORLD` is
`nprow * npcol`. ::

    #include  "c_pexsi_interface.h"
    ...
    {
      /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
      ...;
      /* Initialize PEXSI. 
       * PPEXSIPlan is a handle communicating with the C++ internal data structure */
      PPEXSIPlan   plan;
      
      plan = PPEXSIPlanInitialize( 
          MPI_COMM_WORLD, 
          nprow,
          npcol,
          mpirank, 
          &info );
    
      /* Tuning parameters of PEXSI. The default options is reasonable to
       * start, and the parameters in options can be changed.  */
      PPEXSIOptions  options;
      PPEXSISetDefaultOptions( &options );
      
    
      /* Load the matrix into the internal data structure */
      PPEXSILoadRealHSMatrix( 
          plan, 
          options,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          AnzvalLocal,
          1,     // S is an identity matrix here
          NULL,  // S is an identity matrix here
          &info );
    
      /* Factorize the matrix symbolically */
      PPEXSISymbolicFactorizeRealUnsymmetricMatrix( 
          plan,
          options,
          &info );
    
      /* Main routine for computing selected elements and save into AinvnzvalLocal */
      PPEXSISelInvRealUnsymmetricMatrix (
          plan,
          options,
          AnzvalLocal,
          AinvnzvalLocal,
          &info );
    
      ...;
      /* Post processing AinvnzvalLocal */
      ...; 
    
      PPEXSIPlanFinalize(
          plan,
          &info );
    } 

This routine computes the selected elements of the matrix 
:math:`A^{-1}=(H - z S)^{-1}` in parallel.  The input matrix :math:`H`
follows the :ref:`Distribute CSC format <secDistCSC>`, defined through the variables `colptrLocal`,
`rowindLocal`, `HnzvalLocal`.  The input matrix :math:`S` can be omitted if it
is an identity matrix and by setting `isSIdentity=1`. If :math:`S` is not
an identity matrix, the nonzero sparsity pattern is assumed to be the
same as the nonzero sparsity pattern of :math:`H`.  Both `HnzvalLocal` and
`SnzvalLocal` are double precision arrays.  

.. note:: 

  As discussed in :ref:`selected elements <defSelectedElem>`,
  for general non-symmetric matrices, the selected elements are the
  elements such that :math:`\{A_{j,i}\ne 0\}`.  This means that the matrix
  elements computed corresponding to the sparsity pattern of :math:`A^T`.
  However, storing the matrix elements :math:`\{A^{-1}_{i,j}\vert
  A_{j,i}\ne 0\}` is practically cumbersome, especially in the context of
  distributed computing. Hence we choose to store the selected elements
  for :math:`A^{-T}`, i.e. :math:`\{A^{-T}_{i,j}\vert A_{i,j}\ne 0\}`.
  These are the values obtained from the non-symmetric version of PSelInv.
  


An example is given in `examples/driver_pselinv_real_unsym.c`, which evaluates the
selected elements of the inverse of the matrix saved in
`examples/big.unsym.matrix`.  See also `PPEXSISelInvRealUnsymmetricMatrix`
for detailed information of its usage.



Parallel selected inversion for a complex unsymmetric matrix
=============================================================


The parallel selected inversion routine for a complex unsymmetric matrix
is very similar to the real unsymmetric case. An example is given in
`examples/driver_pselinv_complex_unsym.c`. See also `PPEXSISelInvComplexUnsymmetricMatrix`
for detailed information of its usage.

Similar to the case of real unsymmetric matrices, the values
:math:`\{A^{-T}_{i,j}\vert A_{i,j}\ne 0\}` are the values obtained from
the non-symmetric version of PSelInv.



.. _pageDFT1:

Solving Kohn-Sham density functional theory: I
================================================


The simplest way to use PEXSI to solve Kohn-Sham density functional
theory is to use the `PPEXSIDFTDriver2` routine. This routine uses
built-in heuristics to obtain values of some parameters in PEXSI and
provides a relatively small set of adjustable parameters for users to
tune.  This routine estimates the chemical potential self-consistently
using a combined approach of inertia counting procedure and Newton's
iteration through PEXSI. Some heuristic approach is also implemented in
this routine for dynamic adjustment of the chemical potential and some
stopping criterion.

An example routine is given in `examples/driver_ksdft.c`, which solves a fake DFT
problem by taking a Hamiltonian matrix from `examples/lap2dr.matrix`.

Here is the structure of the code using the simple driver routine. ::

    #include  "c_pexsi_interface.h"
    ...
    {
      /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
      ...;
      /* Initialize PEXSI. 
       * PPEXSIPlan is a handle communicating with the C++ internal data structure */
      PPEXSIPlan   plan;
    
      /* Set the outputFileIndex to be the pole index */
      /* The first processor for each pole outputs information */ 
      if( mpirank % (nprow * npcol) == 0 ){
        outputFileIndex = mpirank / (nprow * npcol);
      }
      else{
        outputFileIndex = -1;
      }
      
      plan = PPEXSIPlanInitialize( 
          MPI_COMM_WORLD, 
          nprow,
          npcol,
          outputFileIndex, 
          &info );
    
      /* Tuning parameters of PEXSI. See PPEXSIOption for explanation of the
       * parameters */
      PPEXSIOptions  options;
      PPEXSISetDefaultOptions( &options );
    
      options.temperature  = 0.019; // 3000K
      options.muPEXSISafeGuard  = 0.2; 
      options.numElectronPEXSITolerance = 0.001;
      /* method = 1: Contour integral ; method = 2: Moussa optimized poles; default is 2*/
      options.method = 2; 
      /* typically 20-30 poles when using method = 2; 40-80 poles when method = 1 */
      options.numPole  = 20; 
      /* 2 points parallelization is set as default. */
      options.nPoints = 2; 
    
      /* Load the matrix into the internal data structure */
      PPEXSILoadRealHSMatrix( 
          plan, 
          options,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          isSIdentity,
          SnzvalLocal,
          &info );
    
      /* Call the simple DFT driver2 using PEXSI */
      PPEXSIDFTDriver2(
          plan,
          options,
          numElectronExact,
          &muPEXSI,                   
          &numElectronPEXSI,         
          &numTotalInertiaIter,   
          &info );
    
      /* Retrieve the density matrix and other quantities from the plan */
      if(mpirank < nprow * npcol ) {

      PPEXSIRetrieveRealDM(
          plan,
          DMnzvalLocal,
          &totalEnergyH,
          &info );

      PPEXSIRetrieveRealEDM(
          plan,
          options,
          EDMnzvalLocal,
          &totalEnergyS,
          &info );
      }

      /* Clean up */
      PPEXSIPlanFinalize(
          plan,
          &info );
    } 
    
.. _pageDFT2:

Solving Kohn-Sham density functional theory: II
================================================


In a DFT calculation, the information of the symbolic factorization can
be reused for different :math:`(H,S)` matrix pencil if the sparsity pattern does
not change.  An example routine is given in `examples/driver2_ksdft.c`, which solves
a fake DFT problem by taking a Hamiltonian matrix from
`examples/lap2dr.matrix`.

Here is the structure of the code using the simple driver routine. ::

    #include  "c_pexsi_interface.h"
    ...
    {
      /* Perform DFT calculation as in the previous note */
    
      /* Update and obtain another set of H and S */
    
      /* Solve the problem once again without symbolic factorization */
      PPEXSILoadRealHSMatrix( 
          plan, 
          options,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          isSIdentity,
          SnzvalLocal,
          &info );
    
      // No need to perform symbolic factorization 
      options.isSymbolicFactorize = 0;
      // Given a good guess of the chemical potential, no need to perform 
      // inertia counting.
      options.isInertiaCount = 0;
      // Optional update mu0, muMin0, muMax0 in PPEXSIOptions
    
      PPEXSIDFTDriver2(
          plan,
          options,
          numElectronExact,
          &muPEXSI,                   
          &numElectronPEXSI,         
          &numTotalInertiaIter,   
          &info );
 
    } 

.. note:: 

  The built-in heuristics in `PPEXSIDFTDriver2` may not be
  optimal. It handles only one :math:`(H,S)` pair at a time, and does
  not accept multiple matrix pairs :math:`\{(H_l,S_l)\}` as in the case of
  spin polarized calculations.  For expert users and developers, it
  should be relatively easy to dig into the driver routine, and only use
  `PEXSI::PPEXSIData::SymbolicFactorizeRealSymmetricMatrix` 
  (for symbolic factorization), 
  `PEXSI::PPEXSIData::CalculateNegativeInertiaReal` 
  (for inertia counting), and
  `PEXSI::PPEXSIData::CalculateFermiOperatorReal` 
  (for one-shot PEXSI calculation) to improve heuristics and extend the
  functionalities.



Parallel computation of the Fermi operator for complex Hermitian matrices
=======================================================================================


The PPEXSIDFTDriver routine and PPEXSIDFTDriver2 routines are standalone
routines for solving the density matrix with the capability of finding
the chemical potential effectively. This can be used for :math:`\Gamma`
point calculation. For electronic structure calculations with k-points,
multiple Hamiltonian operators may be needed to compute the number of
electrons. The PEXSI package provides expert level routines for such
purpose.  See driver_fermi_complex.c for an example of the components.
