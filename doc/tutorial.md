Tutorial              {#pageTutorial}
========

- @subpage pagePselinvComplex
- @subpage pagePEXSISolve


<!-- ************************************************************ -->
@page pagePselinvComplex Parallel selected inversion for a complex matrix
\tableofcontents

The computation of @ref defSelectedElem "selected elements" of an
inverse matrix is a standalone functionality of %PEXSI. For 
C/C++ programmers, the parallel selected inversion routine can be used as
follows.

~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
  ...;
  /* Main routine for computing selected elements and save into AinvnzvalLocal */
  PPEXSISelInvInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      zShift,
      ordering,
      npSymbFact,
      MPI_COMM_WORLD,
      AinvnzvalLocal,
      &info );
  ...;
  /* Post processing AinvnzvalLocal */
  ...; 
} 
~~~~~~~~~~ 

This routine computes the selected elements of the matrix 
\f$A^{-1}=(H - z S)^{-1}\f$ in parallel.  The input matrix \f$H\f$
follows the @ref secDistCSC, defined through the variables `colptrLocal`,
`rowindLocal`, `HnzvalLocal`.  The input matrix \f$S\f$ can be omitted if it
is an identity matrix and by setting `isSIdentity=1`. If \f$S\f$ is not
an identity matrix, the nonzero sparsity pattern is assumed to be the
same as the nonzero sparsity pattern of \f$H\f$.  Both `HnzvalLocal` and
`SnzvalLocal` are double precision arrays.  The output array `AinvnzvalLocal` is a
double array which is twice the size of  `HnzvalLocal` due to the usage
of complex format.

An example is given in driver_pselinv.c. See also @ref
PPEXSISelInvInterface for detailed information of its usage.



<!-- ************************************************************ -->
@page pagePEXSISolve Solving Kohn-Sham density functional theory

%PEXSI provides the following interface routines for electronic
structure calculation based on the Kohn-Sham density functional theory.
An example routine is given in driver_ksdft.c.

1) Estimate the range of chemical potential
-------------------------------------------

Starting from a rough estimate of the chemical potential \f$\mu\f$,
[muMin0, muMax0], obtain a refined interval for the chemical potential
[muMinInertia, muMaxInertia].

~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Setup the input matrix in distributed compressed sparse column (CSC) format */ 
  ...;
  /* Step 1. Estimate the range of chemical potential */
  PPEXSIInertiaCountInterface(
    nrows,  
    nnz,    
    nnzLocal, 
    numColLocal,
    colptrLocal,
    rowindLocal,
    HnzvalLocal,
    isSIdentity,
    SnzvalLocal,
    temperature,
    numElectronExact,
    muMin0,          
    muMax0,         
    numPole,       
    maxIter,      
    numElectronTolerance,
    ordering,           
    npPerPole,         
    npSymbFact,       
    comm,            
    &muMinInertia,                
    &muMaxInertia,               
    &muLowerEdge,               
    &muUpperEdge,              
    &numIter,                 
    &muList,                 
    &numElectronList,       
    &info                  
    );
  ...; 
} 
~~~~~~~~~~ 

See @ref PPEXSIInertiaCountInterface for detailed information of its usage.


@note The main purpose of the step of estimating the range of chemical
potential is to allow a wide range of initial guess of the chemical
potential.  When the guess of the chemical potential is already accurate,
such as the case in consecutive steps of the self-consistent-field (SCF)
iteration, the Newton iteration used in the next step will be efficient
enough.  In such case the step of estimating the range of chemical
potential can be skipped. 


2) Solve electronic structure problem
-------------------------------------

After obtaining the range of the chemical potential and a reasonably
accurate guess for the chemical potential.  The electronic structure
problem for obtaining the selected elements of the Fermi operator
(a.k.a. the single particle density matrix) can be
computed using %PEXSI.  The chemical potential is refined at the same
time using Newton's method.  Besides the Fermi operator, the energy
density matrix (for computing the force) and the free energy density
matrix (for computing the Helmholtz free energy) can be computed at the
same time.

@note The iteration of the chemical potential is usually in the inner
loop of an electronic structure calculation.  The self-consistent field
(SCF) iteration serves as an outer loop.  Because it usually takes
a few SCF steps to converge the electron density, it is not necessary to
achieve very high accuracy in the inner loop, especially during the
initial few SCF steps.

~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Step 2. Solve KSDFT using PEXSI */
  PPEXSISolveInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      temperature,
      numElectronExact,
      muInertia,
      muMinInertia,
      muMaxInertia,
      gap,
      deltaE,
      numPole,
      muMaxIter,
      PEXSINumElectronTolerance,
      ordering,
      npPerPole,
      npSymbFact,
      MPI_COMM_WORLD,
      DMnzvalLocal,
      EDMnzvalLocal,
      FDMnzvalLocal,
      &muPEXSI,
      &numElectron,
      &muMinPEXSI,
      &muMaxPEXSI,
      &muIter,
      muList,
      numElectronList,
      numElectronDrvList,
      &info );
  ...; 
} 
~~~~~~~~~~ 


See @ref PPEXSISolveInterface for detailed information of its usage.


3) Post processing (optional)
-----------------------------

After the chemical potential has converged, post-processing steps can be
performed.  Currently %PEXSI supports the computation of the density of
states (zero temperature) via counting the negative inertia.  

~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Step 3. Post processing */

  /* Compute the density of states (DOS) via inertia counting (without
   * including finite temperature effects) */

  PPEXSISolveInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      temperature,
      numElectronExact,
      muInertia,
      muMinInertia,
      muMaxInertia,
      gap,
      deltaE,
      numPole,
      muMaxIter,
      PEXSINumElectronTolerance,
      ordering,
      npPerPole,
      npSymbFact,
      MPI_COMM_WORLD,
      DMnzvalLocal,
      EDMnzvalLocal,
      FDMnzvalLocal,
      &muPEXSI,
      &numElectron,
      &muMinPEXSI,
      &muMaxPEXSI,
      &muIter,
      muList,
      numElectronList,
      numElectronDrvList,
      &info );
  ...; 
} 
~~~~~~~~~~ 


See @ref PPEXSIRawInertiaCountInterface for detailed information of its usage.

Besides the DOS, the spatially resolved local DOS can also be computed
without computing eigenvalues or eigenvectors.


~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
...
{
  /* Step 3. Post processing */

  /* Compute the local density of states (LDOS).
   *
   * Only the first pole group participates in the computation of the
   * selected inversion for a single shift. */

    PPEXSILocalDOSInterface(
        nrows,
        nnz,
        nnzLocal,
        numColLocal,
        colptrLocal,
        rowindLocal,
        HnzvalLocal,
        isSIdentity,
        SnzvalLocal,
        Energy,
        eta,
        ordering,
        npSymbFact,
        readComm,
        localDOSnzvalLocal,
        &info);

  ...; 
} 
~~~~~~~~~~ 


See @ref PPEXSILocalDOSInterface for detailed information of its usage.

@note When `mpisize > npPerPole`, @ref PPEXSILocalDOSInterface should
not be executed by all processors, but only by a subgroup of processors.
For more information see driver_ksdft.c.
