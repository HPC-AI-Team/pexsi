Tutorial              {#pageTutorial}
========

- @subpage pagePselinvComplex
- @subpage pagePEXSISolve


<!-- ************************************************************ -->
@page pagePselinvComplex Parallel selected inversion for a complex matrix
\tableofcontents

The computation of @ref defSelectedElem "selected elements" selected elements of an inverse matrix is a
standalone functionality of %PEXSI. If you are a C/C++ programmer, the
parallel selected inversion routine can be used as follows.

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
`SnzvalLocal` are double arrays.  The output array `AinvnzvalLocal` is a
double array which is twice the size of  `HnzvalLocal` due to the usage
of complex format.

An example is given in driver_pselinv.c. See also @ref
PPEXSISelInvInterface for detailed information of its usage.



<!-- ************************************************************ -->
@page pagePEXSISolve Electronic structure calculation using PEXSI

**TBD** 
