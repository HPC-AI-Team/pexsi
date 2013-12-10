Core Functionality            {#pageCoreFunction}
==================

If only the driver routines are to be used, then for C, include the
interface file:

~~~~~~~~~~{.c}
#include  "c_pexsi_interface.h"
~~~~~~~~~~

For FORTRAN, there is no interface routines such as
`f_pexsi_interface.F90` yet.  However, the FORTRAN routines can directly
be used.  See @ref pageFORTRAN for more information.

The remaining section is mainly for C++ developers to have more detailed control
of the %PEXSI package.  

For C++ and usage beyond the driver routines, include the following file
~~~~~~~~~~{.cpp}
#include  "ppexsi.hpp"
~~~~~~~~~~


- @subpage pageDataType
- @subpage pagePole
- @subpage pageFactor
- @subpage pageSelInv
- @subpage pageFORTRAN

<!-- ************************************************************ -->
@page pageDataType Data type 
\tableofcontents


Basic data type    {#secBasicDataType}
===============


The basic data types `int` and `double` are redefined as `Int` and
`Real`, in order to improve compatibility for different architecture
especially on 64-bit machines (**not implemented yet**). 
The 64-bit long integer `int64_t` is also redefined as `LongInt`.

The complex arithmetic is treated using the standard C++ `<complex>`
library.  The complex data type is `std::complex<double>`, and is
redefined as `Complex` in the implementation.  


~~~~~~~~~~{.cpp}
typedef    int                   Int;
typedef    int64_t               LongInt;
typedef    double                Real;
typedef    std::complex<double>  Complex; 
~~~~~~~~~~

NumVec, NumMat, NumTns             {#secNumStructure}
======================

The design of %PEXSI tries to eliminate as much as possible the direct
usage of pointers. This helps reducing memory leak.  Commonly used
pointers are wrapped into different classes.

The most commonly used are @ref PEXSI::NumVec "NumVec", @ref
PEXSI::NumMat "NumMat", and @ref PEXSI::NumTns "NumTns", which are
wrappers for 1D array (vector), 2D array (matrix) and 3D array (tensor),
respectively.  The arrays are always saved contiguously in memory.

These wrapper classes can both actually own an array (by specifying
`owndata_=true`), and view an array (by specifying `owndata_=false`).
Vector/Matrix/Tensor elements can be accessed using `()`.  

The underlying pointer can be accessed using the member function `Data()`.

**Commonly used wrapper classes**

`NumVec`:

~~~~~~~~~~{.cpp}
typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;
~~~~~~~~~~

`NumMat`:

~~~~~~~~~~{.cpp}
typedef NumMat<bool>       BolNumMat;
typedef NumMat<Int>        IntNumMat;
typedef NumMat<Real>       DblNumMat;
typedef NumMat<Complex>    CpxNumMat;
~~~~~~~~~~

`NumTns`:

~~~~~~~~~~{.cpp}
typedef NumTns<bool>       BolNumTns;
typedef NumTns<Int>        IntNumTns;
typedef NumTns<Real>       DblNumTns;
typedef NumTns<Complex>    CpxNumTns;
~~~~~~~~~~


Compressed sparse column format    {#secCSC}
===============================

**TBD**

The class for serial CSC format matrix is @ref PEXSI::SparseMatrix
"SparseMatrix".

Distributed compressed sparse column format    {#secDistCSC}
===========================================


We use the following convention for distributed CSC format for saving a
sparse matrix.  We assume the number of processor is \f$P\f$, the number of
rows and columns of the matrix is \f$N\f$.  The class for distributed
memory CSC format matrix is @ref PEXSI::DistSparseMatrix
"DistSparseMatrix".

- `DistSparseMatrix` uses FORTRAN convention (1-based) index, and MPI
  uses C convention (0-based) index. 
- Each processor holds \f$\lfloor N/P \rfloor\f$ consequentive columns,
  with the exception that the last processor (`mpirank == P-1`) holds a
  all the remaining \f$N - (P-1) \lfloor N/P \rfloor\f$ columns. 
  The first column holds by the i-th processor is \f$ i \lfloor N/P \rfloor\f$.
  The number of columns on each local processor is usually denoted by
  `numColLocal`. 
- `colptrLocal`, which is an integer array of type `IntNumVec` of
  dimension `numColLocal + 1`, stores the pointers to the nonzero row
  indices and nonzero values in `rowptrLocal` and `nzvalLocal`,
  respectively.  
- `rowindLocal`, which is an integer array of type `IntNumVec` of
  dimension `nnzLocal`, stores the nonzero row indices in each column.
- `nzvalLocal`, which is an array of flexible type (usually `Real` or
  `Complex`) `NumVec` of dimension `nnzLocal`, stores the nonzero values
  in each column.

2D block cyclic format    {#secCyclic}
======================

**TBD**  

<!-- ************************************************************ -->
@page pagePole Pole expansion   
\tableofcontents

The pole expansion is used to expand Fermi-Dirac functions and other
derived quantities using a number of Green's functions (poles).

> @ref PEXSI::GetPoleDensity "GetPoleDensity" 

Pole expansion for the Fermi-Dirac operator.
This is the most commonly used subroutine for the pole expansion,
and can be used to compute the shifts and weights for calculating
the density matrix, the total energy, and the Hellman-Feynman
force.  This routine obtains the expansion

\f[
  f_{\beta} (z) = \frac{2}{1+e^{\beta z}} \approx 
  \mathrm{Im} \sum_{l=1}^{P} \frac{\omega^{\rho}_l}{z-z_l}
\f]


> @ref PEXSI::GetPoleDensityDrvMu "GetPoleDensityDrvMu" 

Pole expansion for the derivative of the Fermi-Dirac
operator with respect to the chemical potential mu.
This routine can be used to evaluate the derivative of the number
of electrons with respect to the chemical potential for the
Newton step for updating the chemical potential.

Note that \f$f_{\beta}\f$ does not explicitly contain \f$\mu\f$,
so this routine actually computes the expansion

\f[
   -\frac{\partial f_{\beta}}{\partial z} (z) =
   2\beta \frac{e^{\beta z}}{(1+e^{\beta z})^2} 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{\mu}_l}{z-z_l}
\f]


> @ref PEXSI::GetPoleDensityDrvMu "GetPoleDensityDrvMu" 

Pole expansion for the derivative of the Fermi-Dirac
operator with respect to the temperature T \f$(1/\beta)\f$.

This routine can be used to extrapolate the number of electrons
from a finite temperature calculation to a zero temperature
calculation, using the derivative information.  However, this
functionality is not used anymore in the current version of
%PEXSI.
                                                                
                                                                
\f[
   \frac{\partial f_{\beta}}{\partial (1/\beta)} (z) =
   2 \beta^2 z \frac{e^{\beta z}}{(1+e^{\beta z})^2} 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{T}_l}{z-z_l}
\f]



> @ref PEXSI::GetPoleHelmholtz "GetPoleHelmholtz" 

Pole expansion for the Helmholtz free energy function.

This routine can be used to compute the (Helmholtz) free energy
when finite temperature effect exists. This is especially
important for metallic system and other small gapped systems. 
This routine expands the free energy function

\f[
   f^{\mathcal{F}}_{\beta}(z) = -\frac{2}{\beta} \log 
   (1 + e^{-\beta z}) \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{\mathcal{F}}_l}{z-z_l}
\f]

> @ref PEXSI::GetPoleForce "GetPoleForce" 

This routine can be used to compute the Pulay contribution of the
atomic force in electronic structure calculations.  This term is
especially important when basis set is not complete and changes
with atomic positions.
This routine expands the free energy function

\f[
   f^{E}_{\beta}(z) = (z+\mu) f_{\beta}(z) 
   \approx \mathrm{Im} \sum_{l=1}^{P}
   \frac{\omega^{E}_l}{z-z_l}
\f]


<!-- ************************************************************ -->
@page pageFactor Factorization
\tableofcontents



Procedure for factorization       {#secProcedureFactor}
===========================

Before the selected inversion step, the matrix saved in 
[DistSparseMatrix](@ref PEXSI::DistSparseMatrix) format must first be
factorized.  In principle this can be done with any \f$LDL^T\f$
factorization or \f$LU\f$ factorization routines.  In the current
version of %PEXSI, [SuperLU_DIST
v3.3](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) is used for the
\f$LU\f$ factorization.  

@note
To avoid conflict with other routines in %PEXSI, the SuperLU_DIST
routines are encapsulated in superlu_dist_interf.cpp. Access to
SuperLU_DIST routines are made through a wrapper class
[SuperLUMatrix](@ref PEXSI::SuperLUMatrix).

The basic steps for factorization include:

- Convert a `DistSparseMatrix` into the native format (`SuperMatrix` in
  SuperLU_DIST) of the factorization routine.
  
- Symbolic factorization.

- Numerical factorization.

Related structures and subroutines
----------------------------------

> @ref PEXSI::SuperLUGrid "SuperLUGrid"

A thin interface for the mpi grid strucutre in SuperLU.

> @ref PEXSI::SuperLUOptions "SuperLUOptions"

A thin interface for passing parameters to set the SuperLU
options.

> @ref PEXSI::SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc "SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc"

Convert a distributed sparse matrix in compressed sparse
column format into the SuperLU compressed row format.  The output is
saved in the current %SuperLUMatrix.

@note
Although LU factorization is used, the routine
assumes that the matrix is strictly symmetric, and therefore the
compressed sparse row (CSR) format, used by SuperLU_DIST, gives
exactly the same matrix as formed by the compresed sparse column
format (CSC).

> @ref PEXSI::SuperLUMatrix::SymbolicFactorize "SuperLUMatrix::SymbolicFactorize"

This routine factorizes the superlu matrix symbolically.  Symbolic
factorization contains three steps.

- Permute the matrix to reduce fill-in.
- Symbolic factorize the matrix.
- Distribute the matrix into 2D block cyclic format.

This routine is controlled via 
[SuperLUOptions](@ref PEXSI::SuperLUOptions). In particular, the permutation strategy is
controlled by 
[SuperLUOptions::ColPerm](@ref PEXSI::SuperLUOptions::ColPerm).
 
> @ref PEXSI::SuperLUMatrix::NumericalFactorize "SuperLUMatrix::NumericalFactorize"

Performs LU factorization numerically. 


Example
-------


~~~~~~~~~~{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Numerical factorization
  luMat.NumericalFactorize();

  ...;
}
~~~~~~~~~~

Reuse symbolic factorization      {#secSymbolicReuse}
============================

In SuperLU_DIST, the same symbolic factorization can be reused for
factorizing different matrices.  To reuse the symbolic factorization,
one should follow the steps below.

(After symbolic factorization)
- Destroy the `SuperMatrix`.
- Convert another `DistSparseMatrix` into the native format (`SuperMatrix` in
  SuperLU_DIST) of the factorization routine.
- Redistribute the matrix into 2D block cyclic format.
- Numerical factorization.

Related structures and subroutines
----------------------------------

> @ref PEXSI::SuperLUMatrix::DestroyAOnly "SuperLUMatrix::DestroyAOnly"

Releases the data in A but keeps other data, such as LUstruct. 

This allows one to perform factorization of
matrices of the same pattern, such as the option

`fact = SamePattern_SameRowPerm`

in SuperLU_DIST.

> @ref PEXSI::SuperLUMatrix::Distribute "SuperLUMatrix::Distribute"

Distribute redistrbutes the SuperMatrix in parallel so that it is ready
for the numerical factorization.

Example
-------

~~~~~~~~~~{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Destroy the SuperMatrix saved in luMat.
  luMat.DestroyAOnly();


  // Construct another matrix BMat with the same sparsity pattern as A.
  DistSparseMatrix<Complex>  BMat; 
  ...;
  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( BMat );
  // Redistribute into 2D block cyclic format.
  luMat.Distribute();


  // Numerical factorization
  luMat.NumericalFactorize();

  ...;
}
~~~~~~~~~~

Triangular solve and accuracy check    {#secTriangularSolve}
===================================

The triangular solve routines provided by SuperLU_DIST can be used to
check the accuracy of the factorization as well as the selected
inversion.

(After numericl factorization)
- Construct the distributed right hand sides.
- Solve \f$Ax=b\f$. Multiple right hand sides can be solved
  simultaneously.

Related structures and subroutines
----------------------------------

> @ref PEXSI::SuperLUMatrix::SolveDistMultiVector "SuperLUMatrix::SolveDistMultiVector"

Solve A x = b with b overwritten by x for distributed multivector.

> @ref PEXSI::SuperLUMatrix::CheckErrorDistMultiVector "SuperLUMatrix::CheckErrorDistMultiVector"

Print out the error by direct comparison with the true solution in
distributed format.

Example
-------

The following example performs factorization, solves for a series of
right hand sides and compare the accuracy.

~~~~~~~~~~{.cpp}
#include "ppexsi.hpp"
{
  ...;
  // Construct AMat
  DistSparseMatrix<Complex>  AMat;
  ...;

  // Setup SuperLU
  SuperLUGrid g( comm, nprow, npcol );
  SuperLUOptions luOpt;
  luOpt.ColPerm = "MMD_AT_PLUS_A";
  SuperLUMatrix luMat( g );

  // Matrix conversion
  luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

  // Symbolic factorization
  luMat.SymbolicFactorize();

  // Numerical factorization
  luMat.NumericalFactorize();
  
  // Construct a global matrix (for error checking)
  SuperLUMatrix A1( g ), GA( g );
  A1.DistSparseMatrixToSuperMatrixNRloc( AMat );
  A1.ConvertNRlocToNC( GA );
  
  // Construct the distributed right hand sides and the exact solution.
  CpxNumMat xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
  CpxNumMat xTrueLocal, bLocal;
  UniformRandom( xTrueGlobal );
  GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );
  A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
  A1.DistributeGlobalMultiVector( bGlobal,     bLocal );


  // Solve and check the error.
  luMat.SolveDistMultiVector( bLocal, berr );
  luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );

  ...;
}
~~~~~~~~~~


<!-- ************************************************************ -->
@page pageSelInv Selected Inversion
\tableofcontents



Procedure for Selected Inversion     {#secProcedureSelInv}
================================


After factorizing a [SuperLUMatrix](@ref PEXSI::SuperLUMatrix) luMat (See the [factorization](@ref secProcedureFactor) page for
information on how to perform factorization), the parallel selected inversion can be computed.



@note
To provide a layer of abstraction from the matrix format used during the factorization, the [PMatrix](@ref PEXSI::PMatrix) class is used during the selected inversion.

@note
All major operations of [PMatrix](@ref PEXSI::PMatrix), including the selected inversion, are defined directly as member functions of [PMatrix](@ref PEXSI::PMatrix).

The basic steps for selected inversion are:
- Conversion from [SuperLUMatrix](@ref PEXSI::SuperLUMatrix) to [PMatrix](@ref PEXSI::PMatrix).
 
<!-- 
  Symbolic information

      SuperNodeType super; 
      PMatrix PMloc;
      luMat.SymbolicToSuperNode( super );  
  
  Numerical information, both L and U.

      luMat.LUstructToPMatrix( PMloc ); 
-->

- Preparation of communicators and preprocessing.

<!--
  Construct the communication pattern for SelInv.

      PMloc.ConstructCommunicationPattern(); 
  
  Numerical preparation so that SelInv only involves Gemm.

      PMloc.PreSelInv();  
-->

- Parallel selected inversion.

<!--
      PMloc.SelInv_P2p();
-->

- Conversion from [PMatrix](@ref PEXSI::PMatrix) back to [DistSparseMatrix](@ref PEXSI::DistSparseMatrix) format.

<!--
  Get the information in DistSparseMatrix format 

      DistSparseMatrix<Scalar> Ainv;
      PMloc.PMatrixToDistSparseMatrix( Ainv );  
-->







Related structures and subroutines
----------------------------------

> @ref PEXSI::GridType "GridType"

> @ref PEXSI::PMatrix "PMatrix"





GridType    {#secGridType}
========

PMatrix     {#secPMatrix}
=======




<!-- ************************************************************ -->
@page pageFORTRAN FORTRAN interface
\tableofcontents

All the interface routines in C has its FORTRAN version, given in
interface.cpp.  All FORTRAN interface routines start with `f_` and can
be called in FORTRAN directly.  

@note Most FORTRAN compilers mangles the subroutine names in a way that
in the corresponding C code, the subroutine name should end with an
underscore "_".  This is controlled in `make.inc`, by adding `-DAdd_` in
the `COMMONDEFS` variable.  This macro controls the behavior of
`FORTRAN()` as defined in environment.hpp, as well as the internal
routines such as `BLAS()` and `LAPACK()`.
