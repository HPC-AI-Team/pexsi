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

  

<!-- ************************************************************ -->
@page pagePole Pole expansion   
\tableofcontents

The pole expansion is used to expand Fermi-Dirac functions and other
derived quantities using a number of Green's functions (poles).

@ref PEXSI::GetPoleDensity "GetPoleDensity" 

Pole expansion for the Fermi-Dirac operator.
This is the most commonly used subroutine for the pole expansion,
and can be used to compute the shifts and weights for calculating
the density matrix, the total energy, and the Hellman-Feynman
force.  This routine obtains the expansion

\f[
  f_{\beta} (z) = \frac{2}{1+e^{\beta z}} \approx 
  \mathrm{Im} \sum_{l=1}^{P} \frac{\omega^{\rho}_l}{z-z_l}
\f]


@ref PEXSI::GetPoleDensityDrvMu "GetPoleDensityDrvMu" 

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


@ref PEXSI::GetPoleDensityDrvMu "GetPoleDensityDrvMu" 

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



@ref PEXSI::GetPoleHelmholtz "GetPoleHelmholtz" 

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

@ref PEXSI::GetPoleForce "GetPoleForce" 

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
@page pageSelInv Selected Inversion
\tableofcontents




GridType    {#secGridType}
========

PMatrix     {#secPMatrix}
=======

The class PMatrix contains the following subroutines in PEXSI.

@ref PEXSI::PMatrix 

@ref PEXSI::PMatrix "PMatrix"

void @ref PEXSI::PMatrix::PMatrixToDistSparseMatrix2 "PMatrixToDistSparseMatrix2"
( const @ref PEXSI::DistSparseMatrix "DistSparseMatrix"<Scalar>& A, @ref PEXSI::DistSparseMatrix "DistSparseMatrix"<Scalar>& B )

PEXSI::NumMat structure

PEXSI::PMatrix::ConstructCommunicationPattern\_Bcast

`void PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );`

void PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );

void PEXSI::PMatrix::PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );

~~~~~~~~~~{.cpp}
void PEXSI::PMatrix::PMatrixToDistSparseMatrix2( const PEXSI::DistSparseMatrix<Scalar>& A, PEXSI::DistSparseMatrix<Scalar>& B );
~~~~~~~~~~
    Convert a matrix    

<!-- ************************************************************ -->
@page pageFORTRAN FORTRAN interface
\tableofcontents
