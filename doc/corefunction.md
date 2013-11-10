Core Functionality     {#pageCoreFunction}
==================

- @subpage pagePole
- @subpage pageSelInv

@page pagePole Pole expansion   
\tableofcontents

See page @ref pagePselinvComplex for more examples.

The pole expansion is used to expand Fermi-Dirac functions and other
derived quantities using a number of Green's functions (poles).

int @ref PEXSI::GetPoleDensity "GetPoleDensity" (Complex* zshift, Complex* zweight, 
    int Npole, double temp, double gap, double deltaE,
    double mu);

~~~~{.cpp}
int GetPoleDensity (Complex* zshift, Complex* zweight, 
    int Npole, double temp, double gap, double deltaE,
    double mu);
~~~~

Expand the function
\f[
  f_{\beta} (z) = \frac{2}{1+e^{\beta z}} \approx 
  \mathrm{Im} \sum_{l=1}^{P} \frac{\omega^{\rho}_l}{z-z_l}
\f]



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

