For Selected Inversion Only         {#PageSelInv}
===========================
\tableofcontents




GridType    {#SecGridType}
========

PMatrix     {#SecPMatrix}
=======

The class PMatrix contains the following subroutines in PEXSI.

@ref PEXSI::PMatrix structure

PEXSI::NumMat structure

PEXSI::PMatrix::ConstructCommunicationPattern\_Bcast

`void PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );`

void PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );

void PEXSI::PMatrix::PMatrixToDistSparseMatrix2( const DistSparseMatrix<Scalar>& A, DistSparseMatrix<Scalar>& B );

~~~~~~~~~~{.cpp}
void PEXSI::PMatrix::PMatrixToDistSparseMatrix2( const PEXSI::DistSparseMatrix<Scalar>& A, PEXSI::DistSparseMatrix<Scalar>& B );
~~~~~~~~~~
    Convert a matrix    

