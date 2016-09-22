% Generate a test Hermitian matrix, which is the discretized version of
%
% -(\nabla + i k)^2 = -\Delta - 2ik \nabla + k^2
%
% This is implemented using central difference formula for Delta and
% nabla for Hermitian

N = 20;
k = [0.5;0];
e = ones(N,1);
nab1D = 0.5*spdiags([-e e],[-1 1],N,N);
Lap1D = spdiags([e -2*e e],[-1 0 1],N,N);

A = -kron(speye(N),Lap1D)-kron(Lap1D,speye(N)) ...
  -2i*k(1)*kron(speye(N),nab1D)-2i*k(2)*kron(nab1D,speye(N)) ...
  +(k(1)^2+k(2)^2)*speye(N^2);
    
% Input/output of the formatted matrix in complex arithmetic is not
% supported. Therefore output the real and imaginary part separately and
% do the same thing when reading the matrix.

[colptr, rowind, nzvalr] = ...
  WriteSparseMatrixLU(real(A), 'lap2dc_real.matrix',1);

nzvali = ...
  WriteSparseMatrixLU2(imag(A), 'lap2dc_imag.matrix', colptr, rowind, 1);
