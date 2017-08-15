function [L,U]=LU_factor(A)
% LU factorization of an n by n matrix A
% using Gauss elimination without pivoting
% LU_factor.m
% A is factored as A = L*U
% Output:
% L is lower triangular with the main diagonal part = 1s.
% U is upper triangular and is stored in the original mtx A
% and must be zeroed out to get U
% K. Ming Leung, 01/26/03

n = size(A,2);
L=speye(n);
for k=1:n
  if (A(k,k) == 0) Error('Pivoting is needed!'); end
  L(k+1:n,k)=A(k+1:n,k)/A(k,k);
  for j=k+1:n
    A(j,:)=A(j,:)-L(j,k)*A(k,:);
  end
end
U=triu(A);
