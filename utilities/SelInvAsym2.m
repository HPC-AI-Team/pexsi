function Ainv = SelInvAsym2( A )
% Selected inversion routine for asymmetric matrix. 
%
% This is test implementation 2.
%
% Lin Lin
% 02/14/2015

[L,U,p,q] = lu(A,'vector');

B = A(p,q);
% Binv1 = inv(B);

numCol   = length( A );
Binv = sparse( numCol, numCol );
% Symbolic phase.  Figure out the fill in pattern.
% Note: this is not an efficient implementation
symbBinv = sparse( numCol, numCol );
for k = 1 : numCol
	indnnzL = find(abs(L(k:end,k))>eps) + (k-1);
  indnnzU = find(abs(U(k,k:end))>eps) + (k-1);
  indnnz = union(indnnzL, indnnzU);
  symbBinv(indnnz, indnnz) = 1;
end

Binv( numCol, numCol ) = 1.0 / U(numCol,numCol);
for k = numCol-1 : - 1 : 1
	
  indnnz =  find(symbBinv(k+1:end,k)) + k;
%   indnnz =  (k+1):numCol;
  
  Lhat = L(indnnz,k) / L(k,k);
  Uhat = U(k,indnnz) / U(k,k);

	% L
	Binv(indnnz, k) = -Binv(indnnz, indnnz) * Lhat;

	% D
	Binv(k, k) = 1.0 / (L(k,k) * U(k, k)) - Uhat * Binv(indnnz, k);

	% U
	Binv(k, indnnz) = - Uhat * Binv(indnnz, indnnz);
  
end


Ainv2 = sparse(numCol, numCol);
% Note the permutation of the column and row
Ainv2(q,p) = Binv;
Ainv = sparse(numCol, numCol);
% Note compare the nonzero pattern of A'
indA = find(A');
Ainv(indA) = Ainv2(indA);

disp('Checking the accuracy..')
Ainv1=inv(A);norm(full(Ainv(indA)-Ainv1(indA)))