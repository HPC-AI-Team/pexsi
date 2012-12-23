function Ainv = SelInvElem( L, U )
% SelInvElem performs the element-wise selected inversion based on the
% LU-factorization.
%
% This is subroutine to facilitate the debugging process.
%
% Lin Lin
% 12/13/2012

numCol   = length( L );
Ainv = zeros( size( L ) );

Ainv( numCol, numCol ) = 1 / U(numCol,numCol);
for k = numCol-1 : - 1 : 1
	
	indnnz = find(abs(L(k+1:end,k))>eps) + k;

	% L
	Ainv(indnnz, k) = -Ainv(indnnz, indnnz) * L(indnnz,k);

	% D
	Ainv(k, k) = 1 / U(k, k) - transpose(L(indnnz,k)) * ...
		Ainv(indnnz, k);

	% U
	Ainv(k, k+1:end) = transpose( Ainv(k+1:end, k) );

end
