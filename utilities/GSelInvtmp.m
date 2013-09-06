function M = GSelInvElem( G, A, L, U )
% GSelInvElem performs the element-wise generalized selected inversion
% for the matrix M = G A G
%
% The implementation follows my note "Algorithm for generalized selected
% inversion" (2012).
% 
%
% Lin Lin
% 12/20/2012

numCol   = length( A );
if( norm( size(G) - size(A), 'inf' ) > eps )
	error('The sizes of G and A do not match.');
end
if( norm( size(G) - size(L), 'inf' ) > eps )
	error('The sizes of L and A do not match.');
end

M = zeros(size(A));

%% Step 1: Calculate the intermediate W

W = full(A);
for k = 1 : numCol-1
	indnnz = find(abs(L(k+1:end,k))>eps) + k;
	W(indnnz, k) = ( W(indnnz, k) - L(indnnz, k) * W(k, k);
	W(k, indnnz) = W(indnnz, k).';
	W(indnnz, indnnz) = W(indnnz, indnnz) + ...
		L(indnnz, k) * W(k,k) * L(indnnz,k).' - ...
    W(indnnz,k) * L(indnnz, k).' - L(indnnz, k) * W(indnnz, k).';
end

W(end,end)

%% Step 2: Calculate GAG
M( numCol, numCol ) = G(numCol, numCol) * W(numCol, numCol) * ...
	G(numCol, numCol );

for k = numCol-1 : - 1 : 1
	
	indnnz = find(abs(L(k+1:end,k))>eps) + k;

	% Preparation
	tmp1 = G(indnnz, indnnz) * W(indnnz, k) / U(k,k);
	tmp2 = M(indnnz, indnnz) * L(indnnz, k);

	% L
  M(indnnz, k) = tmp1 - tmp2;

	% D
	M(k, k) = ...
		1/(U(k,k)) * W(k,k) * 1/(U(k,k)) - L(indnnz, k).' * tmp1 - ...
		tmp1.' * L(indnnz,k) + L(indnnz,k).' * tmp2;

	% U
	M(k, indnnz) = transpose( M(indnnz, k) );

end
