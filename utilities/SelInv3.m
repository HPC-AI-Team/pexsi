function Ainv = SelInv3( L, U, superPtr )
% SelInv3 performs the selected inversion but in another ordering
% method, that is different from SelInv2. 
%
% superPtr is the same as that in SuperNode structure in PSelInv.hpp,
% and can be obtained from the output file of the test subroutine in
% pexsi.
%
% Lin Lin
% 12/13/2012

numSuper = length(superPtr) - 1;
numCol   = length( L );
Ainv = zeros( size( L ) );
assert( numCol == superPtr( numSuper+1 ) );

% Note: numSuper-th supernode has been treated in PreSelInv.
supInd = superPtr( numSuper ) + 1 : superPtr( numSuper+1 );
Ainv( supInd, supInd ) = inv( L(supInd, supInd) * U(supInd, supInd) );

% FIXME
for ksup = numSuper-1 : -1 : 1
  supInd = superPtr(ksup)+1 : superPtr(ksup+1);
	offdiagInd = superPtr(ksup+1) + 1 : numCol;

	% LUpdateBuf(i,k) <- -\sum_j S(i,j) L(j,k)
	% This is not so efficient since it computes too many elements but let
	% us live with it first
	LBuf = L( offdiagInd, supInd ) * inv( L( supInd, supInd ) );
	LUpdateBuf = - Ainv( offdiagInd, offdiagInd ) * LBuf;
	zeroInd = find( abs(LBuf(:)) < eps );
	LUpdateBuf( zeroInd ) = 0;
	supInd(1)

  % Ainv(k,k) <- APreInv(k,k) - \sum_i transpose(L(i,k)) * LUpdateBuf(i,k)
	DBuf = inv( L(supInd, supInd) * U(supInd, supInd) ) - ...
		LBuf.' * LUpdateBuf;
	Ainv( supInd, supInd ) = (DBuf + DBuf.') / 2;
	% supInd(1)
	% supInd(end)
	% norm(DBuf-DBuf.','inf')

	% U(k,i) <- LUpdateBuf(i,k)^T
	Ainv( supInd, offdiagInd ) = transpose( LUpdateBuf );

	% L(i,k) <- LUpdateBuf(i,k)
	Ainv( offdiagInd, supInd ) = LUpdateBuf;

end




