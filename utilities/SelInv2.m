function Ainv = SelInv2( L, U, superPtr )
% SelInv2 performs the selected inversion but in another ordering
% method.  Hopefully, this is more stable. This turns out to be a
% HORRIBLE scheme.
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
for ksup = numSuper-1 : -1 : numSuper-30
  supInd = superPtr(ksup)+1 : superPtr(ksup+1);
	offdiagInd = superPtr(ksup+1) + 1 : numCol;

	% LUpdateBuf(i,k) <- -\sum_j S(i,j) L(j,k)
	% This is not so efficient since it computes too many elements but let
	% us live with it first
	LBuf = L( offdiagInd, supInd );
	LUpdateBuf = - Ainv( offdiagInd, offdiagInd ) * LBuf;
	zeroInd = find( abs(LBuf(:)) < eps );
	LUpdateBuf( zeroInd ) = 0;

  % Ainv(k,k) <- APreInv(k,k) - \sum_i transpose(L(i,k)) * LUpdateBuf(i,k)
	if(0)
		DBuf = inv(diag(diag( U(supInd, supInd) ))) - ...
			transpose( LBuf ) * LUpdateBuf;

		Lkk = L( supInd, supInd );

		Ainv( supInd, supInd ) = (transpose(Lkk) \ DBuf) / Lkk;

		LUpdateBuf = LUpdateBuf / Lkk;
	end

	if(0)
		DBuf = inv(diag(diag( U(supInd, supInd) ))) - ...
			transpose( LBuf ) * LUpdateBuf;

		invLkk = inv( L( supInd, supInd ) );

		Ainv( supInd, supInd ) = transpose(invLkk) * (DBuf) * invLkk;

		LUpdateBuf = LUpdateBuf * invLkk;
	end

	% U(k,i) <- LUpdateBuf(i,k)^T
	Ainv( supInd, offdiagInd ) = transpose( LUpdateBuf );

	% L(i,k) <- LUpdateBuf(i,k)
	Ainv( offdiagInd, supInd ) = LUpdateBuf;

end




