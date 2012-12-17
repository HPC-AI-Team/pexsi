function Ainv = SelInvRefine( A, APreInv, superPtr )
% SelInvRefine is the test program for selected inversion with
% refinement for each supernode.  For the details of the refinement see
% the OneNote
%
% "Iterative refinement of selected inversion"
%
% superPtr is the same as that in SuperNode structure in PSelInv.hpp,
% and can be obtained from the output file of the test subroutine in
% pexsi.
%
% Lin Lin
% 12/13/2012

numSuper = length(superPtr) - 1;
numCol   = length(APreInv);
Ainv = zeros( size( APreInv ) );
assert( numCol == superPtr( numSuper+1 ) );

% Note: numSuper-th supernode has been treated in PreSelInv.
supInd = superPtr( numSuper ) + 1 : superPtr( numSuper+1 );
Ainv( supInd, supInd ) = APreInv( supInd, supInd );

% FIXME
for ksup = numSuper-1 : -1 : numSuper-20
  supInd = superPtr(ksup)+1 : superPtr(ksup+1);
	offdiagInd = superPtr(ksup+1) + 1 : numCol;

	% LUpdateBuf(i,k) <- -\sum_j S(i,j) L(j,k)
	LBuf = APreInv( offdiagInd, supInd );
	LUpdateBuf = - Ainv( offdiagInd, offdiagInd ) * LBuf;
	zeroInd = find( abs(LBuf(:)) < eps );
	LUpdateBuf( zeroInd ) = 0;

  % Ainv(k,k) <- APreInv(k,k) - \sum_i transpose(L(i,k)) * LUpdateBuf(i,k)
	Ainv( supInd, supInd ) = APreInv( supInd, supInd ) -...
		transpose( LBuf ) * LUpdateBuf;

	% Iterative refinement of LUpdateBuf(i,k) and Ainv(k,k)
	eta = -A( supInd, offdiagInd ) * LUpdateBuf - ...
		A( supInd, supInd 
	
	% U(k,i) <- LUpdateBuf(i,k)^T
	Ainv( supInd, offdiagInd ) = transpose( LUpdateBuf );

	% L(i,k) <- LUpdateBuf(i,k)
	Ainv( offdiagInd, supInd ) = LUpdateBuf;

end




