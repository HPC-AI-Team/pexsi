function Ainv = PreSelInv( L, U, superPtr )
% PreSelInv performs the same procedure as PMatrix::PreSelInv in
% pselinv.cpp.
%
% This is subroutine to facilitate the debugging process.
%
% L and U are obtained from the LU factorization
%
% [L, U] = lu( A );
%
% superPtr is the same as that in SuperNode structure in PSelInv.hpp,
% and can be obtained from the output file of the test subroutine in
% pexsi.
%
% Lin Lin
% 12/9/2012

numSuper = length(superPtr) - 1;
Ainv = zeros( size( L ) );
for ksup = 1 : numSuper
  supInd = superPtr(ksup)+1 : superPtr(ksup+1);

	if( ksup < numSuper )
		% L(i,k) <- L(i,k) * L(k,k)^{-1} 
		diagL  = L(supInd, supInd);
		Ainv(superPtr(ksup+1):end, supInd) = ...
			L(superPtr(ksup+1):end, supInd) * inv( diagL );

		% U(k, i) <- L(i, k)
		Ainv(supInd, superPtr(ksup+1):end) = ...
			transpose( Ainv(superPtr(ksup+1):end, supInd) );
	end

	% L(i,i) <- [L(k,k) * U(k,k)]^{-1}
	Ainv(supInd, supInd) = ...
		inv( L(supInd, supInd) * U(supInd, supInd) );
end




