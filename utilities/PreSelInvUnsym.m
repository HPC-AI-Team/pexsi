function Ainv = PreSelInvUnsym( L, U, superPtr )

numSuper = length(superPtr) - 1;
Ainv = sparse(  size( L,1 ), size(L,2) );
%Ainv2 = sparse(  size( L,1 ), size(L,2) );
%for k = 1 : size(L,1)
%  Ainv2(k,k) = inv(U(k,k));
%  ind = k + find(U(k,k+1:end));  
%  if(~isempty(ind))
%    Ainv2(k,ind) = U(k,ind) * Ainv2(k,k);
%  end
%end


for ksup = 1 : numSuper
  supInd = superPtr(ksup)+1 : superPtr(ksup+1);
  for i = 1: length(supInd)
    irow = supInd(i);
    Ainv(irow,irow) = inv(U(irow,irow));
    
    ind = irow + find(U(irow,irow+1:end));  
    if(~isempty(ind))
      Ainv(irow,ind) = U(irow,ind) * Ainv(irow,irow);
    end
  end


%		diagU  = U(supInd, supInd);
%    ind = supInd(1) + find(U(supInd, superPtr(ksup+1):end));  
%    if(~isempty(ind))
%		  Ainv(supInd, ind) = ...
%			 inv(diagU) * U(supInd, ind);
%    end
%
%	  Ainv(supInd, supInd) = ...
%		  inv(L(supInd, supInd) * U(supInd, supInd) );


end
Ainv = Ainv + tril(L,-1); 
%Ainv2 = Ainv2 + tril(L,-1); 
%Ainv = Ainv2;
%Ainv(:,1)
%Ainv2(:,1)


end
