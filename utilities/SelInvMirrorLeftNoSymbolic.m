function Ainv = SelInvMirrorLeftNoSymbolic( A, chkerr )
% Prototype implementation of the "Mirror left-looking" version of
% selected inversion. 
%
%   Ainv = SelInvMirrorLeft( A, chkerr )
%
% Assume A is pre-ordered, and the L, U factors exist. It first computes
% the column-wise factorization of A, and then compute selected
% inversion restricted to the pattern of (L+U)+(L+U)'.
%
% The "Mirror left-looking", i.e. right-looking selected inversion
% algorithm is implemented. The factorization is right-looking.
%
% This is a column-by-column implementation and is very inefficient.
%
% Lin Lin
% 11/23/2015

numCol   = length( A );
assert( size(A,1) == size(A,2) );
if( nargin < 2 )
  chkerr = 0;
end

%% Symbolic LU factorization
%  This is a very simple implementation.
%  The resulting LU factor is structurally symmetric.
tic;

indA = find(A);
% Symbolic copy of A+A', together with all diagonals
symA = sparse( numCol, numCol ); 
symA(indA) = eps;
symA = symA + symA';

% Allocate the space for L,U  
% Lower triangular: L with unit diagonal.
% Upper triangular: U
Afactor = symA;
Afactor(indA) = A(indA);
timeSymFac = toc;
fprintf('Time for symbolic factorization = %15.2e\n', timeSymFac);

%% Numerical LU factorization
tic;
for k = 1 : numCol
  ind = k+find(Afactor(k+1:end,k));
  if(~isempty(ind))
    Afactor(ind,k) = Afactor(ind,k) / Afactor(k,k);
    Afactor(k,ind) = Afactor(k,ind);
    Afactor(ind,ind) = Afactor(ind,ind) - Afactor(ind,k) * Afactor(k,ind);
  end
end
timeNumFac = toc;
fprintf('Time for numerical factorization = %15.2e\n', timeNumFac);

if(chkerr)
  L = speye(numCol)+tril(Afactor,-1);
  U = triu(Afactor);
  fprintf('norm(A-L*U) = %15.5e\n', normest(A-L*U));
end

%% Selected inversion. Mirror left-looking
tic;
Ainv = Afactor;
% Pre-selected inversion
for k = 1 : numCol
  Ainv(k,k) = inv(Ainv(k,k));
	ind = k + find(symA(k+1:end,k));  
  if(~isempty(ind))
    Ainv(k,ind) = Ainv(k,ind) * Ainv(k,k);
  end
end

for k = numCol : -1 : 1
	ind = k + find(symA(k+1:end,k));  
  if(~isempty(ind))
    Ainv(ind,k) = -Ainv(ind,ind) * Ainv(ind,k);
    Ainv(k,k) = Ainv(k,k) - Ainv(k,ind) * Ainv(ind,k);
    Ainv(k,ind) = -Ainv(k,ind) * Ainv(ind,ind);
  end
end
timeSelInv = toc;
fprintf('Time for selected inversion = %15.2e\n', timeSelInv);

if(chkerr)
  indsymA = find(symA);
  Ainvfull = inv(A);
  fprintf('norm(Ainvfull-Ainv) = %15.5e\n', ...
    norm(Ainvfull(indsymA)-Ainv(indsymA)));
  fprintf('Tr(A*Ainv)-N=%15.5e\n', trace(Ainv*A)-numCol);
end
