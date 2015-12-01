function Ainv = SelInvMirrorRight( A, chkerr )
% Prototype implementation of the "Mirror Right-looking" version of
% selected inversion. 
%
%   Ainv = SelInvMirrorRight( A, chkerr )
%
% Assume A is pre-ordered, and the L, U factors exist. It first computes
% the column-wise factorization of A, and then compute selected
% inversion restricted to the pattern of (L+U)+(L+U)'.
%
% The "Mirror right-looking", i.e. left-looking selected inversion
% algorithm is implemented. The factorization is right-looking.
%
% This is a column-by-column implementation and is very inefficient.
%
% Lin Lin
% Revision: 11/24/2015

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
symA = symA + eps*speye(numCol);
for k = 1 : numCol
  ind = k + find(symA(k+1:end,k));
  symA(ind, ind) = eps;
end
symA(find(symA)) = eps;

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
  ind = k+find(symA(k+1:end,k));
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

%% Selected inversion. Mirror right-looking
tic;
Ainv = symA;
% Pre-selected inversion
for k = 1 : numCol
  Afactor(k,k) = inv(Afactor(k,k));
	ind = k + find(symA(k+1:end,k));  
  if(~isempty(ind))
    Afactor(k,ind) = Afactor(k,ind) * Afactor(k,k);
  end
end

for k = numCol : -1 : 1
  % SelInv: the only work left is the diagonal
  Ainv(k,k) = Afactor(k,k);
	ind = k + find(symA(k+1:end,k));  
  ind2 = find(symA(k,1:k-1));
  Ainv(k,k) = Ainv(k,k) - Afactor(k,ind) * Ainv(ind,k);
  
  % Update: Mirror right-looking. L part
  if(~isempty(ind2))
    % Diagonal
    Ainv(k,ind2) = Ainv(k,ind2) - Ainv(k,k) * Afactor(k,ind2);
  end
  if(~isempty(ind) & ~isempty(ind2))
    for jp = ind2
      indi = ind(find(Afactor(ind,jp)));
      % Outer product 
      Ainv(indi,jp) = Ainv(indi,jp) - Ainv(indi,k)*Afactor(k,jp);
      % Inner product
      Ainv(k,jp) = Ainv(k,jp) - Ainv(k,indi) * Afactor(indi,jp);
    end
  end

  % Update: Mirror right-looking. U part
  if(~isempty(ind2))
    % Diagonal
    Ainv(ind2,k) = Ainv(ind2,k) - Afactor(ind2,k) * Ainv(k,k);
  end
  if(~isempty(ind) & ~isempty(ind2))
    for jp = ind2
      indi = ind(find(Afactor(ind,jp)));
      % Outer product 
      Ainv(jp,indi) = Ainv(jp,indi) - Afactor(jp,k)*Ainv(k,indi);
      % Inner product
      Ainv(jp,k) = Ainv(jp,k) - Afactor(jp,indi)*Ainv(indi,k);
    end

    % % Outer product 
    % Ainv(ind2,ind) = Ainv(ind2,ind) - ...
      % Afactor(ind2,k)*Ainv(k,ind);
    % % Inner product
    % Ainv(ind2,k) = Ainv(ind2,k) - Afactor(ind2,ind)*Ainv(ind,k);
  end
end
timeSelInv = toc;
fprintf('nnz(Ainv) = %g, nnz(Afactor) = %g\n', nnz(Ainv), nnz(Afactor))
fprintf('Time for selected inversion = %15.2e\n', timeSelInv);

if(chkerr)
  indsymA = find(symA);
  Ainvfull = inv(A);
  fprintf('norm(Ainvfull-Ainv) = %15.5e\n', ...
    norm(Ainvfull(indsymA)-Ainv(indsymA)));
  fprintf('Tr(A*Ainv)-N=%15.5e\n', trace(Ainv*A)-numCol);
end
