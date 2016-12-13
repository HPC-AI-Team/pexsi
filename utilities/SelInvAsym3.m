function [Ainv,APreSelinv,Afactor] = SelInvAsym3( A, chkerr, b )
  % Prototype implementation of selected inversion for general asymmetric
  % matrices, without symmetrization and allowing row/column permutation.
  %
  % It also uses the transpose formulation and is the memory efficient
  % version and only uses Afactor.
  %
  % This is a column-by-column implementation and is very inefficient.
  %
  % It returns
  %
  % {A^{-1}_{ij} | A_{ji}\ne 0}
  %
  % Lin Lin
  % Revision: 11/18/2016

  numCol   = length( A );
  assert( size(A,1) == size(A,2) );
  if( nargin < 2 )
    chkerr = 0;
  end

  if( nargin < 3 )
    b = 1;
  end

  %% Symbolic LU factorization
  %  This is a very simple implementation.
  %  The resulting LU factor is structurally symmetric.
  %  However, in MATLAB this step does not seem to be very effective.
  tic;
  p = 1:size(A,1);
  q = 1:size(A,2);
  [L,U] = LU_factor(A);
  % MATLAB's LU factorization for A
  %[L,U,p,q] = lu(A,'vector');

  timeNumFac = toc;
  fprintf('Time for numerical factorization = %15.2e\n', timeNumFac);

  % Allocate the space for Afactor
  Afactor = tril(L,-1)+U;
  if(0)
    figure(1)
    spy(Afactor)
  end

  %% Selected inversion.
  tic;
  % Pre-selected inversion
%  superPtr = 0:numCol;
%  APreSelinv = PreSelInvUnsym(L,U,superPtr);
%  Afactor = APreSelinv;

  for k = 1 : numCol
    Afactor(k,k) = inv(Afactor(k,k));
    ind = k + find(Afactor(k,k+1:end));  
    if(~isempty(ind))
      Afactor(k,ind) = Afactor(k,ind) * Afactor(k,k);
    end
  end
  APreSelinv = Afactor;
 
  %norm(full(APreSelInv - Afactor))
  %pause
  %Afactor 
  %pause

  %Afactor=APreSelInv;

  for k = numCol : -1 : 1
    Afactor(k,k) = transpose(Afactor(k,k));
    indL = k + find(Afactor(k+1:end,k));
    indU = k + find(Afactor(k,k+1:end));  
    if(~isempty(indL) & ~isempty(indU))
      LBuf = transpose(Afactor(indL,k));
      UBuf = transpose(Afactor(k,indU));
      AinvBuf = Afactor(indL,indU); 
      Afactor(indL,k) = -Afactor(indL,indU) * UBuf;
      Afactor(k,indU) = -LBuf * Afactor(indL,indU);
      if( 0 && k==379)
        %AinvBuf


        LDiag = Afactor(k,k)
        LBuf = LBuf
        SinvU = Afactor(indL,k)
        DiagUpdtU = LBuf * Afactor(indL,k)

        UBuf = UBuf
        LSinv = Afactor(k,indU)
        DiagUpdtL = Afactor(k,indU)*UBuf



        Adiag = Afactor(k,k) - LBuf * Afactor(indL,k)
      end

      Afactor(k,k)   = Afactor(k,k)' - LBuf * Afactor(indL,k);
      %Afactor(k,k)   = Afactor(k,k) - Afactor(k,indU)*UBuf;
    end
  end

  timeSelInv = toc;
  fprintf('Time for selected inversion = %15.2e\n', timeSelInv);

  indsymAT = find(A.');
  % AinvT(p,q) = Afactor;
  invp = zeros(size(p));
  invp(p) = 1:length(p);
  invq = zeros(size(q));
  invq(q) = 1:length(q);
  AinvT = Afactor(invp,invq);
  % AinvT = P'*Afactor*Q';

  Ainv = transpose(AinvT);
  if(chkerr)
    Ainvfull = inv(A);
    fprintf('norm(Ainvfull-Ainv,sel) = %15.5e\n', ...
      norm(Ainvfull(indsymAT)-Ainv(indsymAT)));
  end

