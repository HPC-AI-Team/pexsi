N = 500;
A = sprandn(N,N,0.005)+speye(N)*2;
% AMD ordering
perm = symamd(A);
A = A(perm,perm);
fprintf('\nMirror left-looking...\n\n')
Ainv = SelInvMirrorLeft(A,1);
figure(1)
fprintf('\nSparsity of L,U factor...\n\n')
spy(Ainv)
pause
fprintf('\nMirror right-looking...\n\n')
Ainv = SelInvMirrorRight(A,1);
