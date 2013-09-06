% findmu finds the chemical potential using root-finding method.
%
% Lin Lin
% 12/28/2012

eigval = load('eigval.dat');
T = 3000;
K2au = 3.166815d-6;
beta = 1/(T*K2au);
numElecExact = 1221;
% numElecExact = 1272;

ResFunc = @(efermi) sum( fermidirac(eigval, efermi, beta) ) - ...
	numElecExact;
ResDrvFunc = @(efermi) sum( fermidiracdrv( eigval, efermi, beta) );

mu0 = -0.65478; 
mumin = -3;
mumax =  15;
[mu, resval, exitflag, output] = fzero( ResFunc, [mumin, mumax] );
output

[muHist1, fHist1, gHist1] = findzero1( ResFunc, ResDrvFunc, ...
	mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	0, 0, 30, 1e-3);

[muHist2, fHist2, gHist2] = findzero2( ResFunc, ResDrvFunc, ...
	mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	0, 0, 30, 1e-3);

[muHist3, fHist3, gHist3] = findzero3( ResFunc, ResDrvFunc, ...
	mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	0, 0, 30, 1e-3);

[muHist4, fHist4, gHist4] = findzero4( ResFunc, ResDrvFunc, ...
	mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	0, 0, 30, 1e-3);
