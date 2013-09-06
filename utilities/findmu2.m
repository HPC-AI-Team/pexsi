% findmu2 finds the chemical potential for a range of T.
%
% Lin Lin
% 12/28/2012

eigval = load('eigval.dat');
K2au = 3.166815d-6;
T = K2au * [10:10:800 1000 3000 5000 10000];
T = T(end:-1:1);
numElecExact = 1221;
% numElecExact = 1272;

muT = zeros(size(T));
muZeroT = zeros(size(T));
for i = 1 : length(T)
	beta = 1 / T(i);
	ResFunc = @(efermi) sum( fermidirac(eigval, efermi, beta) ) - ...
		numElecExact;
	ResDrvFunc = @(efermi) sum( fermidiracDrvMu( eigval, efermi, beta) );
	ResDrvTFunc = @(efermi) sum( fermidiracDrvT( eigval, efermi, beta) );

	if( i == 1 )
		mu0 = -0.65;
	else
		mu0 = muZeroT(i-1);
	end
	[muT(i), resval, exitflag, output] = fzero( ResFunc, mu0 );
	dMudT = -ResDrvTFunc(muT(i)) / ResDrvFunc(muT(i));
	muZeroT(i) = muT(i) - T(i) / 2 * dMudT;
	fprintf('At T = %10.3f K, mu(T) = %15.5f, mu(T=0) = %15.5f\n', ...
		T(i)/K2au, muT(i), muZeroT(i));
end


% [muHist1, fHist1, gHist1] = findzero1( ResFunc, ResDrvFunc, ...
	% mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	% 0, 0, 30, 1e-3);
% 
% [muHist2, fHist2, gHist2] = findzero2( ResFunc, ResDrvFunc, ...
	% mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	% 0, 0, 30, 1e-3);
% 
% [muHist3, fHist3, gHist3] = findzero3( ResFunc, ResDrvFunc, ...
	% mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	% 0, 0, 30, 1e-3);
% 
% [muHist4, fHist4, gHist4] = findzero4( ResFunc, ResDrvFunc, ...
	% mu0, -5, 15, 0-numElecExact, length(eigval)-numElecExact, ...
	% 0, 0, 30, 1e-3);
