% TTPOLE is an example to test that the pole expansion works for Fermi-Dirac
% distribution.
%
% Lin Lin
% Lastest revision: 09/11/2011

% The unit of energy is hartree!
Npole    = 80;
T        = 300;
Gap      = 0.0;
DeltaE   = 2.0;
% mu       = -6.550428925333580e-01;
mu       = -1.478040987387258e+00;
% mu       = -1.481346204776733e+00;

K2au = 3.166815d-6;
beta = 1/(T*K2au);

disp('Fermi-Dirac');

[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu, @fermidirac);
x = eigval;
f1 = fermidirac(x, mu, beta);
f2 = zeros(size(x));
for l = 1 : Npole
  f2 = f2 + imag(zweight(l)./(x-zshift(l)));
end

fprintf('Number of electrons (exact) = %25.15f\n', sum(f1) * 2.0 );
fprintf('Number of electrons (pole)  = %25.15f\n', sum(f2) * 2.0 );


