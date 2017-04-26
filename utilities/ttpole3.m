% TTPOLE is an example to test that the pole expansion works for Fermi-Dirac
% distribution.
%
% This tests the electron density for complex Hermitian matrices
%
% Lin Lin
% Lastest revision: 09/09/2016

% The unit of energy is hartree!
Npole    = 120;
T        = 300;
Gap      = 0.0;
DeltaE   = 2.0;
mu       = 0.270;
nSpin    = 2;

K2au = 3.166815d-6;
beta = 1/(T*K2au);

disp('Fermi-Dirac');
lap2dcgen;
H = A;
Ns = length(H);
S = speye(Ns);
[Psi,Ev]=eig(full(H));
Ev = diag(Ev);

[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu, nSpin);

occExact = nSpin*fermidirac(Ev, mu, beta);

RhoExact = Psi*diag(occExact)*Psi';

RhoPole = zeros(size(RhoExact));

for l = 1 : Npole
  Gl = inv(full(H - zshift(l)*S));
  % Does not work for complex Hermitian matrices.
  % The number of electrons is correct, but the electron density is not
  % RhoPole = RhoPole + imag(zweight(l)*Gl);
  
  % Correct version: sum and correct
  RhoPole = RhoPole + zweight(l)*Gl;
end
% Correction step
RhoPole = 1/(2i)*(RhoPole-RhoPole');

fprintf('Number of electrons (exact) = %25.15f\n', sum(diag(RhoExact)) );
fprintf('Number of electrons (pole)  = %25.15f\n', sum(diag(RhoPole)) );
fprintf('||RhoExact - RhoPole||_2    = %25.15f\n', ...
  norm(RhoExact-RhoPole));
