function [zshift, zweight] = getpole( Npole, T, Gap, DeltaE, mu, func )
% getpole.m generates the poles and weights for the calculation of electronic
% density. 
% LL: VERY IMPORTANT: This getpole is very different from getpole.m in PARSEC/GetPole!
%
% Input:
%		Npole: the number of poles to be used.
%		T: temperature, unit(K)
%		Gap: 	Energy gap defined to be min(abs(EV-mu)). EV is the eigenvalue set
%			of Hamiltonian, unit(hartree)
%		DeltaE: Spectrum width defined to be max(EV)-min(EV). EV is the eigenvalue
%			set of Hamiltonian, unit(hartree)
%		mu: Chemical potential, unit(hartree)
%
% Output:
%		
%		zshift:  Shift of poles (polepos).
%		zweight: Weight of poles (poleweight).
% Usage of pole information to calculate eletronic density: 
%
%   Rho = zeros(N, 1);  
%   % VERY IMPORTANT!  This is different from previous version
%
%   for i = 1 : npoles
%     Calculate Result = diag( inv( (H - polepos(i))) )  !
%     % VERY IMPORTANT!  Sign of polepos is different from previous version
%     Rho = Rho + imag( poleweight(i) * Result );
%   end 
%

K2au = 3.166815d-6;
au2K = 315774.67;
beta = au2K/(T);

Npolehalf = Npole/2;
M = DeltaE;
mshift = (pi/beta)^2;
m2 = mshift+(Gap)^2;
M2 = M^2;
k = (sqrt(M2/m2)-1)/(sqrt(M2/m2)+1);
L = -log(k)/pi;
[K,Kp] = ellipkkp(L);

t = .5i*Kp - K + (.5:Npolehalf)*2*K/Npolehalf;  %t = .9999999999i*Kp - K + (.5:Npolehalf)*2*K/Npolehalf;
[u cn dn] = ellipjc(t,L);
z = sqrt(m2*M2)*((1/k+u)./(1/k-u));
dzdt = cn.*dn./(1/k-u).^2;

zsqrt = sqrt(z-mshift);


zweight = zeros(Npole, 1);
zshift = zeros(Npole, 1);

% From Eq. (2.10) in 
%   Lin, Lu, Ying and E, " Pole-Based Approximation of the Fermi-Dirac
%   Function",  Chin. Ann. Math. 30B, 729, 2009
%
for j = 1 : Npolehalf
  zshift(j) = mu + zsqrt(j);
  zshift(j+Npolehalf) = mu - zsqrt(j);
  % Old version which does not include the identity in the formulation.
  %
  % zweight(j) = ...
    % 2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / zsqrt(j) * ...
    % dzdt(j) * (-tanh(beta*zsqrt(j)/2));
  % zweight(j+Npolehalf) = ...
    % 2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / (-zsqrt(j)) * ...
    % dzdt(j) * (-tanh(beta*(-zsqrt(j))/2));
  %
  % New version which takes into account the identity.
  zweight(j) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / zsqrt(j) * ...
    dzdt(j) * func(zsqrt(j), 0, beta);
  zweight(j+Npolehalf) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / (-zsqrt(j)) * ...
    dzdt(j) * func(-zsqrt(j), 0, beta);
end
