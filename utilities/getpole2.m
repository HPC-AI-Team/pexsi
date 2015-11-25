function [zshift, zweight] = getpole2( Npole, T, Gap, DeltaE, mu, dmu, nspin )
% getpole2.m obtains the shift evaluated at mu, and weights evaluated at
% mu+dmu.  
%
% Note: this is only accurate when dmu is small.
%
% Lin Lin
% Revision: 11/24/2015

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

fd = @(z) nspin./(1+exp(beta*z));  % NOTE: No spin degeneracy here!

% From Eq. (2.10) in 
%   Lin, Lu, Ying and E, " Pole-Based Approximation of the Fermi-Dirac
%   Function",  Chin. Ann. Math. 30B, 729, 2009
%
for j = 1 : Npolehalf
  zshift(j) = mu + zsqrt(j);
  zshift(j+Npolehalf) = mu - zsqrt(j);
  zweight(j) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / zsqrt(j) * ...
    dzdt(j) * fd(zsqrt(j)-dmu);
  zweight(j+Npolehalf) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / (-zsqrt(j)) * ...
    dzdt(j) * fd(-zsqrt(j)-dmu);
end
