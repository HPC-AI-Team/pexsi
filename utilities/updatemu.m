% updatemu computes the Fermi energy based on the same set of poles.
%
% Lin Lin
% Revision: 11/24/2015

rng(1);
eigval = 10*rand(10000,1);
eigval = sort(eigval);
T = 300;
K2au = 3.166815d-6;
beta = 1/(T*K2au);
numElecExact = 1000;

ResFunc = @(efermi) sum( fermidirac(eigval, efermi, beta) ) - ...
	numElecExact;
ResDrvFunc = @(efermi) sum( fermidiracdrv( eigval, efermi, beta) );

% Exact mu
[occ, muexact] = getocc(eigval, numElecExact, beta);


%% Pole expansion at exact mu, and compute traceG independent of mu.
Npole = 80;
Gap   = 0.0;
DeltaE = 10;
[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, muexact, 1);
traceG = zeros(Npole,1);
for l = 1 : Npole
  traceG(l) = sum(1.0./(eigval - zshift(l)));
end

NePole = imag(sum(zweight.*traceG));

fprintf('muExact     = %25.15e\n', muexact);
fprintf('NePole      = %25.15e\n', NePole);
fprintf('Error of Ne = %25.15e\n', abs(NePole-numElecExact));

%% Ne as a function of mu, comparing the exact solution and the pole
% expansion
muVec = linspace(muexact - T*K2au/2, muexact + T*K2au/2);
NeExactVec = zeros(size(muVec));
NePoleVec  = zeros(size(muVec));
for il = 1 : length(muVec)
  mu = muVec(il);
  dmu = mu - muexact;
  NeExactVec(il) = sum(fermidirac(eigval, mu, beta));
  [~, zweight] = getpole2(Npole, T, Gap, DeltaE, muexact, dmu, 1);
  NePoleVec(il) = imag(sum(zweight.*traceG));
end

figure
plot(muVec, NeExactVec, 'b-o', muVec, NePoleVec, 'r-d');
legend('Exact','Pole');
figure
plot(muVec,NeExactVec-NePoleVec,'k-p');

%% Pole expansion using the same traceG but changing weights
%
% Bisection for finding mu.
mu0 = muexact - 0.5*T*K2au;
[zshift0, zweight] = getpole(Npole, T, Gap, DeltaE, mu0, 1);
traceG0 = zeros(Npole,1);
for l = 1 : Npole
  traceG0(l) = sum(1.0./(eigval - zshift0(l)));
end
fprintf('mu0         = %25.15e\n', mu0);
fprintf('Ne0         = %25.15e\n', imag(sum(zweight.*traceG0)));

muwidth = 5*T*K2au;
mumin = mu0 - muwidth;
mumax = mu0 + muwidth;

mu = mu0;
tol = 1e-3;
for il = 1 : 10
  dmu = mu - mu0;
  [~, zweight] = getpole2(Npole, T, Gap, DeltaE, mu0, dmu, 1);
  Ne = imag(sum(zweight.*traceG0));
  if( abs(Ne - numElecExact) < tol )
    break;
  end
  if( Ne < numElecExact )
    mumin = mu;
    mu = 0.5*(mu+mumax);
  else
    mumax = mu;
    mu = 0.5*(mumin+mu);
  end
end

% Compute the number of electrons using pole expansion at mu
[zshiftmu, zweight] = getpole(Npole, T, Gap, DeltaE, mu, 1);
traceGmu = zeros(Npole,1);
for l = 1 : Npole
  traceGmu(l) = sum(1.0./(eigval - zshift(l)));
end
NePoleMu = imag(sum(zweight.*traceGmu));

fprintf('muFind      = %25.15e\n', mu);
fprintf('NeFind      = %25.15e\n', Ne);
fprintf('NePoleMu    = %25.15e\n', NePoleMu);

