function f = fermidirac(ev,efermi,Tbeta);
% Fermi-Dirac function
dif = ev - efermi;
if( numel(efermi) > 1 )
	error('efermi has to be a scalar');
end
f = zeros(size(ev));
ind1 = find(dif > 0 );
ind2 = find(dif <= 0 );
f(ind1) = exp(-Tbeta * dif(ind1))./(1+exp(-Tbeta*dif(ind1)));
f(ind2) = 1./(1+exp(Tbeta*dif(ind2)));
